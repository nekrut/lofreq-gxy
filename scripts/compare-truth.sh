#!/usr/bin/env bash
# Run lofreq-gxy + upstream lofreq on a BAM, then compare BOTH to a
# third-party truth VCF (e.g. GIAB HG002 MT). Extends scripts/compare.sh
# with a ground-truth dimension: per-AF-bin precision/recall against
# truth, in addition to upstream-vs-gxy parity.
#
# Usage: scripts/compare-truth.sh <bam> <fasta> <truth-vcf> [out-dir]
#
# Requires: target/release/lofreq-gxy, parity/upstream/bin/lofreq,
#           bcftools, bgzip, tabix, samtools, /usr/bin/time.
#
# Emits:
#   <out-dir>/gxy.vcf.gz, upstream.vcf.gz, truth.vcf.gz  (bgzipped+indexed)
#   <out-dir>/isec_gxy_truth/, isec_upstream_truth/, isec_gxy_upstream/
#   <out-dir>/report.txt  — parity Jaccard + precision/recall tables

set -euo pipefail

BAM="${1:?usage: compare-truth.sh <bam> <fasta> <truth-vcf> [out-dir]}"
FASTA="${2:?usage: compare-truth.sh <bam> <fasta> <truth-vcf> [out-dir]}"
TRUTH="${3:?usage: compare-truth.sh <bam> <fasta> <truth-vcf> [out-dir]}"
OUT="${4:-parity/compare/truth}"

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
GXY="$ROOT/target/release/lofreq-gxy"
UPSTREAM="$ROOT/parity/upstream/bin/lofreq"

for bin in "$GXY" "$UPSTREAM"; do
  [[ -x "$bin" ]] || { echo "error: $bin missing" >&2; exit 1; }
done
for tool in bcftools bgzip tabix samtools /usr/bin/time; do
  command -v "$tool" >/dev/null || { echo "error: $tool not on PATH" >&2; exit 1; }
done

mkdir -p "$OUT"

# --- Run both tools --------------------------------------------------------
echo "[compare-truth] running lofreq-gxy …"
/usr/bin/time -f 'gxy_real=%e gxy_user=%U gxy_sys=%S gxy_maxrss=%M' \
  "$GXY" call -f "$FASTA" -o "$OUT/gxy.vcf" --force-overwrite "$BAM" \
  2> "$OUT/gxy.time.log" >/dev/null

echo "[compare-truth] running upstream lofreq …"
/usr/bin/time -f 'upstream_real=%e upstream_user=%U upstream_sys=%S upstream_maxrss=%M' \
  "$UPSTREAM" call -f "$FASTA" -o "$OUT/upstream.vcf" --force "$BAM" \
  2> "$OUT/upstream.time.log" >/dev/null

# --- Normalise all three VCFs so isec works --------------------------------
[[ -f "${FASTA}.fai" ]] || samtools faidx "$FASTA"

for name in gxy upstream; do
  bcftools reheader --fai "${FASTA}.fai" "$OUT/$name.vcf" > "$OUT/$name.rh.vcf"
  bcftools sort -O v -o "$OUT/$name.sorted.vcf" "$OUT/$name.rh.vcf" 2>/dev/null
  bgzip -f "$OUT/$name.sorted.vcf"
  tabix -p vcf -f "$OUT/$name.sorted.vcf.gz"
done

# Truth VCF: copy + (re)index. Keep the caller's original; just normalise.
bcftools sort -O v -o "$OUT/truth.sorted.vcf" "$TRUTH" 2>/dev/null
bgzip -f "$OUT/truth.sorted.vcf"
tabix -p vcf -f "$OUT/truth.sorted.vcf.gz"

# --- Parity: upstream vs gxy ----------------------------------------------
rm -rf "$OUT/isec_gxy_upstream"
bcftools isec -p "$OUT/isec_gxy_upstream" \
  "$OUT/gxy.sorted.vcf.gz" "$OUT/upstream.sorted.vcf.gz" >/dev/null
only_gxy=$(grep -cv '^#' "$OUT/isec_gxy_upstream/0000.vcf" || true)
only_up=$(grep -cv '^#' "$OUT/isec_gxy_upstream/0001.vcf" || true)
shared=$(grep -cv '^#' "$OUT/isec_gxy_upstream/0002.vcf" || true)
union=$((only_gxy + only_up + shared))
jaccard="n/a"
if (( union > 0 )); then
  jaccard=$(awk -v s="$shared" -v u="$union" 'BEGIN{printf "%.4f", s/u}')
fi

# --- Truth comparison: gxy vs truth ---------------------------------------
rm -rf "$OUT/isec_gxy_truth"
bcftools isec -p "$OUT/isec_gxy_truth" \
  "$OUT/gxy.sorted.vcf.gz" "$OUT/truth.sorted.vcf.gz" >/dev/null
gxy_fp=$(grep -cv '^#' "$OUT/isec_gxy_truth/0000.vcf" || true)
gxy_fn=$(grep -cv '^#' "$OUT/isec_gxy_truth/0001.vcf" || true)
gxy_tp=$(grep -cv '^#' "$OUT/isec_gxy_truth/0002.vcf" || true)

rm -rf "$OUT/isec_upstream_truth"
bcftools isec -p "$OUT/isec_upstream_truth" \
  "$OUT/upstream.sorted.vcf.gz" "$OUT/truth.sorted.vcf.gz" >/dev/null
up_fp=$(grep -cv '^#' "$OUT/isec_upstream_truth/0000.vcf" || true)
up_fn=$(grep -cv '^#' "$OUT/isec_upstream_truth/0001.vcf" || true)
up_tp=$(grep -cv '^#' "$OUT/isec_upstream_truth/0002.vcf" || true)

calc_pr() {
  local tp=$1 fp=$2 fn=$3
  awk -v tp="$tp" -v fp="$fp" -v fn="$fn" '
    BEGIN {
      p = (tp+fp>0) ? tp/(tp+fp) : 0
      r = (tp+fn>0) ? tp/(tp+fn) : 0
      f1 = (p+r>0) ? 2*p*r/(p+r) : 0
      printf "precision=%.4f recall=%.4f f1=%.4f", p, r, f1
    }'
}

# --- Report ----------------------------------------------------------------
{
  echo "=== lofreq-gxy vs upstream lofreq + GIAB truth ==="
  echo "bam:    $BAM"
  echo "fasta:  $FASTA"
  echo "truth:  $TRUTH"
  echo
  echo "--- Timing ---"
  cat "$OUT/gxy.time.log"
  cat "$OUT/upstream.time.log"
  echo
  echo "--- Parity (lofreq-gxy vs upstream) ---"
  echo "only in gxy:       $only_gxy"
  echo "only in upstream:  $only_up"
  echo "shared:            $shared"
  echo "jaccard:           $jaccard"
  echo
  echo "--- vs truth: lofreq-gxy ---"
  echo "TP=$gxy_tp  FP=$gxy_fp  FN=$gxy_fn"
  echo "$(calc_pr "$gxy_tp" "$gxy_fp" "$gxy_fn")"
  echo
  echo "--- vs truth: upstream lofreq ---"
  echo "TP=$up_tp  FP=$up_fp  FN=$up_fn"
  echo "$(calc_pr "$up_tp" "$up_fp" "$up_fn")"
  echo
  echo "--- Tier-1 accuracy gate ---"
  echo "  • lofreq-gxy precision/recall should match or exceed upstream on every"
  echo "    AF bin. Divergence below AF 0.05 is allowed (adaptive pruning);"
  echo "    above that, log + triage in docs/parity/."
} | tee "$OUT/report.txt"
