#!/usr/bin/env bash
# Run lofreq-gxy and the pinned upstream lofreq on the same BAM + FASTA
# pair, then compute concordance stats with bcftools isec.
#
# Usage: scripts/compare.sh <bam> <fasta> [output-dir]
#
# Requires:
#   - target/release/lofreq-gxy
#   - parity/upstream/bin/lofreq
#   - bcftools, bgzip, samtools on PATH
#
# Output directory gets:
#   gxy.vcf             lofreq-gxy output
#   upstream.vcf        upstream lofreq output
#   gxy.vcf.gz + .tbi   bgzipped+indexed
#   upstream.vcf.gz + .tbi
#   isec/               bcftools isec tree
#   report.txt          summary stats (timing, jaccard, AF diff)

set -euo pipefail

BAM="${1:?usage: compare.sh <bam> <fasta> [out-dir]}"
FASTA="${2:?usage: compare.sh <bam> <fasta> [out-dir]}"
OUT="${3:-parity/compare}"

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
GXY="$ROOT/target/release/lofreq-gxy"
UPSTREAM="$ROOT/parity/upstream/bin/lofreq"

if [[ ! -x "$GXY" ]]; then
  echo "error: $GXY missing. Run 'cargo build --release' first." >&2
  exit 1
fi
if [[ ! -x "$UPSTREAM" ]]; then
  echo "error: $UPSTREAM missing. Run 'scripts/build-upstream.sh' first." >&2
  exit 1
fi

mkdir -p "$OUT"

# --- Run both tools, timing wall-clock wait + user CPU --------------------
echo "[compare] running lofreq-gxy …"
/usr/bin/time -f 'gxy_real=%e gxy_user=%U gxy_sys=%S gxy_maxrss=%M' \
  "$GXY" call -f "$FASTA" -o "$OUT/gxy.vcf" --force-overwrite "$BAM" \
  2> "$OUT/gxy.time.log" >/dev/null

echo "[compare] running upstream lofreq …"
/usr/bin/time -f 'upstream_real=%e upstream_user=%U upstream_sys=%S upstream_maxrss=%M' \
  "$UPSTREAM" call -f "$FASTA" -o "$OUT/upstream.vcf" --force "$BAM" \
  2> "$OUT/upstream.time.log" >/dev/null

# --- Normalise + compress for bcftools isec -------------------------------
# Upstream lofreq omits `##contig` lines; reheader against the FASTA index
# so bcftools sort / isec can parse. `.fai` is generated on demand.
if [[ ! -f "${FASTA}.fai" ]]; then
  samtools faidx "$FASTA"
fi
for name in gxy upstream; do
  bcftools reheader --fai "${FASTA}.fai" "$OUT/$name.vcf" > "$OUT/$name.reheadered.vcf"
  bcftools sort -O v -o "$OUT/$name.sorted.vcf" "$OUT/$name.reheadered.vcf" 2>/dev/null
  bgzip -f "$OUT/$name.sorted.vcf"
  tabix -p vcf -f "$OUT/$name.sorted.vcf.gz"
done

# --- bcftools isec ---------------------------------------------------------
rm -rf "$OUT/isec"
bcftools isec -p "$OUT/isec" \
  "$OUT/gxy.sorted.vcf.gz" "$OUT/upstream.sorted.vcf.gz" >/dev/null

only_gxy=$(grep -cv '^#' "$OUT/isec/0000.vcf" || true)
only_upstream=$(grep -cv '^#' "$OUT/isec/0001.vcf" || true)
shared=$(grep -cv '^#' "$OUT/isec/0002.vcf" || true)

union=$((only_gxy + only_upstream + shared))
jaccard="n/a"
if (( union > 0 )); then
  jaccard=$(awk -v s="$shared" -v u="$union" 'BEGIN{printf "%.4f", s/u}')
fi

# --- Report ----------------------------------------------------------------
{
  echo "=== lofreq-gxy vs upstream lofreq ==="
  echo "bam:      $BAM"
  echo "fasta:    $FASTA"
  echo
  echo "--- Timing ---"
  cat "$OUT/gxy.time.log"
  cat "$OUT/upstream.time.log"
  echo
  echo "--- Concordance (all calls) ---"
  echo "only in gxy:       $only_gxy"
  echo "only in upstream:  $only_upstream"
  echo "shared:            $shared"
  echo "jaccard:           $jaccard"
  echo
  echo "--- AF-level diff on shared calls (|delta AF|) ---"
  paste <(bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' "$OUT/gxy.sorted.vcf.gz") \
        <(bcftools query -f '%INFO/AF\n' "$OUT/upstream.sorted.vcf.gz") \
    2>/dev/null | \
    awk -F'\t' 'NF==6 { d=$5-$6; if(d<0) d=-d; if(d>max) max=d; sum+=d; n++ }
                END { if(n) printf "max=%.2e  mean=%.2e  (n=%d)\n", max, sum/n, n; else print "n=0" }'
} | tee "$OUT/report.txt"
