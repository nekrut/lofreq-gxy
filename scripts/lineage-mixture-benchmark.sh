#!/usr/bin/env bash
# Tier 2.3 — lineage-mixture benchmark with absolute truth.
#
# Injects the Alpha (B.1.1.7) lineage-signature SNV set into the
# SARS-CoV-2 Wuhan reference at a range of mixture fractions (AFs)
# and measures precision / recall on the resulting calls against
# the known truth. Unlike parity-against-upstream, this produces
# first-class precision/recall at each injected AF bin — the
# metric the E. coli result showed we need.
#
# Usage:
#   scripts/lineage-mixture-benchmark.sh [depth]
#
# Default depth 1000× (enough to see ~1 % with confidence).
#
# Requires: target/release/gxy-make-fixture + lofreq-gxy-test conda env
# + parity/upstream/bin/lofreq.
#
# Alpha signature set (18 SNVs in NC_045512.2 coordinates):
#   nsp3:         C3267T C5388A T6954C
#   nsp12:        C14676T C15279T T16176C
#   S:            A23063T C23271A C23604A C23709T T24506G G24914C
#   ORF8:         C27972T G28048T A28111G
#   N (GAT→ACA):  G28881A G28882A G28883C

set -euo pipefail

DEPTH="${1:-1000}"
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
GXY="$ROOT/target/release/lofreq-gxy"
UPSTREAM="$ROOT/parity/upstream/bin/lofreq"
FIXTURE="$ROOT/target/release/gxy-make-fixture"
REF="$ROOT/parity/fixtures/sars_cov_2.fa"
OUTROOT="$ROOT/parity/compare/lineage_mixture"

mkdir -p "$OUTROOT"
for tool in samtools bcftools bgzip tabix; do
  command -v "$tool" >/dev/null || { echo "error: $tool not on PATH" >&2; exit 1; }
done
[[ -x "$GXY" ]] || { echo "error: $GXY missing"; exit 1; }
[[ -x "$UPSTREAM" ]] || { echo "error: $UPSTREAM missing"; exit 1; }
[[ -x "$FIXTURE" ]] || { echo "error: $FIXTURE missing"; exit 1; }

# Alpha signature SNVs: pos:alt (reference base implied by the ref).
ALPHA_POS=(3267 5388 6954 14676 15279 16176 23063 23271 23604 23709 24506 24914 27972 28048 28111 28881 28882 28883)
ALPHA_ALT=(T    A    C    T     T     C     T     A     A     T     G     C     T     T     G     A     A     C)
N_SIGS=${#ALPHA_POS[@]}

AFS=(0.01 0.025 0.05 0.10 0.25 0.50)

echo "Alpha signature SNV count: $N_SIGS"
echo "Depths per mixture: $DEPTH"
echo "Target AFs: ${AFS[*]}"
echo

SUMMARY="$OUTROOT/summary.tsv"
printf "af\ttool\ttp\tfp\tfn\tprecision\trecall\tf1\tmean_delta_af\n" > "$SUMMARY"

for af in "${AFS[@]}"; do
  mix_dir="$OUTROOT/af_${af}"
  mkdir -p "$mix_dir"
  bam="$mix_dir/mix.sorted.bam"

  # 1. Generate the mixture BAM.
  if [[ ! -s "$bam" ]]; then
    echo "[af=$af] generating mixture at depth $DEPTH"
    variant_args=()
    for i in $(seq 0 $((N_SIGS-1))); do
      variant_args+=(--variant "${ALPHA_POS[$i]}:${ALPHA_ALT[$i]}:$af")
    done
    "$FIXTURE" --ref "$REF" --out "$mix_dir/mix.bam" \
      --depth "$DEPTH" --seed 1 "${variant_args[@]}" >/dev/null
    samtools sort "$mix_dir/mix.bam" -o "$bam"
    samtools index "$bam"
    rm "$mix_dir/mix.bam"
  fi

  # 2. Build the truth VCF for this AF (one record per signature position).
  truth="$mix_dir/truth.vcf"
  {
    echo "##fileformat=VCFv4.2"
    echo "##INFO=<ID=AF,Number=1,Type=Float,Description=\"Injected AF\">"
    awk '/^>/{next} {seq=seq $0} END{print seq}' "$REF" > "$mix_dir/.seq"
    echo "##contig=<ID=NC_045512.2,length=$(wc -c < "$mix_dir/.seq" | awk '{print $1-1}')>"
    echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
    for i in $(seq 0 $((N_SIGS-1))); do
      pos=${ALPHA_POS[$i]}
      alt=${ALPHA_ALT[$i]}
      ref_base=$(awk -v p="$pos" 'NR==1{next} {s=s $0} END{print toupper(substr(s,p,1))}' "$REF")
      echo -e "NC_045512.2\t${pos}\t.\t${ref_base}\t${alt}\t.\tPASS\tAF=${af}"
    done
  } > "$truth"
  rm -f "$mix_dir/.seq"
  bcftools sort -O v "$truth" 2>/dev/null | bgzip -f > "$truth.gz"
  tabix -f "$truth.gz"

  # 3. Run both callers.
  for tool in gxy upstream; do
    vcf="$mix_dir/$tool.vcf"
    if [[ "$tool" == "gxy" ]]; then
      "$GXY" call -f "$REF" -o "$vcf" --force-overwrite "$bam" >/dev/null 2>&1
    else
      "$UPSTREAM" call -f "$REF" -o "$vcf" --force "$bam" >/dev/null 2>&1
    fi
    # Prepare for isec.
    bcftools reheader --fai "${REF}.fai" "$vcf" \
      | bcftools sort -O v 2>/dev/null \
      | bgzip -f > "$vcf.gz"
    tabix -f "$vcf.gz"

    # 4. Compute TP/FP/FN and mean |ΔAF| on TPs.
    isec_dir="$mix_dir/${tool}_isec"
    rm -rf "$isec_dir"
    bcftools isec -p "$isec_dir" "$vcf.gz" "$truth.gz" >/dev/null 2>&1
    fp=$(grep -cv '^#' "$isec_dir/0000.vcf" || true)
    fn=$(grep -cv '^#' "$isec_dir/0001.vcf" || true)
    tp=$(grep -cv '^#' "$isec_dir/0002.vcf" || true)
    if (( tp + fp > 0 )); then
      prec=$(awk -v t="$tp" -v f="$fp" 'BEGIN{printf "%.4f", t/(t+f)}')
    else
      prec="n/a"
    fi
    if (( tp + fn > 0 )); then
      rec=$(awk -v t="$tp" -v f="$fn" 'BEGIN{printf "%.4f", t/(t+f)}')
    else
      rec="n/a"
    fi
    if [[ "$prec" != "n/a" && "$rec" != "n/a" ]]; then
      f1=$(awk -v p="$prec" -v r="$rec" 'BEGIN{if(p+r>0) printf "%.4f", 2*p*r/(p+r); else print "n/a"}')
    else
      f1="n/a"
    fi
    # mean |ΔAF| on TP calls (compare called AF to target AF).
    mean_daf=$(bcftools query -f '%INFO/AF\n' "$isec_dir/0002.vcf" 2>/dev/null \
      | awk -v tgt="$af" 'BEGIN{s=0; n=0} {d=$1-tgt; if(d<0)d=-d; s+=d; n++} END{if(n) printf "%.4f", s/n; else print "n/a"}')
    printf "%s\t%s\t%d\t%d\t%d\t%s\t%s\t%s\t%s\n" "$af" "$tool" "$tp" "$fp" "$fn" "$prec" "$rec" "$f1" "$mean_daf" >> "$SUMMARY"
  done
done

echo
echo "=== Summary ==="
column -t -s $'\t' "$SUMMARY"
echo
echo "Wrote $SUMMARY"
