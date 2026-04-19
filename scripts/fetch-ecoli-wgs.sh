#!/usr/bin/env bash
# Fetch + subsample + align an E. coli K-12 MG1655 WGS Illumina
# dataset for Tier 2.4 parity testing. Downsamples to an approximate
# target coverage (default 30×) via `seqtk sample`, which is
# deterministic given a seed.
#
# Usage:
#   scripts/fetch-ecoli-wgs.sh [target_depth]
#
# Data source (per TEST-PLAN §2.4):
#   - ENA ERR022075 — E. coli K-12 MG1655, Illumina GAIIx paired-end,
#     ~22.7M pairs = ~990× raw. We subsample to the requested depth.
#   - Reference: NCBI NC_000913.3 (U00096.3) — E. coli K-12 MG1655.
#
# Produces, under parity/fixtures/ecoli_wgs/:
#   reads_R1.fastq.gz / reads_R2.fastq.gz     subsampled pairs
#   NC_000913.3.fa{.fai,.bwt,...}              reference + bwa index
#   ecoli_wgs.sorted.bam (+ .bai)              aligned sorted BAM
#
# Requires: bwa, samtools, seqtk, curl — all via lofreq-gxy-test conda env.
# Idempotent: re-running skips already-fetched / already-built files.

set -euo pipefail

TARGET_DEPTH="${1:-30}"
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
OUT="$ROOT/parity/fixtures/ecoli_wgs"
REF_URL="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_000913.3&rettype=fasta&retmode=text"
R1_URL="https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR022/ERR022075/ERR022075_1.fastq.gz"
R2_URL="https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR022/ERR022075/ERR022075_2.fastq.gz"

# Target coverage arithmetic. E. coli K-12 genome ≈ 4,641,652 bp.
# Raw dataset is ~4.59 Gb total, so ~990× coverage. Subsample fraction
# scales linearly: 30× / 990× ≈ 0.030.
GENOME_BP=4641652
RAW_BASES=4589460200
FRACTION=$(awk -v td="$TARGET_DEPTH" -v bp="$GENOME_BP" -v rb="$RAW_BASES" \
  'BEGIN{ printf "%.4f", (td*bp)/rb }')

mkdir -p "$OUT"
cd "$OUT"

for tool in bwa samtools seqtk curl; do
  command -v "$tool" >/dev/null || { echo "error: $tool not on PATH" >&2; exit 1; }
done

echo "[ecoli-wgs] target dir: $OUT"
echo "[ecoli-wgs] target depth: ${TARGET_DEPTH}× → subsample fraction $FRACTION"

# --- reference --------------------------------------------------------
if [[ ! -s NC_000913.3.fa ]]; then
  echo "[ecoli-wgs] fetching E. coli K-12 MG1655 reference (NC_000913.3)"
  curl -fSL --max-time 60 "$REF_URL" -o NC_000913.3.fa
  # Normalise the header to a short single token — bwa & samtools prefer it
  # short and stable across retrievals.
  awk 'NR==1{print ">NC_000913.3"; next}{print}' NC_000913.3.fa > NC_000913.3.fa.tmp \
    && mv NC_000913.3.fa.tmp NC_000913.3.fa
fi
[[ -s NC_000913.3.fa.fai ]] || samtools faidx NC_000913.3.fa
[[ -s NC_000913.3.fa.bwt ]] || bwa index NC_000913.3.fa

# --- raw FASTQ (full) -------------------------------------------------
for side in 1 2; do
  url_var="R${side}_URL"
  url="${!url_var}"
  dest="raw_R${side}.fastq.gz"
  if [[ ! -s "$dest" ]]; then
    echo "[ecoli-wgs] fetching R${side}: $url (~2 GB)"
    curl -fSL --max-time 1800 "$url" -o "$dest"
  fi
done

# --- subsample with deterministic seed --------------------------------
SEED=1
if [[ ! -s reads_R1.fastq.gz ]]; then
  echo "[ecoli-wgs] seqtk sample -s $SEED fraction=$FRACTION"
  seqtk sample -s "$SEED" raw_R1.fastq.gz "$FRACTION" | gzip > reads_R1.fastq.gz
  seqtk sample -s "$SEED" raw_R2.fastq.gz "$FRACTION" | gzip > reads_R2.fastq.gz
fi

# --- align ------------------------------------------------------------
BAM="ecoli_wgs.sorted.bam"
if [[ ! -s "$BAM" ]]; then
  echo "[ecoli-wgs] bwa mem → samtools sort"
  bwa mem -t "$(nproc)" -R "@RG\tID:ecoli_wgs\tSM:ecoli_wgs\tPL:ILLUMINA" \
    NC_000913.3.fa reads_R1.fastq.gz reads_R2.fastq.gz 2> bwa.log \
    | samtools sort -@ "$(nproc)" -o "$BAM" -
  samtools index "$BAM"
fi

# --- summary ----------------------------------------------------------
n_reads=$(samtools view -c -F 0x100 -F 0x400 "$BAM")
mapped=$(samtools view -c -F 0x104 -F 0x400 "$BAM")
mean_depth=$(samtools depth -a "$BAM" 2>/dev/null \
  | awk '{s+=$3; n++} END{if(n) printf "%.1f", s/n; else print "0"}')
echo
echo "[ecoli-wgs] ready:"
ls -lh "$BAM" NC_000913.3.fa | awk '{print "    "$0}'
printf "Reads (primary, non-dup): %d  mapped: %d  mean depth: %s\n" \
  "$n_reads" "$mapped" "$mean_depth"
echo
echo "Next:"
echo "  parity/upstream/bin/lofreq indelqual --dindel -f NC_000913.3.fa -o ecoli_wgs.iq.bam ${BAM}"
echo "  samtools index ecoli_wgs.iq.bam"
echo "  scripts/compare.sh \\"
echo "    $OUT/ecoli_wgs.iq.bam \\"
echo "    $OUT/NC_000913.3.fa \\"
echo "    parity/compare/ecoli_wgs_30x"
