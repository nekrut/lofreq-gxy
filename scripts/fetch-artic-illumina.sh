#!/usr/bin/env bash
# Fetch + align a paired-end Illumina SARS-CoV-2 ARTIC(-style) amplicon
# sample for Tier 2.1 parity testing.
#
# Usage:
#   scripts/fetch-artic-illumina.sh <name> <r1_url> <r2_url>
#
# Produces, under parity/fixtures/artic_<name>/:
#   reads_R1.fastq.gz       downloaded R1
#   reads_R2.fastq.gz       downloaded R2
#   NC_045512.2.fa{.fai,.bwt,...}  SARS-CoV-2 reference + bwa index
#   artic_<name>.sorted.bam coordinate-sorted, indexed BAM
#   artic_<name>.sorted.bam.bai
#
# Requires: bwa, samtools, curl (all via the lofreq-gxy-test conda env).
# Idempotent: re-running skips already-fetched / already-built files.
#
# The reference is NC_045512.2 (Wuhan-Hu-1). The alignment is the minimum
# needed to feed compare.sh — no primer trimming, no dedup. That matches
# what `lofreq-gxy` is designed to see as input; any primer bias is a
# property of the dataset and will show up in both tools identically.

set -euo pipefail

NAME="${1:?usage: fetch-artic-illumina.sh <name> <r1_url> <r2_url>}"
R1_URL="${2:?usage: fetch-artic-illumina.sh <name> <r1_url> <r2_url>}"
R2_URL="${3:?usage: fetch-artic-illumina.sh <name> <r1_url> <r2_url>}"

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
OUT="$ROOT/parity/fixtures/artic_$NAME"
REF_URL="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_045512.2&rettype=fasta"

mkdir -p "$OUT"
cd "$OUT"

for tool in bwa samtools curl; do
  command -v "$tool" >/dev/null || { echo "error: $tool not on PATH" >&2; exit 1; }
done

echo "[artic-$NAME] target dir: $OUT"

# --- reads ---------------------------------------------------------------
for side in 1 2; do
  url_var="R${side}_URL"
  url="${!url_var}"
  # ENA 'ftp.sra.ebi.ac.uk/...' entries need an explicit https:// prefix.
  [[ "$url" == ftp.* ]] && url="https://$url"
  dest="reads_R${side}.fastq.gz"
  if [[ ! -s "$dest" ]]; then
    echo "[artic-$NAME] fetching R${side}: $url"
    curl -fSL --max-time 600 "$url" -o "$dest"
  fi
done

# --- reference + bwa index ----------------------------------------------
if [[ ! -s NC_045512.2.fa ]]; then
  echo "[artic-$NAME] fetching NC_045512.2 reference"
  curl -fSL --max-time 60 "$REF_URL" -o NC_045512.2.fa
fi
[[ -s NC_045512.2.fa.fai ]] || samtools faidx NC_045512.2.fa
[[ -s NC_045512.2.fa.bwt ]] || bwa index NC_045512.2.fa

# --- align ---------------------------------------------------------------
BAM="artic_$NAME.sorted.bam"
if [[ ! -s "$BAM" ]]; then
  echo "[artic-$NAME] bwa mem → samtools sort"
  bwa mem -t "$(nproc)" -R "@RG\tID:artic_$NAME\tSM:artic_$NAME\tPL:ILLUMINA" \
    NC_045512.2.fa reads_R1.fastq.gz reads_R2.fastq.gz 2> bwa.log \
    | samtools sort -@ "$(nproc)" -o "$BAM" -
  samtools index "$BAM"
fi

# --- summary -------------------------------------------------------------
n_reads=$(samtools view -c -F 0x100 -F 0x400 "$BAM")
mapped=$(samtools view -c -F 0x104 -F 0x400 "$BAM")
echo
echo "[artic-$NAME] ready:"
ls -lh "$BAM" NC_045512.2.fa | awk '{print "    "$0}'
printf "Reads (primary, non-dup): %d  mapped: %d\n" "$n_reads" "$mapped"
echo
echo "Next:"
echo "  scripts/compare.sh \\"
echo "    $OUT/$BAM \\"
echo "    $OUT/NC_045512.2.fa \\"
echo "    parity/compare/artic_$NAME"
