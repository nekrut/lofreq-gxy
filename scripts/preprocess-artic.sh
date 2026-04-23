#!/usr/bin/env bash
# Preprocess an ARTIC-style amplicon BAM for fair parity comparison:
# primer-trim with ivar + add indel quality scores with lofreq indelqual
# (the `--dindel` preprocessing that upstream lofreq's docs recommend
# before `lofreq call` on indel-bearing datasets).
#
# Usage:
#   scripts/preprocess-artic.sh <name> <primer_bed>
#
# Reads:
#   parity/fixtures/artic_<name>/artic_<name>.sorted.bam
#   parity/fixtures/artic_<name>/NC_045512.2.fa
# Writes:
#   parity/fixtures/artic_<name>/artic_<name>.trim.iq.sorted.bam (+ .bai)
#
# The preprocessed BAM is then fed to scripts/compare.sh so both
# lofreq-gxy and upstream lofreq see the same primer-trimmed,
# indelqual-augmented input.

set -euo pipefail

NAME="${1:?usage: preprocess-artic.sh <name> <primer_bed>}"
BED="${2:?usage: preprocess-artic.sh <name> <primer_bed>}"

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
DIR="$ROOT/parity/fixtures/artic_$NAME"
IN_BAM="$DIR/artic_$NAME.sorted.bam"
REF="$DIR/NC_045512.2.fa"
TRIM_PREFIX="$DIR/artic_$NAME.trim"
IQ_BAM="$DIR/artic_$NAME.trim.iq.sorted.bam"
UPSTREAM_LOFREQ="$ROOT/parity/upstream/bin/lofreq"

for tool in ivar samtools; do
  command -v "$tool" >/dev/null || { echo "error: $tool not on PATH" >&2; exit 1; }
done
[[ -x "$UPSTREAM_LOFREQ" ]] || { echo "error: $UPSTREAM_LOFREQ missing" >&2; exit 1; }
[[ -s "$IN_BAM" ]] || { echo "error: $IN_BAM missing — run fetch-artic-illumina.sh first" >&2; exit 1; }
[[ -s "$BED" ]] || { echo "error: primer BED $BED missing" >&2; exit 1; }

echo "[preprocess-$NAME] primer-trim with ivar"
# -e keeps reads with no primers (don't drop); -q 0 -m 30 = permissive
# defaults that match upstream ARTIC V3 Illumina pipelines.
ivar trim \
  -i "$IN_BAM" \
  -b "$BED" \
  -p "$TRIM_PREFIX" \
  -e -q 0 -m 30 2> "$DIR/ivar.log" >/dev/null

echo "[preprocess-$NAME] sort trimmed BAM"
samtools sort -@ "$(nproc)" -o "$TRIM_PREFIX.sorted.bam" "$TRIM_PREFIX.bam"
rm "$TRIM_PREFIX.bam"

echo "[preprocess-$NAME] lofreq indelqual --dindel"
"$UPSTREAM_LOFREQ" indelqual --dindel -f "$REF" \
  -o "$IQ_BAM" "$TRIM_PREFIX.sorted.bam"
samtools index "$IQ_BAM"
rm "$TRIM_PREFIX.sorted.bam"

echo
echo "[preprocess-$NAME] ready:"
ls -lh "$IQ_BAM" | awk '{print "    "$0}'
echo
echo "Next:"
echo "  scripts/compare.sh \\"
echo "    $IQ_BAM \\"
echo "    $REF \\"
echo "    parity/compare/artic_${NAME}_pp"
