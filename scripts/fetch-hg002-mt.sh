#!/usr/bin/env bash
# Fetch the HG002 mitochondrial smoke-test fixture (Tier 1.2).
#
# MT is effectively haploid with real heteroplasmy — it's the only GIAB
# subset whose truth set maps cleanly onto a low-AF haploid caller.
# Produces, under parity/fixtures/hg002_mt/:
#
#   hg002_mt.bam          chrM reads from HG002 30x HiSeq (~5–30 MB)
#   hg002_mt.bam.bai
#   rCRS.fa               rCRS / GRCh38 chrM (~16 kb)
#   rCRS.fa.fai
#   hg002_mt_truth.vcf    GIAB HG002 MT v1.00 truth VCF
#
# Requires: samtools, bcftools, curl, wget, NCBI/NIH reachable.
# Idempotent: re-running skips already-downloaded files.
#
# Data sources (pinned where possible; URLs may shift over time —
# check TEST-PLAN.md for updates):
#   BAM:    NIST/GIAB HG002 30x HiSeq Illumina, aligned to GRCh38
#   Truth:  GIAB HG002 CMRG v1.00 — MT subset
#   Ref:    rCRS (NC_012920.1) — the GRCh38 chrM

set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
OUT="$ROOT/parity/fixtures/hg002_mt"
mkdir -p "$OUT"
cd "$OUT"

# --- Source URLs (pinned) --------------------------------------------------
# HG002 GRCh38 BAM — 50x NovaSeq Illumina from the DeepVariant case-study
# corpus. Full BAM is ~118 GB, but its colocated .bai means
# `samtools view -b <url> chrM` pulls only chrM reads via HTTP Range.
BAM_URL="${HG002_BAM_URL:-https://storage.googleapis.com/deepvariant/case-study-testdata/HG002_NIST_150bp_50x.bam}"

# GIAB HG002 CMRG v1.00 benchmark. CMRG is "Challenging Medically Relevant
# Genes" — v1.00 includes mitochondrial variants as part of the small-
# variant set.
CMRG_VCF_URL="${HG002_CMRG_URL:-https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/CMRG_v1.00/GRCh38/SmallVariant/HG002_GRCh38_CMRG_smallvar_v1.00.vcf.gz}"

# Mitochondrial reference. The exact sequence depends on the assembly the
# BAM was aligned against:
#
#   GRCh37 MT contig   = NC_012920.1 exactly (rCRS).
#   GRCh38 chrM contig = NC_012920.1 with 'N' at pos 3107 (UCSC ships this).
#
# We auto-detect which contig name the BAM uses below and pick the right
# reference. NCBI is the fallback; UCSC is preferred for GRCh38.
CHRM_UCSC_URL="${CHRM_UCSC_URL:-https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chrM.fa.gz}"
CHRM_NCBI_URL="${CHRM_NCBI_URL:-https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_012920.1&rettype=fasta&retmode=text}"

# --- Preflight -------------------------------------------------------------
for tool in samtools bcftools curl wget; do
  command -v "$tool" >/dev/null || { echo "error: $tool not on PATH" >&2; exit 1; }
done

echo "[hg002-mt] target dir: $OUT"

# --- Discover the MT contig name the BAM uses (chrM vs MT) -----------------
# We do this first so the reference we fetch matches the BAM's header.
if [[ ! -s bam_header.sam ]]; then
  echo "[hg002-mt] probing BAM header for MT contig name"
  samtools view -H "$BAM_URL" > bam_header.sam 2>/dev/null \
    || { echo "error: BAM header unreachable" >&2; exit 1; }
fi
MT_CONTIG=$(
  grep -E '^@SQ' bam_header.sam \
    | grep -oE 'SN:(chrM|MT|chrMT)[[:space:]]' \
    | head -1 \
    | sed 's/^SN://;s/[[:space:]]$//'
)
if [[ -z "$MT_CONTIG" ]]; then
  echo "error: BAM has no MT contig (checked chrM / MT / chrMT)" >&2
  exit 1
fi
echo "[hg002-mt]   BAM uses contig name: $MT_CONTIG"

# --- Reference (mt sequence matching the BAM's contig name) ---------------
if [[ ! -s chrM.fa ]]; then
  rm -f chrM.fa.gz chrM.fa.raw
  # For GRCh38 (`chrM`) prefer UCSC — ships rCRS with 'N' at pos 3107.
  # For GRCh37 (`MT`), NCBI NC_012920.1 is an exact match.
  if [[ "$MT_CONTIG" == "chrM" ]]; then
    echo "[hg002-mt] fetching GRCh38 chrM from UCSC"
    if curl -fSL --max-time 60 "$CHRM_UCSC_URL" -o chrM.fa.gz 2>/dev/null \
        && [[ -s chrM.fa.gz ]]; then
      gunzip -f chrM.fa.gz
    else
      rm -f chrM.fa.gz
      echo "[hg002-mt]   UCSC unreachable; falling back to NCBI NC_012920.1"
      echo "[hg002-mt]   WARNING: coords near pos 3107 will NOT match the BAM."
      curl -fSL --max-time 60 "$CHRM_NCBI_URL" -o chrM.fa.raw
      [[ -s chrM.fa.raw ]] || { echo "error: both reference sources failed" >&2; exit 1; }
      awk 'NR==1{print ">chrM"; next} {print}' chrM.fa.raw > chrM.fa
      rm chrM.fa.raw
    fi
  else
    echo "[hg002-mt] fetching MT from NCBI NC_012920.1"
    curl -fSL --max-time 60 "$CHRM_NCBI_URL" -o chrM.fa.raw
    [[ -s chrM.fa.raw ]] || { echo "error: NCBI unreachable" >&2; exit 1; }
    awk -v name="$MT_CONTIG" 'NR==1{print ">"name; next} {print}' chrM.fa.raw > chrM.fa
    rm chrM.fa.raw
  fi
  samtools faidx chrM.fa
fi

# --- BAM (MT only via remote-indexed samtools view) ------------------------
if [[ ! -s hg002_mt.bam ]]; then
  echo "[hg002-mt] extracting $MT_CONTIG reads from remote HG002 BAM"
  echo "[hg002-mt]   (pulls only $MT_CONTIG via HTTP Range; no full BAM download)"
  samtools view -b "$BAM_URL" "$MT_CONTIG" > hg002_mt.bam.unsorted
  n_reads=$(samtools view -c hg002_mt.bam.unsorted 2>/dev/null || echo 0)
  if [[ "$n_reads" -eq 0 ]]; then
    rm -f hg002_mt.bam.unsorted
    echo "error: zero reads extracted — BAM_URL probably unreachable or index missing" >&2
    exit 1
  fi
  samtools sort -o hg002_mt.bam hg002_mt.bam.unsorted
  rm hg002_mt.bam.unsorted
  samtools index hg002_mt.bam
fi

# --- Truth VCF (optional) --------------------------------------------------
# HISTORICAL NOTE: At the time of writing, GIAB's CMRG v1.00 release does
# *not* include MT variants — CMRG's scope is autosomal medically-relevant
# genes, and the GIAB v4.2.1 small-variant benchmark excludes chrM by
# design. A dedicated HG002 MT truth set is in discussion but not yet
# released.
#
# We attempt to fetch CMRG in case a later release adds MT. If it has no
# MT variants (or is unreachable), the harness falls back to parity-only
# mode (compare.sh), which is still valuable — it catches regressions
# against upstream on real Illumina MT reads with real heteroplasmy.
if [[ ! -s hg002_mt_truth.vcf ]]; then
  echo "[hg002-mt] attempting CMRG v1.00 truth fetch (may have no MT variants)"
  if curl -fSL --max-time 120 "$CMRG_VCF_URL" -o cmrg.vcf.gz 2>/dev/null \
      && [[ -s cmrg.vcf.gz ]]; then
    curl -fSL --max-time 30 "${CMRG_VCF_URL}.tbi" -o cmrg.vcf.gz.tbi 2>/dev/null || true
    [[ -s cmrg.vcf.gz.tbi ]] || tabix -p vcf cmrg.vcf.gz
    bcftools view -r "$MT_CONTIG" cmrg.vcf.gz 2>/dev/null > hg002_mt_truth.vcf \
      || rm -f hg002_mt_truth.vcf
    rm -f cmrg.vcf.gz cmrg.vcf.gz.tbi
    if [[ ! -s hg002_mt_truth.vcf ]] || [[ $(grep -cv '^#' hg002_mt_truth.vcf) -eq 0 ]]; then
      echo "[hg002-mt]   CMRG v1.00 has no $MT_CONTIG variants (expected for this release)."
      echo "[hg002-mt]   Harness will run in parity-only mode."
      rm -f hg002_mt_truth.vcf
    fi
  else
    echo "[hg002-mt]   CMRG fetch failed; parity-only mode."
    rm -f hg002_mt_truth.vcf
  fi
fi

# --- Summary ---------------------------------------------------------------
truth_count=0
[[ -s hg002_mt_truth.vcf ]] && truth_count=$(grep -cv '^#' hg002_mt_truth.vcf)

echo
echo "[hg002-mt] ready:"
ls -lh hg002_mt.bam chrM.fa 2>/dev/null | awk '{print "    "$0}'
[[ -s hg002_mt_truth.vcf ]] && ls -lh hg002_mt_truth.vcf | awk '{print "    "$0}'
echo
echo "Reads in MT BAM ($MT_CONTIG): $(samtools view -c hg002_mt.bam)"
echo "Truth variants:              $truth_count"
echo
if (( truth_count > 0 )); then
  echo "Next: truth-aware comparison"
  echo "  scripts/compare-truth.sh \\"
  echo "    $OUT/hg002_mt.bam \\"
  echo "    $OUT/chrM.fa \\"
  echo "    $OUT/hg002_mt_truth.vcf \\"
  echo "    parity/compare/hg002_mt"
else
  echo "Next: parity-only comparison (no truth VCF available)"
  echo "  scripts/compare.sh \\"
  echo "    $OUT/hg002_mt.bam \\"
  echo "    $OUT/chrM.fa \\"
  echo "    parity/compare/hg002_mt"
fi
