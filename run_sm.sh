#!/usr/bin/env bash
set -euo pipefail

usage() {
cat <<EOF2

Single-sample cross-test somatic pipeline

Required:
  --sample      Sample ID
  --callable    Callable bases
  --somatic     Somatic VCF
  --cross       Cross-test somatic VCF
  --germS       Germline VCF (sample)
  --germT       Germline VCF (test)
  --fpref       VCF used to define pseudo FP
  --out         Output directory

Optional:
  --filter      FILTER string, e.g. PASS
  --ref         Reference FASTA for left-alignment
  --reflen      Reference genome length [default: 3137300923; convenience value only, not recommended for routine use]
  -h, --help    Show this message

Example:
  bash run_sm.sh \
    --sample ID \
    --callable 2500000000 \
    --somatic som.vcf.gz \
    --cross cross.vcf.gz \
    --germS gS.vcf.gz \
    --germT gT.vcf.gz \
    --fpref ref_for_fp.vcf.gz \
    --out outdir \
    --filter PASS \
    --ref reference.fa

Note:
  The example values for --callable and --reflen are placeholders.
  Users should replace them with study-specific callable-base counts and
  reference lengths.

EOF2
}

sample=""
callable=""
som=""
cross=""
germS=""
germT=""
fpref=""
out=""
flt=""
ref=""
refLen="3137300923"

if [[ $# -eq 0 ]]; then usage; exit 1; fi

while [[ $# -gt 0 ]]; do
  case "$1" in
    --sample)   sample="$2"; shift 2 ;;
    --callable) callable="$2"; shift 2 ;;
    --somatic)  som="$2"; shift 2 ;;
    --cross)    cross="$2"; shift 2 ;;
    --germS)    germS="$2"; shift 2 ;;
    --germT)    germT="$2"; shift 2 ;;
    --fpref)    fpref="$2"; shift 2 ;;
    --out)      out="$2"; shift 2 ;;
    --filter)   flt="$2"; shift 2 ;;
    --ref)      ref="$2"; shift 2 ;;
    --reflen)   refLen="$2"; shift 2 ;;
    -h|--help)  usage; exit 0 ;;
    *) echo "Unknown parameter: $1"; usage; exit 1 ;;
  esac
done

missing=0
[[ -z "$sample"   ]] && echo "Missing --sample"   && missing=1
[[ -z "$callable" ]] && echo "Missing --callable" && missing=1
[[ -z "$som"      ]] && echo "Missing --somatic"  && missing=1
[[ -z "$cross"    ]] && echo "Missing --cross"    && missing=1
[[ -z "$germS"    ]] && echo "Missing --germS"    && missing=1
[[ -z "$germT"    ]] && echo "Missing --germT"    && missing=1
[[ -z "$fpref"    ]] && echo "Missing --fpref"    && missing=1
[[ -z "$out"      ]] && echo "Missing --out"      && missing=1

if [[ "$missing" -eq 1 ]]; then
  echo ""
  usage
  exit 1
fi

need() { command -v "$1" >/dev/null 2>&1 || { echo "Missing: $1" >&2; exit 2; }; }
need bcftools
need Rscript

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
plot_script="$script_dir/plot_sm.R"

[[ -f "$plot_script" ]] || { echo "Missing: $plot_script" >&2; exit 2; }

mkdir -p "$out"/{vcf,tmp,tables,plots}

vview() {
  if [[ -n "$flt" ]]; then
    bcftools view -f "$flt" "$1"
  else
    bcftools view "$1"
  fi
}

vnorm() {
  local in="$1" outvcf="$2"
  if [[ -n "$ref" ]]; then
    vview "$in" | bcftools norm -f "$ref" -m -any -Oz -o "$outvcf"
  else
    vview "$in" | bcftools norm -m -any -Oz -o "$outvcf"
  fi
  bcftools index -t "$outvcf"
}

vsplt() {
  local in="$1" snv="$2" ind="$3"
  bcftools view -v snps -Oz -o "$snv" "$in"
  bcftools index -t "$snv"
  bcftools view -v indels -Oz -o "$ind" "$in"
  bcftools index -t "$ind"
}

cnt_vcf_records() {
  local f="$1"
  if [[ -f "$f" ]]; then
    bcftools view -H "$f" | wc -l | awk '{print $1}'
  else
    echo 0
  fi
}

vnorm "$som"   "$out/vcf/$sample.som.vcf.gz"
vnorm "$cross" "$out/vcf/$sample.cross.vcf.gz"
vnorm "$germS" "$out/vcf/$sample.gS.vcf.gz"
vnorm "$germT" "$out/vcf/$sample.gT.vcf.gz"
vnorm "$fpref" "$out/vcf/$sample.fpref.vcf.gz"

vsplt "$out/vcf/$sample.som.vcf.gz"   "$out/vcf/$sample.som.snv.vcf.gz"   "$out/vcf/$sample.som.ind.vcf.gz"
vsplt "$out/vcf/$sample.cross.vcf.gz" "$out/vcf/$sample.cross.snv.vcf.gz" "$out/vcf/$sample.cross.ind.vcf.gz"
vsplt "$out/vcf/$sample.fpref.vcf.gz" "$out/vcf/$sample.fpref.snv.vcf.gz" "$out/vcf/$sample.fpref.ind.vcf.gz"

bcftools view -g het -v snps -Oz -o "$out/vcf/$sample.gS.hetSNV.vcf.gz" "$out/vcf/$sample.gS.vcf.gz"
bcftools index -t "$out/vcf/$sample.gS.hetSNV.vcf.gz"

bcftools view -g het -v snps -Oz -o "$out/vcf/$sample.gT.hetSNV.vcf.gz" "$out/vcf/$sample.gT.vcf.gz"
bcftools index -t "$out/vcf/$sample.gT.hetSNV.vcf.gz"

d="$out/tmp/isec_guniq"
rm -rf "$d"
mkdir -p "$d"
bcftools isec -p "$d" -n=1 -w1 "$out/vcf/$sample.gS.vcf.gz" "$out/vcf/$sample.gT.vcf.gz" >/dev/null

bcftools view -Oz -o "$out/vcf/$sample.guniq.vcf.gz" "$d/0000.vcf"
bcftools index -t "$out/vcf/$sample.guniq.vcf.gz"

vsplt "$out/vcf/$sample.guniq.vcf.gz" "$out/vcf/$sample.guniq.snv.vcf.gz" "$out/vcf/$sample.guniq.ind.vcf.gz"

cnt_uniq() {
  local a="$1" b="$2" tag="$3"
  local dd="$out/tmp/isec_$tag"
  rm -rf "$dd"
  mkdir -p "$dd"
  bcftools isec -p "$dd" -n=1 -w1 "$a" "$b" >/dev/null
  cnt_vcf_records "$dd/0000.vcf"
}

pureS=$(cnt_uniq "$out/vcf/$sample.gS.hetSNV.vcf.gz" "$out/vcf/$sample.gT.hetSNV.vcf.gz" "pureS")
pureT=$(cnt_uniq "$out/vcf/$sample.gT.hetSNV.vcf.gz" "$out/vcf/$sample.gS.hetSNV.vcf.gz" "pureT")

cnt_tp() {
  local a="$1" b="$2" tag="$3"
  local dd="$out/tmp/isec_tp_$tag"
  rm -rf "$dd"
  mkdir -p "$dd"
  bcftools isec -p "$dd" -n=2 "$a" "$b" >/dev/null
  cnt_vcf_records "$dd/0002.vcf"
}

cnt_fp() {
  local a="$1" b="$2" tag="$3"
  local dd="$out/tmp/isec_fp_$tag"
  rm -rf "$dd"
  mkdir -p "$dd"
  bcftools isec -p "$dd" -n=1 -w1 "$a" "$b" >/dev/null
  cnt_vcf_records "$dd/0000.vcf"
}

ct_snv_tp=$(cnt_tp "$out/vcf/$sample.cross.snv.vcf.gz" "$out/vcf/$sample.guniq.snv.vcf.gz" "snv")
ct_ind_tp=$(cnt_tp "$out/vcf/$sample.cross.ind.vcf.gz" "$out/vcf/$sample.guniq.ind.vcf.gz" "ind")

ct_snv_fp=$(cnt_fp "$out/vcf/$sample.som.snv.vcf.gz" "$out/vcf/$sample.fpref.snv.vcf.gz" "snv")
ct_ind_fp=$(cnt_fp "$out/vcf/$sample.som.ind.vcf.gz" "$out/vcf/$sample.fpref.ind.vcf.gz" "ind")

ct_snv_tot=$(bcftools view -H "$out/vcf/$sample.cross.snv.vcf.gz" | wc -l | awk '{print $1}')
ct_ind_tot=$(bcftools view -H "$out/vcf/$sample.cross.ind.vcf.gz" | wc -l | awk '{print $1}')

gHetSNV=$(bcftools view -g het -v snps -H "$out/vcf/$sample.gS.vcf.gz" | wc -l | awk '{print $1}')
gHetIND=$(bcftools view -g het -v indels -H "$out/vcf/$sample.gS.vcf.gz" | wc -l | awk '{print $1}')

somSNV=$(bcftools view -H "$out/vcf/$sample.som.snv.vcf.gz" | wc -l | awk '{print $1}')
somIND=$(bcftools view -H "$out/vcf/$sample.som.ind.vcf.gz" | wc -l | awk '{print $1}')

xT="$out/tables/$sample.x.tsv"
mT="$out/tables/$sample.m.tsv"

{
  echo -e "id\tcb\trefLen\tct_snv_tot\tct_snv_tp\tct_snv_fp\tct_ind_tot\tct_ind_tp\tct_ind_fp\tgHetSNV\tgHetIND\tsomSNV\tsomIND"
  echo -e "$sample\t$callable\t$refLen\t$ct_snv_tot\t$ct_snv_tp\t$ct_snv_fp\t$ct_ind_tot\t$ct_ind_tp\t$ct_ind_fp\t$gHetSNV\t$gHetIND\t$somSNV\t$somIND"
} > "$xT"

{
  echo -e "id\tcb\trefLen\tgHetSNV\tpureS\tpureT"
  echo -e "$sample\t$callable\t$refLen\t$gHetSNV\t$pureS\t$pureT"
} > "$mT"

Rscript "$plot_script" "$xT" "$mT" "$out"

echo "OK"
