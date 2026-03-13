#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 1 ]]; then
  cat >&2 <<USAGE
Usage:
  bash run_sequencing_metrics_pipeline.sh <sample_id_file> [threads] [base_dir] [out_dir] [ref_fasta] [vaf_site_file]

Arguments:
  sample_id_file   Text file with one sample ID per line
  threads          Number of threads (default: 8; users should set this explicitly)
  base_dir         Base directory containing sample subdirectories (default: ../lab_data)
  out_dir          Output directory for merged tables and plots (default: ./sequencing_metrics_results)
  ref_fasta        Reference FASTA (default convenience path only; users should provide their own reference)
  vaf_site_file    VAF site file relative to each sample directory, or an absolute path.
                   Supported formats: VCF/VCF.GZ or BED/BED.GZ with columns:
                   chrom, start, end, ref, alt. Default: hSNP.bed
USAGE
  exit 1
fi

sample_file="$1"
threads="${2:-8}"
base_dir="${3:-../lab_data}"
out_dir="${4:-./sequencing_metrics_results}"
ref_fasta="${5:-$HOME/ref/hg38/references-hg38-v0-Homo_sapiens_assembly38.fasta}"
vaf_site_spec="${6:-hSNP.bed}"
fai="${ref_fasta}.fai"
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
r_script="${script_dir}/plot_sequencing_metrics.R"
bridge_dir="${out_dir}/bridge_tsv"
plot_dir="${out_dir}/plots"
genome_size="3137300923"

for exe in samtools bcftools bam-lorenz-coverage Rscript awk sed bc; do
  command -v "$exe" >/dev/null 2>&1 || { echo "Missing dependency: $exe" >&2; exit 1; }
done

[[ -f "$sample_file" ]] || { echo "Sample ID file not found: $sample_file" >&2; exit 1; }
[[ -f "$ref_fasta" ]] || { echo "Reference FASTA not found: $ref_fasta" >&2; exit 1; }
[[ -f "$fai" ]] || { echo "Reference index not found: $fai" >&2; exit 1; }
[[ -f "$r_script" ]] || { echo "R script not found: $r_script" >&2; exit 1; }

mkdir -p "$bridge_dir" "$plot_dir"

if [[ $# -lt 5 ]]; then
  echo "Warning: the default reference FASTA path is a convenience placeholder and is not recommended for routine use." >&2
  echo "Warning: provide the study-specific reference FASTA explicitly." >&2
fi

if [[ $# -lt 6 ]]; then
  echo "Warning: the default VAF site file is hSNP.bed in each sample directory." >&2
  echo "Warning: confirm that this default matches the study design before running the pipeline." >&2
fi

resolve_sample_path() {
  local samp="$1"
  local sdir="$2"
  local spec="$3"
  local resolved="${spec//\{sample_id\}/$samp}"
  if [[ "$resolved" = /* ]]; then
    printf '%s\n' "$resolved"
  else
    printf '%s/%s\n' "$sdir" "$resolved"
  fi
}

build_vaf_targets_vcf() {
  local src="$1"
  local out_vcf="$2"
  local src_lc="${src,,}"

  if [[ "$src_lc" == *.vcf || "$src_lc" == *.vcf.gz ]]; then
    bcftools view -Ov "$src" | bcftools norm -m -any -Ov -o "$out_vcf"
    return 0
  fi

  if [[ "$src_lc" == *.bed || "$src_lc" == *.bed.gz ]]; then
    {
      printf '##fileformat=VCFv4.2\n'
      printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
      if [[ "$src_lc" == *.gz ]]; then
        gzip -cd "$src"
      else
        cat "$src"
      fi | awk '
        BEGIN{FS=OFS="\t"}
        /^#/ {next}
        NF < 5 {next}
        $2 !~ /^[0-9]+$/ {next}
        $3 !~ /^[0-9]+$/ {next}
        $4 == "" || $4 == "." {next}
        $5 == "" || $5 == "." {next}
        {print $1, $2 + 1, ".", toupper($4), toupper($5), ".", "PASS", "."}
      '
    } > "$out_vcf"

    if [[ "$(awk 'END{print NR}' "$out_vcf")" -le 2 ]]; then
      echo "No usable VAF target records were found in: $src" >&2
      return 1
    fi
    return 0
  fi

  echo "Unsupported VAF site file format: $src" >&2
  echo "Supported formats are .vcf, .vcf.gz, .bed, and .bed.gz." >&2
  return 1
}

run_one_sample() {
  local samp="$1"
  local sdir="${base_dir}/${samp}"
  local bam="${sdir}/${samp}_WGNS_filtered.bam"
  local cov_txt="${sdir}/cov_${samp}.txt"
  local vaf_site_file
  local vaf_targets_vcf="${sdir}/${samp}_vaf_targets.vcf"
  local vaf_ad_vcf="${sdir}/${samp}_vaf_AD.vcf"
  local vaf_values_tsv="${sdir}/${samp}_vaf_values.tsv"

  [[ -d "$sdir" ]] || { echo "Sample directory not found: $sdir" >&2; return 1; }
  [[ -f "$bam" ]] || { echo "BAM not found: $bam" >&2; return 1; }
  vaf_site_file="$(resolve_sample_path "$samp" "$sdir" "$vaf_site_spec")"
  [[ -f "$vaf_site_file" ]] || { echo "VAF site file not found: $vaf_site_file" >&2; return 1; }

  date
  echo "== ${samp}: START Depth =="
  local reads frac chr chr_len start end down_bam
  reads=$(samtools view -@ "$threads" -c "$bam")
  : > "$cov_txt"
  for j in $(seq 20000000 20000000 400000000); do
    frac=$(awk -v j="$j" -v reads="$reads" 'BEGIN{printf "%.9f", j/reads}')
    down_bam="${sdir}/downsample_${j}.bam"
    samtools view -@ "$threads" -s "$frac" -b "$bam" > "$down_bam"
    echo "$samp $j $(samtools depth -@ "$threads" "$down_bam" | awk 'END{print NR+0}')" >> "$cov_txt"
    rm -f "$down_bam"
  done
  echo "== ${samp}: END Depth =="
  date

  echo "== ${samp}: START Lorenz =="
  bam-lorenz-coverage \
    -l "${sdir}/${samp}_lorenz.txt" \
    -c "${sdir}/${samp}_coverage.txt" \
    -L "${sdir}/${samp}_lorenz.jpg" \
    -C "${sdir}/${samp}_coverage.jpg" \
    "$bam"
  echo "== ${samp}: END Lorenz =="
  date

  echo "== ${samp}: START Count =="
  samtools depth -@ "$threads" -g 0x400 -a "$bam" \
    | awk '{count[$3]++} END {for (num in count) print num, count[num]}' \
    > "${sdir}/${samp}_distribution.txt"
  echo "== ${samp}: END Count =="
  date

  echo "== ${samp}: START Coverage =="
  samtools index -@ "$threads" "$bam"
  : > "${sdir}/chrom_cov_table.txt"
  for j in {1..22} X; do
    chr="chr${j}"
    chr_len=$(awk -v c="$chr" '$1==c{print $2}' "$fai")
    [[ -n "$chr_len" ]] || continue
    for ((start=1; start<=chr_len; start+=1000000)); do
      end=$((start + 1000000 - 1))
      [[ "$end" -gt "$chr_len" ]] && end="$chr_len"
      samtools coverage -r "${chr}:${start}-${end}" "$bam" | sed -n '2p' >> "${sdir}/chrom_cov_table.txt"
    done
  done
  echo "== ${samp}: END Coverage =="
  date

  echo "== ${samp}: START VAF =="
  build_vaf_targets_vcf "$vaf_site_file" "$vaf_targets_vcf"
  bcftools mpileup -f "$ref_fasta" -R "$vaf_targets_vcf" -a FORMAT/AD -Ov "$bam" \
    | bcftools norm -m -any -Ov -o "$vaf_ad_vcf"
  bcftools query -f '%REF\t%ALT[\t%AD]\n' "$vaf_ad_vcf" \
    | awk '
      BEGIN{FS=OFS="\t"}
      function classify_variant(ref, alt) {
        if (length(ref) == 1 && length(alt) == 1) return "SNV"
        return "indel"
      }
      NF >= 3 {
        n = split($3, ad, ",")
        if (n < 2) next
        if (ad[1] == "." || ad[2] == ".") next
        denom = ad[1] + ad[2]
        if (denom <= 0) next
        print classify_variant($1, $2), ad[2] / denom
      }
    ' > "$vaf_values_tsv"

  echo "== ${samp}: END VAF =="
  date
}

collect_tables() {
  local depth_curve_tsv="${bridge_dir}/depth_curve.tsv"
  local lorenz_tsv="${bridge_dir}/lorenz_curve.tsv"
  local depth_dist_tsv="${bridge_dir}/depth_distribution.tsv"
  local chrom_cov_tsv="${bridge_dir}/chrom_coverage_1Mb.tsv"
  local vaf_tsv="${bridge_dir}/vaf_values.tsv"

  printf "sample_id\treads_number\tcovered_bases\tgenome_size\tcoverage_frac\n" > "$depth_curve_tsv"
  printf "sample_id\tx_frac\ty_frac\n" > "$lorenz_tsv"
  printf "sample_id\tdepth\tcount\n" > "$depth_dist_tsv"
  printf "sample_id\tchrom\tstart\tend\tcoverage\tmean_depth\tbin_mb\n" > "$chrom_cov_tsv"
  printf "sample_id\tvariant_class\tvaf\n" > "$vaf_tsv"

  while IFS= read -r samp || [[ -n "$samp" ]]; do
    [[ -n "$samp" ]] || continue
    local sdir="${base_dir}/${samp}"
    local cov_txt="${sdir}/cov_${samp}.txt"
    local lorenz_txt="${sdir}/${samp}_lorenz.txt"
    local dist_txt="${sdir}/${samp}_distribution.txt"
    local chrom_txt="${sdir}/chrom_cov_table.txt"
    local vaf_values="${sdir}/${samp}_vaf_values.tsv"

    [[ -s "$cov_txt" ]] && awk -v OFS='\t' -v gs="$genome_size" 'NF>=3{print $1,$2,$3,gs,$3/gs}' "$cov_txt" >> "$depth_curve_tsv"

    [[ -s "$lorenz_txt" ]] && awk -v OFS='\t' -v samp="$samp" '
      BEGIN{ix=0;iy=0}
      NR==1{
        for(i=1;i<=NF;i++){
          if($i ~ /X-fraction-sequenced-bases/) ix=i
          if($i ~ /Y-fraction-genome-covered/) iy=i
        }
        next
      }
      ix>0 && iy>0 {print samp,$ix,$iy}
    ' "$lorenz_txt" >> "$lorenz_tsv"

    [[ -s "$dist_txt" ]] && awk -v OFS='\t' -v samp="$samp" 'NF>=2{print samp,$1,$2}' "$dist_txt" >> "$depth_dist_tsv"

    [[ -s "$chrom_txt" ]] && awk -v OFS='\t' -v samp="$samp" '
      NF>=7{
        bin=int(($2-1)/1000000)+1
        print samp,$1,$2,$3,$6,$7,bin
      }
    ' "$chrom_txt" >> "$chrom_cov_tsv"

    [[ -s "$vaf_values" ]] && awk -v OFS='\t' -v samp="$samp" '
      NF >= 2 && $2 != "NA" && $2 != "nan" && $2 != "NaN" {print samp,$1,$2}
    ' "$vaf_values" >> "$vaf_tsv"
  done < "$sample_file"
}

while IFS= read -r samp || [[ -n "$samp" ]]; do
  [[ -n "$samp" ]] || continue
  run_one_sample "$samp"
done < "$sample_file"

collect_tables
Rscript "$r_script" "$bridge_dir" "$plot_dir"

echo "[OK] Pipeline finished."
echo "Bridge tables: $bridge_dir"
echo "Plots: $plot_dir"
