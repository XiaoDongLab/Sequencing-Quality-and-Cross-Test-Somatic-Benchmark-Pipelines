# Sequencing Quality and Cross-Test Somatic Benchmark Pipelines

This repository contains two independent workflows:

1. `Sequencing Metrics Pipeline` for sequencing quality assessment from BAM and VCF inputs.
2. `Cross-Test Somatic Benchmark Pipeline` for estimating somatic calling performance from paired call sets.

# Environment Requirements

For the sequencing metrics workflow:

- Command-line tools: `samtools`, `bcftools`, `bam-lorenz-coverage`, `Rscript`, `awk`, `sed`, `bc`
- R packages: `ggplot2`, `dplyr`, `tidyr`, `readr`, `patchwork`, `stringr`, `ggbeeswarm`, `scales`

Conda setup:

```bash
conda env create -f sequencing_metrics.yml
conda activate sequencing_metrics
```

# Repository Layout

```text
run_sequencing_metrics_pipeline.sh
plot_sequencing_metrics.R
run_sm.sh
plot_sm.R
sequencing_metrics.yml
README.md
```

# Pipeline 1: Sequencing Metrics

Purpose: compute reads-versus-coverage, Lorenz/Gini, depth-distribution, chromosome-window coverage, and VAF summaries across multiple samples.

Inputs:

- A sample list file with one sample ID per line
- Per-sample directory: `../lab_data/<sample_id>/`
- Required files per sample:
  `<sample_id>_WGNS_filtered.bam`, `ssnv.results.muts.vcf.gz`, `sindel.results.indel.vcf.gz`
- Reference files: `reference.fa`, `reference.fa.fai`

Run:

```bash
bash run_sequencing_metrics_pipeline.sh sample_ids.txt
```

Full form:

```bash
bash run_sequencing_metrics_pipeline.sh \
  sample_ids.txt \
  8 \
  ../lab_data \
  ./sequencing_metrics_results \
  reference.fa
```

Outputs:

- Per sample: coverage tables, Lorenz tables, depth distribution, chromosome coverage windows, filtered VCFs, and AD-tagged VCFs
- Aggregated tables: `sequencing_metrics_results/bridge_tsv/`
- Plots: `sequencing_metrics_results/plots/`

Key plot files:

- `reads_coverage_curve.pdf`
- `lorenz_curve_and_gini.pdf`
- `depth_distribution_0_100.pdf`
- `chrom_window_coverage_all_samples.pdf`
- `vaf_raw_including_zero.pdf`
- `vaf_raw_excluding_zero.pdf`

# Pipeline 2: Cross-Test Somatic Benchmark

Purpose: compare somatic call sets across tests and estimate pseudo-true-positive rates, false-positive rates, mixture behavior, and burden correction.

Rationale: the workflow uses cross-test agreement to identify reproducible signal and uses germline-derived expectations plus a pseudo-false-positive reference to quantify likely noise. This makes the resulting precision and burden estimates more interpretable than using raw call counts alone.

Inputs:

- Required VCFs: `somatic.vcf.gz`, `cross_test_somatic.vcf.gz`, `germline_sample.vcf.gz`, `germline_test.vcf.gz`, `pseudo_fp_reference.vcf.gz`
- Optional reference: `reference.fa`
- All VCFs must be compressed and indexed

Run:

```bash
bash run_sm.sh \
  --sample SAMPLE \
  --callable 2500000000 \
  --somatic som.vcf.gz \
  --cross cross.vcf.gz \
  --germS gS.vcf.gz \
  --germT gT.vcf.gz \
  --fpref ref.vcf.gz \
  --out results
```

Optional arguments: `--filter PASS`, `--ref reference.fa`, `--reflen 3137300923`

Outputs:

- Tables: `<sample>.x.tsv`, `<sample>.m.tsv`, `<sample>.snv_xy.tsv`, `<sample>.ind_xy.tsv`, `<sample>.mix.tsv`, `<sample>.burden.tsv`
- Plots: `<sample>.snv_tp_vs_fp.pdf`, `<sample>.ind_tp_vs_fp.pdf`, `<sample>.mix_inferred.pdf`, `<sample>.burden_snv.pdf`, `<sample>.burden_ind.pdf`
- Temporary `bcftools isec` files: `results/tmp/`

# License

For academic research use.
