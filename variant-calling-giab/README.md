# Germline variant calling on GIAB HG002 chr20: bcftools vs GATK HaplotypeCaller vs DeepVariant (in progress)

Tracking issue: [#14](https://github.com/shoo99/sbmlab-benchmarks/issues/14). No post yet.

**Status: running.** The scripts are in place. `happy/bcftools_35x.extended.csv` and `timings.tsv` hold the first result (bcftools). GATK HaplotypeCaller and DeepVariant are running, and this README gets the full results when they finish.

- **Data**: GIAB HG002 Illumina NovaSeq PCR-free WGS, 35x, BWA-MEM to GRCh38 no-alt, chr20. This is the DeepVariant case-study BAM.
- **Truth**: GIAB v4.2.1 benchmark VCF and regions.
- **Strata**: GIAB stratifications v3.3. They are all difficult regions, not in any difficult region, tandem repeats and homopolymers, segmental duplications, and low mappability.
- **Scoring**: hap.py 0.3.15 with the vcfeval engine.
- **Scripts**:
  - `download.sh` fetches the inputs.
  - `run.sh` runs the steps prepare, bcftools, gatk, deepvariant and happy.
  - `score.py` writes `scores.csv`, and `figure.py` draws the figure.
