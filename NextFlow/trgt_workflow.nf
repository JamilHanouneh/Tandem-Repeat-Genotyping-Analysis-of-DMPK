#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// ---------- Params ----------
params.bam_dir      = "/home/jamil/Documents/NextFLow/PBMM2-BAM-Input-For-IGV-And-TRGT"
params.ref_genome   = "reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
params.repeats_bed  = "reference/pathogenic_repeats.hg38.bed"
params.out_dir      = "results"
params.trgt_cpus    = 2   // override with --trgt_cpus N if you have more cores

log.info """
--- TRGT Genotyping Workflow ---
BAM Directory : ${params.bam_dir}
Reference     : ${params.ref_genome}
Repeats       : ${params.repeats_bed}
Output        : ${params.out_dir}
TRGT threads  : ${params.trgt_cpus}
----------------------------------
"""

// ---------- Workflow ----------
workflow {
  Channel.fromPath(params.ref_genome).set { genome_ch }
  Channel.fromPath(params.repeats_bed).set { repeats_ch }

  // Build FASTA index once
  index_ch = INDEX_GENOME(genome_ch)

  // Discover BAMs (don’t trust existing .bai — rebuild)
  Channel
    .fromPath("${params.bam_dir}/*.pbmm2.repeats.bam")
    .map { Path bam ->
      def sample = bam.getBaseName().replace('.pbmm2.repeats','')
      log.info "Queued sample: ${sample}"
      tuple(sample, bam)
    }
    .set { bam_ch }

  // Ensure a fresh .bai for each BAM
  indexed_bam_ch = ENSURE_BAI(bam_ch)

  // Run TRGT on each sample
  RUN_TRGT(genome_ch, index_ch, repeats_ch, indexed_bam_ch)
}

// ---------- Processes ----------

process INDEX_GENOME {
  tag "faidx"
  publishDir "${params.out_dir}/reference_index", mode: 'copy'

  input:
    path genome_fasta

  output:
    path "${genome_fasta}.fai"

  cpus 1
  conda 'bioconda::samtools'

  script:
  """
  samtools faidx ${genome_fasta}
  """
}

process ENSURE_BAI {
  tag { sample }
  publishDir "${params.out_dir}/reindexed_bam", mode: 'copy'

  input:
    tuple val(sample), path(bam)

  output:
    tuple val(sample), path(bam), path("${bam}.bai")

  cpus 1
  conda 'bioconda::samtools'

  script:
  """
  # Rebuild a clean index (overwrites if present)
  samtools quickcheck -v ${bam} || echo "samtools quickcheck reported issues, attempting reindex" >&2
  samtools index -b ${bam} ${bam}.bai
  """
}

process RUN_TRGT {
  tag { sample }
  publishDir "${params.out_dir}/trgt", mode: 'copy'

  // Let other samples continue even if one fails
  errorStrategy 'finish'
  maxRetries 1

  input:
    path genome_fasta
    path genome_index
    path repeats_bed
    tuple val(sample), path(bam), path(bai)

  output:
    path "${sample}.vcf.gz",       emit: vcfgz
    path "${sample}.spanning.bam", emit: spans
    path "${sample}.trgt.log",     emit: logs

  cpus params.trgt_cpus as int
  conda 'bioconda::trgt bioconda::samtools'

  script:
  """
  set -euo pipefail

  # Run TRGT and capture full log
  trgt genotype \
      --genome ${genome_fasta} \
      --repeats ${repeats_bed} \
      --reads ${bam} \
      --output-prefix ${sample} \
      --threads ${task.cpus} \
      &> ${sample}.trgt.log

  # Be explicit if nothing was produced
  if [[ ! -s ${sample}.vcf.gz ]]; then
    echo "[TRGT] No VCF was produced for ${sample}. See details above." >> ${sample}.trgt.log
    exit 2
  fi
  """
}
