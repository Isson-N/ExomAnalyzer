log.info """
============================================
  PIPELINE: Анализ экзома
  VERSION: 1.0.0
  AUTHOR: Isson
  DESCRIPTION: Всесторонний анализ
============================================
"""
params.input = ''
params.output = '~/'
params.reference = ''
params.dbsnp = ''
params.bed = ''


process fastqcAnalyze {
    publishDir "${params.output}/QC", pattern: '*.{html,zip}', mode: 'copy'
    container 'biocontainers/fastqc:v0.11.9_cv8'
    containerOptions '--user $(id -u):$(id -g)'
    
    input:
    path reads
    
    output:
    path "*.{html,zip}"
    
    script:
    """
    fastqc -o . "$reads"
    """
}

process indexReference {
    container 'biocontainers/bwa:v0.7.17_cv1'
    containerOptions '--user $(id -u):$(id -g)'
    
    input:
    path reference
    
    output:
    path "*"
    
    script:
    """
    bwa index $reference
    """
}



process mappingSequence {
    publishDir "${params.output}/Alignment", mode: 'copy'
    container 'biocontainers/bwa:v0.7.17_cv1'
    containerOptions '--user $(id -u):$(id -g)'
    
    input:
    path reference
    path reads
    path index
    
    output:
    path "*.sam"
    
    script:
    """
    bwa mem $reference $reads > output.sam
    """

}


process convertationOfMapping {
    publishDir "${params.output}/Converted", mode: 'copy'
    container "biocontainers/samtools:v1.9-4-deb_cv1"
    containerOptions '--user $(id -u):$(id -g)'

    input:
      path alignment
      
    output:
      path "*"
      
    script:
    """
    samtools view -bS $alignment > alignment.bam
    samtools sort alignment.bam -o alignmentSort.bam 
    rm alignment.bam
    """
}


process faidxIReference {
    container "biocontainers/samtools:v1.9-4-deb_cv1"
    containerOptions '--user $(id -u):$(id -g)'
    
    input:
      path reference
      path convert
      
    output:
      path "*"
      
    script:
    """
    samtools faidx $reference
    samtools index $convert
    """
}


process deduplicationAndRecalibration {
    publishDir "${params.output}/ExtraProcessingBam", mode: 'copy'
    container "broadinstitute/gatk:latest"
    containerOptions '--user $(id -u):$(id -g)'

    input:
      path convert
      path reference
      path sites
      path dop_indexes
      path bed_file
    
    output:
      path "*.bam"
    
    script:
    """
    gatk AddOrReplaceReadGroups -I $convert -O alignmentSortRg.bam -RGID "Sample1_L1" -RGSM "Sample1" -RGPL "ILLUMINA" -RGLB "Lib1" -RGPU "L001"
    gatk MarkDuplicates -I alignmentSortRg.bam -O alignmentSortRgMarked.bam --METRICS_FILE duplicates_metrics.txt
    samtools index alignmentSortRgMarked.bam
    gatk CreateSequenceDictionary -R $reference -O GRCh38.p14.genome.dict
    gatk IndexFeatureFile -I $sites
    gatk BaseRecalibrator -I alignmentSortRgMarked.bam -R $reference --known-sites $sites -L $bed_file -O data.table
    gatk ApplyBQSR -I alignmentSortRgMarked.bam -R $reference --bqsr-recal-file data.table --allow-missing-read-group -O alignmentSortRgMarkedQual.bam
    rm alignmentSortRg.bam alignmentSortRgMarked.bam
    """   
}


process variantCalling {
    publishDir "${params.output}/VCF", mode: 'copy'
    container "broadinstitute/gatk:latest"
    containerOptions '--user $(id -u):$(id -g)'
    
    input:
      path reference
      path align
      path bed_file
      path dop_indexes
      
    output:
      path "*.vcf"
    
    script:
    """
    gatk CreateSequenceDictionary -R $reference -O GRCh38.p14.genome.dict
    
    gatk HaplotypeCaller -R $reference -I $align -O results.vcf -L $bed_file
    gatk SelectVariants -R $reference -V results.vcf -select-type SNP -O raw_snps.vcf
    gatk SelectVariants -R $reference -V results.vcf -select-type INDEL -O raw_indels.vcf
    
    gatk VariantFiltration -R $reference -V raw_snps.vcf -O snps.vcf --filter-name "QD_filter" --filter-expression "QD < 2.0" --filter-name "FS_filter" --filter-expression "FS > 60.0" --filter-name "MQ_filter" --filter-expression "MQ < 40.0" --filter-name "SOR_filter" --filter-expression "SOR > 3.0" --filter-name "ReadPosRankSum_filter" --filter-expression "ReadPosRankSum < -8.0"
    
    gatk VariantFiltration -R $reference -V raw_indels.vcf -O indels.vcf --filter-name "QD_filter" --filter-expression "QD < 2.0" --filter-name "FS_filter" --filter-expression "FS > 200.0" --filter-name "SOR_filter" --filter-expression "SOR > 10.0"
    
    gatk MergeVcfs -I snps.vcf -I indels.vcf -O filtVar.vcf
    gatk SelectVariants -R $reference -V filtVar.vcf --exclude-filtered true -O passedVar.vcf
    rm results.vcf raw_snps.vcf raw_indels.vcf snps.vcf indels.vcf filtVar.vcf
    """
}


workflow {
    reads = Channel.fromPath(params.input)
    reference = file(params.reference)
    sites = file(params.dbsnp)
    bed_file = file(params.bed)
    
    fastqc_analyze = fastqcAnalyze(reads)
    index = indexReference(reference)
    alignment = mappingSequence(reference, reads, index)
    convert = convertationOfMapping(alignment)
    dop_indexes = faidxIReference(reference, convert)
    extra_proc = deduplicationAndRecalibration(convert, reference, sites, dop_indexes, bed_file)
    variant_calling = variantCalling(reference, extra_proc, bed_file, dop_indexes)
    
}

