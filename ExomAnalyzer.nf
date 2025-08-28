params.input = ''
params.output = '~/'
params.reference = ''
params.dbsnp = ''
params.bed = ''
params.index = ''


log.info """
============================================
  PIPELINE: Analyze Exom
  VERSION: 1.0.0
  AUTHOR: Isson
  DESCRIPTION: SNP-calling of the exome from the Sociogenetic Engineering project
  
  EXOM: ${params.input}
  OUTPUT_DIR: ${params.output}
  REFERENCE: ${params.reference}
  KNOWN_SITES: ${params.dbsnp}
  TARGET_REGIONS: ${params.bed}
============================================
"""

// This functuion return error, if user does not specify mandatory parameters. 
// It does not check params.index, because it is not necessary.
def checkParams() {
    def errors = []
    
    if (!params.input) {
        errors << "--input (path to input files)"
    }
    if (!params.output) {
        errors << "--output (output directory)"
    }
    if (!params.reference){
        errors << "--reference (reference genome)"
    }
    if (!params.dbsnp){
        errors << "--dbsnp (path to file with known sites)"
    }
    if (!params.bed){
        errors << "--bed (target regions of exom)"
    }
    if (errors) {
        log.error "ОШИБКА: Не указаны обязательные параметры:\n" + 
                 errors.collect { "• $it" }.join("\n")
        System.exit(1)
    }
}
checkParams()


// Make quality analyze of sequence
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
    ls -lh
    fastqc -o . $reads
    """
}


process multiqcAnalyze {
    publishDir "${params.output}/QC", mode: 'copy'
    container 'multiqc/multiqc:latest'
    containerOptions '--user $(id -u):$(id -g)'
    
    input:
    path qual
    
    output:
    path "*.{html}"
    
    script:
    """
    multiqc .
    """
}


// If the user does not specify a path to their own index, the program will do it.
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
    samtools view -bS "$alignment" | samtools sort -o alignmentSort.bam
    """
}


process alignmentQC {
    publishDir "${params.output}/Converted", mode: 'copy'
    container "broadinstitute/gatk:latest"
    containerOptions '--user $(id -u):$(id -g)'

    input:
      path converted
      
    output:
      path "*"
      
    script:
    """
    samtools stats $converted > stats.txt
    """
}

process coverageQC {
    publishDir "${params.output}/Converted", mode: 'copy'
    container "gfanz/mosdepth:latest"
    containerOptions '--user $(id -u):$(id -g)'

    input:
      path converted
      path bed_file
      
    output:
      path "*"
      
    script:
    """
    mosdepth -b $bed_file -n sample alignmentSort $converted
    """
}


process faidxIReference {
    container "broadinstitute/gatk:latest"
    containerOptions '--user $(id -u):$(id -g)'
    
    input:
      path reference
      
    output:
      path "*"
      
    script:
    """
    samtools faidx $reference
    gatk CreateSequenceDictionary -R $reference -O GRCh38.p14.genome.dict
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
    samtools index $convert
    gatk AddOrReplaceReadGroups -I $convert -O alignmentSortRg.bam -RGID "Sample1_L1" -RGSM "Sample1" -RGPL "ILLUMINA" -RGLB "Lib1" -RGPU "L001"
    gatk MarkDuplicates -I alignmentSortRg.bam -O alignmentSortRgMarked.bam --METRICS_FILE duplicates_metrics.txt
    samtools index alignmentSortRgMarked.bam
    gatk IndexFeatureFile -I $sites
    gatk BaseRecalibrator -I alignmentSortRgMarked.bam -R $reference --known-sites $sites -L $bed_file -O data.table
    gatk ApplyBQSR -I alignmentSortRgMarked.bam -R $reference --bqsr-recal-file data.table -O alignmentSortRgMarkedQual.bam
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
    reads = Channel.fromPath(params.input).collect()
    reference = file(params.reference)
    sites = file(params.dbsnp)
    bed_file = file(params.bed)
    
    
    fastqc_analyze = fastqcAnalyze(reads)
    multiqc_analyze = multiqcAnalyze(fastqc_analyze)
    
    
    if (!params.index) {
      index = indexReference(reference)
    }
    else {
      index = Channel.fromPath("${params.index}/*").collect()
    }
    
    
    alignment = mappingSequence(reference, reads, index)
    convert = convertationOfMapping(alignment)
    
    
    alignment_qc = alignmentQC(convert)
    coverage_qc = coverageQC(convert, bed_file)
    
    dop_indexes = faidxIReference(reference)
    extra_proc = deduplicationAndRecalibration(convert, reference, sites, dop_indexes, bed_file)
    variant_calling = variantCalling(reference, extra_proc, bed_file, dop_indexes)
    
}

