root="/Users/shutao/Desktop/nextflow_project"
params.index = "${root}/index.csv"
ref = "${root}/Ref/genome.fa"


Channel
    .fromPath(params.index)
    .splitCsv(header:true)
    .map{ row-> tuple(row.sampleId, file(row.read1), file(row.read2)) }
    .set { samples_ch }

//Run bwa v0.7.15 mem
//outbwa = file("$root/bwa_docker/")
//outbwa.mkdir()

process run_bwamem {
    publishDir path: "${root}/bwa_docker_v2/", mode: 'copy'
    
    input:
    set sampleId, file(read1), file(read2) from samples_ch
    
    output:
    val sampleId into sampleidpass
    file "*.sam" into bwa_sam
    
    script:
    """
    bwa mem -t 2 \
            -M \
            -R "@RG\tID:${sampleId}.Seq1\tCN:\tLB:${sampleId}_NoIndex\tPL:ILLUMINA\tPU:${sampleId}_NoIndex_L003\tSM:${sampleId}" \
            $ref \
            $read1 \
            $read2 \
            > ${sampleId}_NoIndex_L003_output.sam
    """
}

process run_sam_to_bam {
    publishDir path: "${root}/bwa_docker_v2/", mode: 'copy'

    input:
    val sampleId from sampleidpass
    file x from bwa_sam

    output:
    val sampleId into sampleidpassV2
    file "*.bam" into bwa_raw_bam
    
    script:
    """
    samtools view -@ 4 \
                  -Sb \
                  $x \
                  > ${sampleId}_NoIndex_L003_output.bam
            
    """
}


picard_tmpdir=file("${root}/bwa_docker_v2/picard_tmpdir/")
picard_tmpdir.mkdir()

process run_sort_sam {
    publishDir path: "${root}/bwa_docker_v2/", mode: 'copy'
    
    input:
    val sampleId from sampleidpassV2
    file x from bwa_raw_bam

    output:
    val sampleId into sampleidpassV3
    file "*.bam" into bwa_sorted_bam
    
    script:
    """
    java -Xmx2g -Djava.io.tmpdir=${picard_tmpdir} \
                -jar ${root}/picard-tools-1.130/picard.jar SortSam \
                VALIDATION_STRINGENCY=LENIENT \
                I=$x \
                O=${sampleId}_NoIndex_L003_output_sorted.bam \
                SORT_ORDER=coordinate
    """
}

process run_merge_sam_files {
    publishDir path: "${root}/bwa_docker_v2/", mode: 'copy'
    
    input:
    val sampleId from sampleidpassV3
    file x from bwa_sorted_bam

    output:
    val sampleId into sampleidpassV4
    file "*.bam" into bwa_noindex_merged_bam
    
    script:
    """
    java -Xmx2g -Djava.io.tmpdir=${picard_tmpdir} \
                -jar ${root}/picard-tools-1.130/picard.jar MergeSamFiles \
                USE_THREADING=true \
                VALIDATION_STRINGENCY=LENIENT \
                INPUT=$x \
                OUTPUT=TCEB1_RCC_${sampleId}_${sampleId}_NoIndex.merged.bam
    """
}

process run_markduplicates_buildbamindex {
    publishDir path: "${root}/bwa_docker_v2/", mode: 'copy'
    
    input:
    val sampleId from sampleidpassV4
    file x from bwa_noindex_merged_bam

    output:
    val sampleId into sampleidpassV5
    file "*.merged.bam" into bwa_merged_bam_forRTC
    file "*.merged.bam" into bwa_merged_bam_forIR
    file "*.merged.bai" into bwa_merged_bai
    file "*.merged.bam.metrics" into bwa_merged_bam_metrics
    file "*.merged.bam.bai" into bwa_buildbamindex
    
    script:
    """
    java -Xmx2g -Djava.io.tmpdir=${picard_tmpdir} \
                -jar ${root}/picard-tools-1.130/picard.jar MarkDuplicates \
                VALIDATION_STRINGENCY=LENIENT \
                INPUT=$x \
                OUTPUT=TCEB1_RCC_${sampleId}.merged.bam \
                METRICS_FILE=TCEB1_RCC_${sampleId}.merged.bam.metrics \
                PROGRAM_RECORD_ID=MarkDuplicates \
                CREATE_INDEX=true

    java -Xmx2g -Djava.io.tmpdir=${picard_tmpdir} \
                -jar ${root}/picard-tools-1.130/picard.jar BuildBamIndex \
                VALIDATION_STRINGENCY=LENIENT \
                INPUT=$x \
                OUTPUT=TCEB1_RCC_${sampleId}.merged.bam.bai
    """
}


//gatk_bundle_mills = "${root}/gatk_bundle/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
//gatk_bundle_known_indels = "${root}/gatk_bundle/Homo_sapiens_assembly38.known_indels.vcf.gz"
gatk_tmpdir=file("${root}/bwa_docker_v2/gatk_tmpdir/")
gatk_tmpdir.mkdir()

process run_realigned_target_creator {
    publishDir path: "${root}/output_gatk/GATK/3.7.0/TCEB1_RCC/phaseI_3.7.0/interval/", mode: 'copy'
    
    input:
    val sampleId from sampleidpassV5
    file x from bwa_merged_bam_forRTC

    output:
    val sampleId into sampleidpassV7
    file "TCEB1_RCC.chr1.intervals" into chr1_intervals
    
    script:
    """
    java -Xmx6g -Djava.io.tmpdir=${gatk_tmpdir} \
                -jar /src/GATK-3.7-0-gcfedb67/GenomeAnalysisTK.jar \
                -T RealignerTargetCreator \
                -I $x \
                -R /src/GATK-3.7-0-gcfedb67/genome.fa \
                -known /src/GATK-3.7-0-gcfedb67/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
                -known /src/GATK-3.7-0-gcfedb67/Homo_sapiens_assembly38.known_indels.vcf.gz \
                -L chr1 \
                -o TCEB1_RCC.chr1.intervals \
                --allow_potentially_misencoded_quality_scores \
                -nt 2
    """
}

process run_indel_realigner {
    publishDir path: "${root}/bwa_docker_v2/", mode: 'copy'
    
    input:
    val sampleId from sampleidpassV7
    file x from bwa_merged_bam_forIR
    file "TCEB1_RCC.chr1.intervals" from chr1_intervals

    output:
    val sampleId into sampleidpassV8
    file "TCEB1_RCC.indelrealigned.chr1.bam" into indelrealigned_chr1_bam_forbqsr
    file "TCEB1_RCC.indelrealigned.chr1.bam" into indelrealigned_chr1_bam_forprintreads
    file "TCEB1_RCC.indelrealigned.chr1.bai" into indelrealigned_chr1_bai
    
    script:
    """
    java -Xmx6g -Djava.io.tmpdir=${gatk_tmpdir} \
                -jar /src/GATK-3.7-0-gcfedb67/GenomeAnalysisTK.jar \
                -T IndelRealigner \
                -I $x \
                -R /src/GATK-3.7-0-gcfedb67/genome.fa \
                -compress 0 \
                -known /src/GATK-3.7-0-gcfedb67/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
                -known /src/GATK-3.7-0-gcfedb67/Homo_sapiens_assembly38.known_indels.vcf.gz \
                --allow_potentially_misencoded_quality_scores \
                -targetIntervals TCEB1_RCC.chr1.intervals \
                -o TCEB1_RCC.indelrealigned.chr1.bam \
                -L chr1 \
                -L unmapped
    """
}


//dbsnp138 = "${root}/gatk_bundle/resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf.gz"
//    dbsnp138_vcf_gz = Channel.fromPath(dbsnp138)
//dbsnp138_tbi = "${root}/gatk_bundle/resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf.gz.tbi"
//    dbsnp138_vcf_gz_tbi = Channel.fromPath(dbsnp138_tbi)

params.dbsnp = Channel.fromPath("${root}/gatk_bundle/resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf.gz").getVal()
params.dbsnptbi = Channel.fromPath("${root}/gatk_bundle/resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf.gz.tbi").getVal()

process run_bqsr {
    publishDir path: "${root}/output_gatk/GATK/3.7.0/TCEB1_RCC/phaseI_3.7.0/recalibration_table/", mode: 'copy'
    
    input:
    val sampleId from sampleidpassV8
    file "TCEB1_RCC.indelrealigned.chr1.bam" from indelrealigned_chr1_bam_forbqsr
    tuple file(dbsnp), file(dbsnptbi) from Channel.value([params.dbsnp, params.dbsnptbi])

    output:
    val sampleId into sampleidpassV9
    file "TCEB1_RCC_recalibration_table.grp" into recalibration_table_grp
    
    script:
    """
    java -Xmx8g -Djava.io.tmpdir=${gatk_tmpdir} \
                -jar /src/GATK-3.7-0-gcfedb67/GenomeAnalysisTK.jar \
                -T BaseRecalibrator \
                -I TCEB1_RCC.indelrealigned.chr1.bam \
                -R /src/GATK-3.7-0-gcfedb67/genome.fa \
                -l INFO \
                -nct 1 \
                -knownSites /src/GATK-3.7-0-gcfedb67/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
                -knownSites /src/GATK-3.7-0-gcfedb67/Homo_sapiens_assembly38.known_indels.vcf.gz \
                -knownSites $dbsnp \
                --allow_potentially_misencoded_quality_scores \
                --solid_nocall_strategy PURGE_READ \
                --solid_recal_mode SET_Q_ZERO_BASE_N \
                --covariate ReadGroupCovariate \
                --covariate QualityScoreCovariate \
                --covariate CycleCovariate \
                --covariate ContextCovariate \
                -o TCEB1_RCC_recalibration_table.grp
    """
}

process run_print_reads {
    publishDir path: "${root}/bwa_docker_v2/", mode: 'copy'
    
    input:
    val sampleId from sampleidpassV9
    file "TCEB1_RCC.indelrealigned.chr1.bam" from indelrealigned_chr1_bam_forprintreads
    file "TCEB1_RCC_recalibration_table.grp" from recalibration_table_grp

    output:
    val sampleId into sampleidpassV10
    file "*.recalibrated.chr1.bam" into recalibrated_chr1_bam
    file "*.recalibrated.chr1.bai" into recalibrated_chr1_bai
    
    script:
    """
    java -Xmx6g -Djava.io.tmpdir=${gatk_tmpdir} \
                -jar /src/GATK-3.7-0-gcfedb67/GenomeAnalysisTK.jar \
                -T PrintReads \
                -I TCEB1_RCC.indelrealigned.chr1.bam \
                -R /src/GATK-3.7-0-gcfedb67/genome.fa \
                -BQSR TCEB1_RCC_recalibration_table.grp \
                -nct 1 \
                -o ${sampleId}.recalibrated.chr1.bam \
                --sample_name $sampleId \
                -L chr1 \
                -L unmapped
    """
}

process run_reheader {
    publishDir path: "${root}/bwa_docker_v2/", mode: 'copy'
    
    input:
    val sampleId from sampleidpassV10
    file x from recalibrated_chr1_bam

    output:
    val sampleId into sampleidpassV11
    file "*.recalibrated.chr1.bam_reheadered.bam" into recalibrated_chr1_bam_reheadered_bam_forindex
    
    script:
    """
    samtools view -H $x | \
    samtools reheader - $x \
    > ${sampleId}.recalibrated.chr1.bam_reheadered.bam
    """
}

process run_build_bam_index {
    publishDir path: "${root}/bwa_docker_v2/", mode: 'copy'
    
    input:
    val sampleId from sampleidpassV11
    file x from recalibrated_chr1_bam_reheadered_bam_forindex

    output:
    val sampleId into sampleidpassV12
    tuple file(x), file("*.recalibrated.chr1.bam_reheadered.bam.bai") into recalibrated_chr1_bam_reheadered_bam_bai
    
    script:
    """
    java -Xmx8g -Djava.io.tmpdir=${picard_tmpdir} \
                -jar ${root}/picard-tools-1.130/picard.jar BuildBamIndex \
                INPUT=$x \
                OUTPUT=${sampleId}.recalibrated.chr1.bam_reheadered.bam.bai \
                VALIDATION_STRINGENCY=LENIENT
    """
}

params.bam = Channel.fromPath("${root}/bwa_docker_v2/*.recalibrated.chr1.bam_reheadered.bam").getVal()
params.bam_bai = Channel.fromPath("${root}/bwa_docker_v2/*.recalibrated.chr1.bam_reheadered.bam.bai").getVal()

process run_gatk_hc {
    publishDir path: "${root}/output_gatk/GATK/3.7.0/TCEB1_RCC/phaseII_3.7.0/raw_vcf/", mode: 'copy'
    
    input:
    val sampleId from sampleidpassV12
    tuple file(dbsnp), file(dbsnptbi), file(bam), file(bam_bai) from Channel.value([params.dbsnp, params.dbsnptbi, params.bam, params.bam_bai])
 //   tuple file(x), file("*.recalibrated.chr1.bam_reheadered.bai") from recalibrated_chr1_bam_reheadered_bam_bai
 //   tuple file(bam), file(bam_bai) from Channel.value([params.bam, params.bam_bai])

    output:
    val sampleId into sampleidpassV13
    file "TCEB1_RCC_chr1.vcf" into chr1_vcf
    
    script:
    """
    java -Xmx8g -Djava.io.tmpdir=${gatk_tmpdir} \
                 -jar /src/GATK-3.7-0-gcfedb67/GenomeAnalysisTK.jar \
                 -T HaplotypeCaller \
                 -R /src/GATK-3.7-0-gcfedb67/genome.fa \
                 -l INFO \
                 --output_mode EMIT_VARIANTS_ONLY \
                 --dbsnp $dbsnp \
                 -I $bam \
                 --sample_ploidy 2 \
                 -stand_call_conf 50 \
                 --allow_potentially_misencoded_quality_scores \
                 -nct 2 \
                 -L chr1 \
                 -o TCEB1_RCC_chr1.vcf
    """
}
