process GATK_HARD_FILTER {

    publishDir "${publishDirectory}", mode: "link", pattern: "${finalVcf}.gz"
    publishDir "${publishDirectory}", mode: "link", pattern: "${finalVcf}.gz.tbi"
    publishDir "${publishDirectory}", mode: "link", pattern: "${finalVcf}.gz.md5sum"
    publishDir "${publishDirectory}", mode: "link", pattern: "${finalVcf}.gz.tbi.md5sum"

    input:
        path inputGvcf
        val projectName
        val referenceGenome
        val sampleIds
        val publishDirectory

    output:
        tuple path("${finalVcf}.gz"), path("${finalVcf}.gz.tbi")
        path "${finalVcf}.gz.md5sum"
        path "${finalVcf}.gz.tbi.md5sum"

    script:
        javaOpts ='-Xmx14g'
        sansCohortVcf = "${projectName}.genotyped.sans_cohort.vcf"
        finalVcf = "${projectName}.HF.final.vcf"

        String sampleIdsString = ""
        for (int i = 0; i < sampleIds.size() - 1; i++) {
            sampleIdsString += sampleIds[i] + ","
        }
        sampleIdsString += sampleIds[sampleIds.size() - 1]

        """
        sed '/GATKCommandLine=/d' ${inputGvcf} | bcftools view -c 1 -s ${sampleIdsString} > ${sansCohortVcf}

        # Hard Filter SNVs
        time java -Xmx10g \
            -jar \$MOD_GSGATK_DIR/GenomeAnalysisTK.jar \
            -T VariantFiltration \
            -R ${referenceGenome} \
            --filterName QDFilter -filter "QD < 5.0" \
            --filterName QUALFilter -filter "QUAL <= 50.0" \
            --filterName ABFilter -filter "ABHet > 0.75" \
            --filterName SBFilter -filter "SB >= 0.10" \
            --filterName HRunFilter -filter "HRun > 4.0" \
            --clusterSize 3 \
            --clusterWindowSize 10 \
            -V:vcf ${sansCohortVcf} \
            -o ${finalVcf}

        ## BGZIP the vcf
        bgzip -f ${finalVcf}
        tabix -p vcf -f ${finalVcf}.gz

        md5sum ${finalVcf}.gz | awk '{print \$1}' > ${finalVcf}.gz.md5sum
        md5sum ${finalVcf}.gz.tbi | awk '{print \$1}' > ${finalVcf}.gz.tbi.md5sum
        """
}