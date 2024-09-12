process GATK_COMBINE_GVCFS {

    input:
        tuple path(inputGvcfs), path(inputTbis), val(chromosome), val(batchId)
        val chrTag
        val projectName
        val referenceGenome

    output:
        tuple val(chromosome), path("${vcfFile}.gz"), emit: gvcf
        tuple val(chromosome), path("${vcfFile}.gz.tbi"), emit: tbi

    script:
        javaOpts = '-Xmx14g'
        vcfFile = "${projectName}.combined.${batchId}.${chromosome}.g.vcf"

        String inputGvcfsString = ""
        for (gvcfFile: inputGvcfs) {
            inputGvcfsString += "-V " + gvcfFile + " \\\n"
        }

        """
        java ${javaOpts} \
            -jar \$MOD_GSGATK_DIR/GenomeAnalysisTK.jar \
            -T CombineGVCFs \
            -L ${chrTag}${chromosome} \
            -R ${referenceGenome} \
            ${inputGvcfsString} \
            -o ${vcfFile}

        ## Zip and tabix the vcf file        
        bgzip -f ${vcfFile}
        tabix -p vcf -f ${vcfFile}.gz
        """
}