process GATK_CONCAT_CHROMOSOME_GVCFS {

    input:
        path inputGvcfs
        val projectName
        val referenceGenome

    output:
        path gvcfFile

    script:
        javaOpts = '-Xmx8g'
        gvcfFile = "${projectName}.genotyped.vcf"

        String inputGvcfsString = ""
        for (gvcfFile: inputGvcfs) {
            inputGvcfsString += "-V " + gvcfFile + " \\\n"
        }

        """
        java ${javaOpts} \
            -cp \$MOD_GSGATK_DIR/GenomeAnalysisTK.jar org.broadinstitute.gatk.tools.CatVariants \
            -R ${referenceGenome} \
            --assumeSorted \
            ${inputGvcfsString} \
            -out ${gvcfFile}

        """
}