process GATK_GENOTYPE_GVCFS {

    input:
        tuple path(inputGvcfs), path(inputTbis), val(chromosome), val(genotypingTarget)
        val projectName
        val referenceGenome
        val dbSnp

    output:
        path vcfFile

    script:
        javaOpts = '-Xmx46g'
        vcfFile = "${projectName}.genotyped.${chromosome}.vcf"

        String inputGvcfsString = ""
        for (gvcfFile: inputGvcfs) {
            inputGvcfsString += "-V " + gvcfFile + " \\\n"
        }

        """
        java ${javaOpts} \
            -jar \$MOD_GSGATK_DIR/GenomeAnalysisTK.jar \
            -T GenotypeGVCFs \
            -R ${referenceGenome} \
            -L ${genotypingTarget} \
            --dbsnp ${dbSnp} \
            ${inputGvcfsString} \
            -o ${vcfFile}
        """
}