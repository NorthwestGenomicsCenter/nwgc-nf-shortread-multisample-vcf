process PICARD_SORT_VCF {

    input:
        path gvcf
        val projectName
        val referenceDictionary

    output:
        path gvcfName

    script:
        gvcfName = "${projectName}.sorted.genotyped.vcf"

        """
        java \
            -XX:InitialRAMPercentage=80 \
            -XX:MaxRAMPercentage=85 \
            -jar \$PICARD_DIR/picard.jar SortVcf \
            I=${gvcf} \
            O=${gvcfName} \
            SEQUENCE_DICTIONARY=${referenceDictionary}
        """
}