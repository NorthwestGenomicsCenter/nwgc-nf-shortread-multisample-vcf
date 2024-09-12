process BCFTOOLS_QUERY_SAMPLE_NAME {

    input:
        path gvcf
    
    output:
        env SAMPLE_NAME

    script:
        """
        SAMPLE_NAME=\$(bcftools query -l ${gvcf})
        """
}