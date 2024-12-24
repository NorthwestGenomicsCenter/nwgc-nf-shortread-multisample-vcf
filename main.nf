include { BCFTOOLS_QUERY_SAMPLE_NAME } from "./modules/bcftools_query_sample_name.nf"
include { GATK_COMBINE_GVCFS } from "./modules/gatk_combine_gvcfs.nf"
include { GATK_CONCAT_CHROMOSOME_GVCFS } from "./modules/gatk_concat_chromosome_gvcfs.nf"
include { GATK_GENOTYPE_GVCFS } from "./modules/gatk_genotype_gvcfs.nf"
include { GATK_HARD_FILTER } from "./modules/gatk_hard_filter.nf"
include { PICARD_SORT_VCF } from "./modules/picard_sort_vcf.nf"


workflow {
    // Define needed variables from parameters
    def projectName = params.projectName
    def referenceGenome = params.referenceGenome
    def referenceDictionary = params.referenceDictionary
    def dbSnp = params.dbSnp
    def chrTag = params.isGRC38 ? "chr" : ""
    def chr25 = params.isGRC38 ? "M" : "MT"
    chromosomes = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, "X", "Y", chr25]
    ch_chromosomes = Channel.fromList(chromosomes)
    def genotypingTargetPrefix = chrTag
    def genotypingTargetSuffix = ""
    if (!params.sequencingTarget.equals('WHOLE_GENOME')) {
        genotypingTargetPrefix = "${params.sequencingTargetIntervalsDirectory}/${params.sequencingTarget}.chr"
        genotypingTargetSuffix = ".list"
    }
    def useCRDRExomeControls = params.useCRDRExomeControls
    def useCRDRWholeGenomeControls = params.useCRDRWholeGenomeControls
    def publishDirectory = params.publishDirectory

    // Create a list of all input gvcfs in the params file
    def gvcfs = []

    if (params.gvcfs) {
        gvcfs = gvcfs + params.gvcfs
    }

    if (params.sampleIds) {
        def sampleGvcfPaths = []
        for (sampleId: params.sampleIds) {
            sampleGvcfPaths.add("${params.sampleRootDirectory}/${sampleId}/${params.samplePipelineDirectoryName}/${sampleId}.${params.sequencingTarget}.main.g.vcf.gz")
        }
        
        gvcfs = gvcfs + sampleGvcfPaths
    }

    tbis = []
    for (int i=0; i < gvcfs.size(); i++) {
        tbis[i] = "${gvcfs[i]}.tbi"
    }

    // Make a channel of the gvcf and tbi files taken as input
    ch_gvcfs = Channel.fromList(gvcfs)
    ch_tbis = Channel.fromList(tbis)

    // gets a list of all input sample names
    BCFTOOLS_QUERY_SAMPLE_NAME(ch_gvcfs)
    ch_sampleNames = BCFTOOLS_QUERY_SAMPLE_NAME.out.collect()

    ch_genotypeTargetInput = Channel.empty()
    // If there are 200 or more GVCFs GATK GenotypeGVCFS wont work so we need to batch them first.
    // We also need to batch the GVCFs if using CRDR control samples as those are batched.
    int numSamples = gvcfs.size();
    int maxInputsForHaplotypeCaller = 200;
    int batchSize = 5;
    if (numSamples > 200 || params.useCRDRWholeGenomeControls || params.useCRDRExomeControls) {
        // Calculate needed number of batches / batch size
        int numBatches = 1
        if (numSamples > 200) {
            numBatches = (int) Math.ceil(numSamples/batchSize);
            while (numBatches > maxInputsForHaplotypeCaller) {
                batchSize += 5;
                numBatches = (int) Math.ceil(numSamples/batchSize);
            }
        }

        gvcfBatches = []
        for (int batchNum = 1; batchNum <= numBatches; batchNum++) {
            int gvcfStartId = (batchNum - 1) * batchSize;
            int gvcfEndId = gvcfStartId + batchSize;

            // This is in case the #gvcfs is not evenly divided among the batches
            if (batchNum == numBatches)
                gvcfEndId = numSamples;
            
            // Creates a tuple with the batch gvcfs, batch tbis, batch id, and each of the 25 chromsomes.
            for (def chr : chromosomes) {
                gvcfBatches.add([gvcfs.subList(gvcfStartId, gvcfEndId), tbis.subList(gvcfStartId, gvcfEndId), chr, batchNum])
            }
        }
        ch_gvcfBatches = Channel.fromList(gvcfBatches)

        GATK_COMBINE_GVCFS(ch_gvcfBatches, chrTag, projectName, referenceGenome)
        // Override gvcf and tbi channels with the batch gvcfs / tbis (by chromosome)
        ch_gvcfs = GATK_COMBINE_GVCFS.out.gvcf.groupTuple()
        ch_tbis = GATK_COMBINE_GVCFS.out.tbi.groupTuple()

        if (useCRDRWholeGenomeControls || useCRDRExomeControls) {
            // crdr controls
            ch_crdrGvcfControls = ch_chromosomes.map { def chr -> 
                def gvcfControls = []
                for (int crdrBatchId = 1; crdrBatchId <= 10; crdrBatchId++) {
                    if (useCRDRWholeGenomeControls) {
                        gvcfControls.add("/net/nwgc/vol1/references/human/grc38/multisample_controls/whole_genome/crdr_control_genomes.merged.${crdrBatchId}.chr${chr}.gvcf.gz")
                    }
                    else if (useCRDRExomeControls) {
                        gvcfControls.add("/net/nwgc/vol1/references/human/grc38/multisample_controls/exome/crdr_control_exomes.merged.${crdrBatchId}.chr${chr}.gvcf.gz")
                    }
                }
                return [chr, gvcfControls]
            }

            ch_crdrTbiControls = ch_chromosomes.map { def chr -> 
                def tbiControls = []
                for (int crdrBatchId = 1; crdrBatchId <= 10; crdrBatchId++) {
                    if (useCRDRWholeGenomeControls) {
                        tbiControls.add("/net/nwgc/vol1/references/human/grc38/multisample_controls/whole_genome/crdr_control_genomes.merged.${crdrBatchId}.chr${chr}.gvcf.gz.tbi")
                    }
                    else if (useCRDRExomeControls) {
                        tbiControls.add("/net/nwgc/vol1/references/human/grc38/multisample_controls/exome/crdr_control_exomes.merged.${crdrBatchId}.chr${chr}.gvcf.gz.tbi")
                    }
                }
                return [chr, tbiControls]
            }

            ch_gvcfs
            | combine(ch_crdrGvcfControls, by: 0)
            | map { def chr, def sampleGvcfs, def controls -> [chr, sampleGvcfs + controls] }
            | set { ch_gvcfs }

            ch_tbis
            | combine(ch_crdrTbiControls, by: 0)
            | map { def chr, def sampleTbis, def controls -> [chr, sampleTbis + controls] }
            | set { ch_tbis }
        }


        ch_gvcfs
        | combine(ch_tbis, by: 0) // combines channels so that now we have [chromosome, [gvcfs], [tbis]]
        | map { def chromosome, def combinedGvcfs, def combinedTbis -> [combinedGvcfs, combinedTbis, chromosome, "${genotypingTargetPrefix}${chromosome}${genotypingTargetSuffix}"]}
        | set { ch_genotypeTargetInput }
    }
    else {
        ch_genotypeTargetInput = ch_chromosomes.map { def chr -> [gvcfs, tbis, chr, "${genotypingTargetPrefix}${chromosome}${genotypingTargetSuffix}"] } 
    }

    
    GATK_GENOTYPE_GVCFS(ch_genotypeTargetInput, projectName, referenceGenome, dbSnp)
    GATK_CONCAT_CHROMOSOME_GVCFS(GATK_GENOTYPE_GVCFS.out.collect(), projectName, referenceGenome)
    PICARD_SORT_VCF(GATK_CONCAT_CHROMOSOME_GVCFS.out, projectName, referenceDictionary)
    GATK_HARD_FILTER(PICARD_SORT_VCF.out.collect(), projectName, referenceGenome, ch_sampleNames, publishDirectory)
}
