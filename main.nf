#!/usr/bin/env nextflow

process DRY_RUN {

    input:
    tuple path(ref), path(ref_fai)
    tuple val(family_id), val(sample_name_child), path(reads_child), path(reads_bai_child), val(sample_name_parent1), path(reads_parent1), path(reads_bai_parent1), val(sample_name_parent2), path(reads_parent2), path(reads_bai_parent2), val(regions), val(model_type)

    output:
    tuple path(ref), path(ref_fai)
    tuple val(family_id), val(sample_name_child), path(reads_child), path(reads_bai_child), val(sample_name_parent1), path(reads_parent1), path(reads_bai_parent1), val(sample_name_parent2), path(reads_parent2), path(reads_bai_parent2), val(regions), env(make_examples_cs_args), env(make_examples_calling_args), env(call_variants_child_args), env(call_variants_parent1_args), env(call_variants_parent2_args)
    
    script:
    """
    run_deeptrio \
        --model_type=$model_type \
        --ref=$ref \
        --sample_name_child=$sample_name_child \
        --sample_name_parent1=$sample_name_parent1 \
        --sample_name_parent2=$sample_name_parent2 \
        --reads_child=$reads_child \
        --reads_parent1=$reads_parent1 \
        --reads_parent2=$reads_parent2 \
        --output_vcf_child=child.vcf.gz \
        --output_vcf_parent1=parent1.vcf.gz \
        --output_vcf_parent2=parent2.vcf.gz \
        --output_gvcf_child=child.g.vcf.gz \
        --output_gvcf_parent1=parent1.g.vcf.gz \
        --output_gvcf_parent2=parent2.g.vcf.gz \
        --dry_run=true > commands.txt

    make_examples_cs_args=\$(grep "/opt/deepvariant/bin/deeptrio/make_examples --mode candidate_sweep" commands.txt | awk -F'/opt/deepvariant/bin/deeptrio/make_examples' '{print \$2}' | sed 's/--ref "[^"]*"//g' | sed 's/--sample_name "[^"]*"//g' | sed 's/--reads "[^"]*"//g' | sed 's/--sample_name_parent1 "[^"]*"//g' | sed 's/--reads_parent1 "[^"]*"//g' | sed 's/--sample_name_parent2 "[^"]*"//g' | sed 's/--reads_parent2 "[^"]*"//g' | sed 's/--examples "[^"]*"//g' | sed 's/--candidate_positions "[^"]*"//g' | sed 's/--gvcf "[^"]*"//g')
    make_examples_calling_args=\$(grep "/opt/deepvariant/bin/deeptrio/make_examples --mode calling" commands.txt | awk -F'/opt/deepvariant/bin/deeptrio/make_examples' '{print \$2}' | sed 's/--ref "[^"]*"//g' | sed 's/--sample_name "[^"]*"//g' | sed 's/--reads "[^"]*"//g' | sed 's/--sample_name_parent1 "[^"]*"//g' | sed 's/--reads_parent1 "[^"]*"//g' | sed 's/--sample_name_parent2 "[^"]*"//g' | sed 's/--reads_parent2 "[^"]*"//g' | sed 's/--examples "[^"]*"//g' | sed 's/--candidate_positions "[^"]*"//g' | sed 's/--gvcf "[^"]*"//g')
    call_variants_child_args=\$(grep "/opt/deepvariant/bin/call_variants" commands.txt | grep "child" | awk -F'/opt/deepvariant/bin/call_variants' '{print \$2}' | sed 's/--outfile "[^"]*"//g' | sed 's/--examples "[^"]*"//g')
    call_variants_parent1_args=\$(grep "/opt/deepvariant/bin/call_variants" commands.txt | grep "parent1" | awk -F'/opt/deepvariant/bin/call_variants' '{print \$2}' | sed 's/--outfile "[^"]*"//g' | sed 's/--examples "[^"]*"//g')
    call_variants_parent2_args=\$(grep "/opt/deepvariant/bin/call_variants" commands.txt | grep "parent2" | awk -F'/opt/deepvariant/bin/call_variants' '{print \$2}' | sed 's/--outfile "[^"]*"//g' | sed 's/--examples "[^"]*"//g')
    """

}

process MAKE_EXAMPLES {

    input:
    tuple path(ref), path(ref_fai)
    tuple val(family_id), val(sample_name_child), path(reads_child), path(reads_bai_child), val(sample_name_parent1), path(reads_parent1), path(reads_bai_parent1), val(sample_name_parent2), path(reads_parent2), path(reads_bai_parent2), val(regions), val(make_examples_cs_args), val(make_examples_calling_args), val(call_variants_child_args), val(call_variants_parent1_args), val(call_variants_parent2_args)

    output:
    tuple val(family_id), val(sample_name_child), path('make_examples_child.*.gz'), path('gvcf_child.*.gz'), path('*.example_info.json'), val(call_variants_child_args)         , emit: child
    tuple val(family_id), val(sample_name_parent1), path('make_examples_parent1.*.gz'), path('gvcf_parent1.*.gz'), path('*.example_info.json'), val(call_variants_parent1_args) , emit: parent1
    tuple val(family_id), val(sample_name_parent2), path('make_examples_parent2.*.gz'), path('gvcf_parent2.*.gz'), path('*.example_info.json'), val(call_variants_parent2_args) , emit: parent2
    
    script:
    def regions_arg = regions ? "--regions ${regions}" : ""
    """
    seq 0 ${task.cpus - 1} | parallel -q --halt 2 --line-buffer make_examples \\
        --ref "${ref}" --sample_name "${sample_name_child}" --reads "${reads_child}" --sample_name_parent1 "${sample_name_parent1}" --reads_parent1 "${reads_parent1}" \\
        --sample_name_parent2 "${sample_name_parent2}" --reads_parent2 "${reads_parent2}" ${regions_arg} --examples "make_examples.tfrecord@${task.cpus}.gz" --gvcf "gvcf.tfrecord@${task.cpus}.gz" --candidate_positions "candidate_positions@${task.cpus}.gz" ${make_examples_cs_args}
    
    seq 0 ${task.cpus - 1} | parallel -q --halt 2 --line-buffer make_examples \\
        --ref "${ref}" --sample_name "${sample_name_child}" --reads "${reads_child}" --sample_name_parent1 "${sample_name_parent1}" --reads_parent1 "${reads_parent1}" \\
        --sample_name_parent2 "${sample_name_parent2}" --reads_parent2 "${reads_parent2}" ${regions_arg} --examples "make_examples.tfrecord@${task.cpus}.gz" --gvcf "gvcf.tfrecord@${task.cpus}.gz" --candidate_positions "candidate_positions@${task.cpus}.gz" ${make_examples_calling_args}
    """

}

process CALL_VARIANTS {

    input:
    tuple val(family_id), val(sample_name), path(make_examples), val(gvcf), path(example_info), val(call_variants_args)

    output:
    tuple val(family_id), val(sample_name), path('*.gz'), val(gvcf)

    script:
    def matcher = make_examples[0].baseName =~ /^(.+)-\d{5}-of-(\d{5})$/
    def make_examples_name = matcher[0][1]
    def make_examples_num_shards = matcher[0][2] as int
    """
    call_variants --outfile "call_variants_output.tfrecord.gz" --examples "${make_examples_name}@${make_examples_num_shards}.gz" ${call_variants_args}
    """

}

process POSTPROCESS_VARIANTS {

    publishDir "${params.output_dir}/${family_id}/${sample_name}" , mode: 'copy'

    input:
    tuple path(ref), path(ref_fai)
    tuple val(family_id), val(sample_name), path(call_variants), path(gvcf)

    output:
    path("${sample_name}.*")

    script:
    def matcher = gvcf[0].baseName =~ /^(.+)-\d{5}-of-(\d{5})$/
    def gvcf_name = matcher[0][1]
    def gvcf_num_shards = matcher[0][2] as int
    """
    postprocess_variants --ref "${ref}" --sample_name "${sample_name}" --infile "call_variants_output.tfrecord.gz" --nonvariant_site_tfrecord_path "${gvcf_name}@${gvcf_num_shards}.gz" --cpus "${task.cpus}" --outfile "${sample_name}.vcf.gz" --gvcf_outfile "${sample_name}.g.vcf.gz"
    """

}


workflow {

    // Create a tuple of the reference FASTA and its index file. Throw an error if the index file is not found.
    ref_fai_path = file("${params.ref}.fai")
    if (!ref_fai_path.exists()) {
        throw new RuntimeException("Reference FASTA index file not found, expected at: ${ref_fai_path}")
    }
    ch_ref = [ file(params.ref), ref_fai_path ]

    // Create a channel that holds tuples for each sample, containing the sample name, the BAM file and its index file. Throw an error if the index file is not found.
    ch_samples = Channel.of(params.samples).flatMap()
    ch_samples = ch_samples.map { family_id, sample_name_child, reads_child, sample_name_parent1, reads_parent1, sample_name_parent2, reads_parent2, regions, model_type ->
        def reads_bai_child = file("${reads_child}.bai")
        if (!reads_bai_child.exists()) {
            throw new RuntimeException("BAM index file not found for sample: ${sample_name_child}, expected at: ${reads_bai_child}")
        }
        def reads_bai_parent1 = file("${reads_parent1}.bai")
        if (!reads_bai_parent1.exists()) {
            throw new RuntimeException("BAM index file not found for sample: ${sample_name_parent1}, expected at: ${reads_bai_parent1}")
        }
        def reads_bai_parent2 = file("${reads_parent2}.bai")
        if (!reads_bai_parent2.exists()) {
            throw new RuntimeException("BAM index file not found for sample: ${sample_name_parent2}, expected at: ${reads_bai_parent2}")
        }
        
        def regions_val = regions == '' ? [] : regions
        [ family_id, sample_name_child, file(reads_child), reads_bai_child, sample_name_parent1, file(reads_parent1), reads_bai_parent1, sample_name_parent2, file(reads_parent2), reads_bai_parent2, regions_val, model_type ]
    }

    // Do a dry run of DeepTrio to extract the arguments for MAKE_EXAMPLES and CALL_VARIANTS stages
    DRY_RUN(ch_ref, ch_samples)

    // Run the MAKE_EXAMPLES stage
    MAKE_EXAMPLES(DRY_RUN.out) 

    // Create a channel for the outputs of MAKE_EXAMPLES stage (child's, parent1's and parent2's) and run the CALL_VARIANTS stage
    ch_call_variants = MAKE_EXAMPLES.out.child.mix(MAKE_EXAMPLES.out.parent1, MAKE_EXAMPLES.out.parent2)    
    CALL_VARIANTS(ch_call_variants) 

    // Run the combined POSTPROCESS_VARIANTS and VCF_STATS_REPORT stages with the outputs of CALL_VARIANTS stage
    POSTPROCESS_VARIANTS(ch_ref, CALL_VARIANTS.out)

}