nextflow.enable.dsl = 2

// Group reads by integration site across all timepoints for each patient.
// PCR duplicates share sample_id and PCR_ID.

process PREPARE_PATIENT {
    tag "${patient_id}"
    container 'python:3.12-slim-bookworm'
    cpus 1
    memory '3 GB'

    publishDir "${params.outdir}/merged_fastas", mode: 'copy', pattern: '*_merged.fasta', overwrite: true
    publishDir "${params.outdir}/unintegrated_sequences", mode: 'copy', pattern: '*_unintegrated.fasta', overwrite: true
    publishDir "${params.outdir}/review/unmatched_fastas", mode: 'copy', pattern: '*.unmatched.fasta', overwrite: true

    input:
    tuple val(patient_id), val(sample_ids), val(timepoints), val(original_ref_paths), path(fastas, stageAs: 'fastas??????/*'), path(site_files, stageAs: 'sites??????/*'), path(ref_files, stageAs: 'refs??????/*')

    output:
    path "${patient_id}_merged.fasta", emit: merged_fasta
    tuple val(patient_id), path('clone_inputs'), path('patient_reference.fasta'), path('clone_info.json'), path("${patient_id}.counts.csv"), path("${patient_id}.unmatched.csv"), path("${patient_id}.assignments.csv"), path("${patient_id}.unused_sites.csv"), path("${patient_id}.reference_selection.csv"), path("${patient_id}_unintegrated.fasta"), emit: prepared

    script:
    def fs = fastas instanceof List ? fastas : [fastas]
    def ss = site_files instanceof List ? site_files : [site_files]
    def rs = ref_files instanceof List ? ref_files : [ref_files]
    def entries = (0..<sample_ids.size()).collect { i ->
        [sample_id:sample_ids[i], timepoint:timepoints[i],
         fasta:fs[i].toString(), integration_site:ss[i].toString(),
         reference:rs[i].toString(), original_best_hiv_ref_filepath:original_ref_paths[i]]
    }
    def contents = groovy.json.JsonOutput.toJson(entries)
    """
    set -euo pipefail
    cat > staged_samples.json <<'NF_SAMPLES_EOF'
    ${contents}
    NF_SAMPLES_EOF
    prepare_patient.py --patient '${patient_id}' --samples staged_samples.json --check-orientation ${params.check_orientation}
    """
}

process ALIGN_CLONES {
    tag "${patient_id}"
    container 'quay.io/biocontainers/mafft:7.505--hec16e2b_0'
    cpus 2
    memory '4 GB'
    time '8h'

    input:
    tuple val(patient_id), path(clone_inputs), path(reference), path(clone_info), path(counts), path(unmatched), path(assignments), path(unused_sites), path(reference_selection), path(unintegrated)

    output:
    tuple val(patient_id), path('aligned_clones'), path(clone_info), path(counts), path(unmatched), path(assignments), path(unused_sites), path(reference_selection), path(unintegrated), emit: aligned

    script:
    """
    set -euo pipefail
    mkdir -p aligned_clones
    for fa in clone_inputs/*.fasta; do
        [[ -f "\$fa" ]] || continue
        name=\$(basename "\$fa")
        n=\$(grep -c '^>' "\$fa")
        if [[ "\$n" -eq 1 ]]; then
            cp "\$fa" "aligned_clones/\$name"
        else
            mafft --auto --thread ${task.cpus} \\
                --addfragments "\$fa" '${reference}' > "aligned_clones/\$name"
        fi
    done
    """
}

process FINALIZE_PATIENT {
    tag "${patient_id}"
    container 'python:3.12-slim-bookworm'
    cpus 1
    memory '3 GB'

    publishDir "${params.outdir}/clone_consensus", mode: 'copy', pattern: '*_consensus.fasta', overwrite: true
    publishDir "${params.outdir}/clone_alignments", mode: 'copy', pattern: '*', saveAs: { name -> (name.endsWith('.fasta') || name.endsWith('.csv')) ? null : name }, overwrite: true

    input:
    tuple val(patient_id), path(aligned_clones), path(clone_info), path(counts), path(unmatched), path(assignments), path(unused_sites), path(reference_selection), path(unintegrated)

    output:
    path "${patient_id}_consensus.fasta", emit: final_fasta
    path "${patient_id}", emit: clone_alignments
    path "${patient_id}.clone_summary.csv", emit: clone_summaries
    path "${patient_id}.patient_stats.csv", emit: stats
    path "${patient_id}.counts.csv", emit: counts
    path "${patient_id}.unmatched.csv", emit: unmatched
    path "${patient_id}.assignments.csv", emit: assignments
    path "${patient_id}.unused_sites.csv", emit: unused_sites
    path "${patient_id}.reference_selection.csv", emit: references

    script:
    """
    set -euo pipefail
    finalize_patient.py \\
        --patient '${patient_id}' \\
        --aligned-clones '${aligned_clones}' \\
        --clone-info '${clone_info}' \\
        --counts '${counts}' \\
        --unmatched '${unmatched}' \\
        --unintegrated '${unintegrated}' \\
        --reference-selection '${reference_selection}' \\
        --min-base-fraction ${params.min_base_fraction} \\
        --min-family-fraction ${params.min_family_fraction} \\
        --min-families-per-call ${params.min_families_per_call}
    """
}

process MAKE_SUMMARIES {
    container 'python:3.12-slim-bookworm'
    cpus 1
    memory '1 GB'

    publishDir params.outdir, mode: 'copy', pattern: 'patient_sequence_counts.csv', overwrite: true
    publishDir params.outdir, mode: 'copy', pattern: 'clone_summary.csv', overwrite: true
    publishDir params.outdir, mode: 'copy', pattern: 'patient_clone_statistics.csv', overwrite: true
    publishDir params.outdir, mode: 'copy', pattern: 'reference_selection.csv', overwrite: true
    publishDir "${params.outdir}/review", mode: 'copy', pattern: 'unmatched_fasta_reads.csv', overwrite: true
    publishDir "${params.outdir}/review", mode: 'copy', pattern: 'read_assignments.csv', overwrite: true
    publishDir "${params.outdir}/review", mode: 'copy', pattern: 'unused_integration_records.csv', overwrite: true

    input:
    path counts
    path clone_summaries
    path stats
    path unmatched
    path assignments
    path unused_sites
    path references

    output:
    path 'patient_sequence_counts.csv'
    path 'clone_summary.csv'
    path 'patient_clone_statistics.csv'
    path 'reference_selection.csv'
    path 'unmatched_fasta_reads.csv'
    path 'read_assignments.csv'
    path 'unused_integration_records.csv'

    script:
    """
    set -euo pipefail
    make_summaries.py \\
      --counts '*.counts.csv' \\
      --clone-summaries '*.clone_summary.csv' \\
      --stats '*.patient_stats.csv' \\
      --unmatched '*.unmatched.csv' \\
      --assignments '*.assignments.csv' \\
      --unused-sites '*.unused_sites.csv' \\
      --references '*.reference_selection.csv'
    """
}

workflow {
    if (!params.input) error 'Required: --input /path/to/annotation.csv'
    if (!params.outdir) error 'Required: --outdir results'

    def manifest = file(params.input, checkIfExists: true)
    def manifestDir = manifest.parent

    // Read the sample manifest, then collect all samples for each patient.
    patient_files = Channel.fromPath(manifest.toString(), checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            def sid = (row.sample_id ?: row['\uFEFFsample_id'])?.toString()?.trim()
            def pid = row.patient_id?.toString()?.trim()
            def tp  = row.timepoint?.toString()?.trim()
            def fastaPath = row.hiv_fasta_filepath?.toString()?.trim()
            def sitePath  = row.integration_site_filepath?.toString()?.trim()
            def refPath   = row.best_hiv_ref_filepath?.toString()?.trim()
            if (![sid, pid, tp, fastaPath, sitePath, refPath].every { it }) {
                error 'Input CSV must specify sample_id, patient_id, timepoint, hiv_fasta_filepath, integration_site_filepath, best_hiv_ref_filepath'
            }
            [sid, pid].each { value ->
                if (!(value ==~ /[A-Za-z0-9][A-Za-z0-9._-]*/))
                    error "Unsafe sample_id/patient_id: ${value}"
            }
            def fasta = manifestDir.resolve(fastaPath).normalize()
            def site  = manifestDir.resolve(sitePath).normalize()
            def ref   = manifestDir.resolve(refPath).normalize()
            [fasta, ref].each { fp ->
                if (!(fp.toString().toLowerCase() ==~ /.*\.(fa|fasta|fna)(\.gz)?/))
                    error "Unexpected FASTA extension: ${fp}"
                if (!fp.toFile().isFile()) error "FASTA file not found: ${fp}"
            }
            if (!site.toFile().isFile()) error "Integration-site table not found: ${site}"
            tuple(pid, sid, tp, fasta, site, ref)
        }
        .groupTuple()
        .map { pid, sids, tps, fas, sites, refs ->
            def samples = (0..<sids.size()).collect { i ->
                [id:sids[i], timepoint:tps[i], fasta:fas[i], site:sites[i], reference:refs[i]]
            }.sort { a,b -> a.id <=> b.id }
            if (samples.collect { it.id }.unique().size() != samples.size()) {
                error "Duplicate sample_id for patient ${pid}"
            }
            tuple(pid, samples.collect { it.id }, samples.collect { it.timepoint },
                  samples.collect { it.reference.toString() }, samples.collect { it.fasta },
                  samples.collect { it.site }, samples.collect { it.reference })
        }

    PREPARE_PATIENT(patient_files)
    // Align each clone to the selected patient reference before calling consensus.
    ALIGN_CLONES(PREPARE_PATIENT.out.prepared)
    FINALIZE_PATIENT(ALIGN_CLONES.out.aligned)

    // Combine patient-level reports.
    MAKE_SUMMARIES(
        FINALIZE_PATIENT.out.counts.collect(),
        FINALIZE_PATIENT.out.clone_summaries.collect(),
        FINALIZE_PATIENT.out.stats.collect(),
        FINALIZE_PATIENT.out.unmatched.collect(),
        FINALIZE_PATIENT.out.assignments.collect(),
        FINALIZE_PATIENT.out.unused_sites.collect(),
        FINALIZE_PATIENT.out.references.collect()
    )
}
