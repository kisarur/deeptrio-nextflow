# Nextflow Pipeline for DeepTrio

This repository contains a Nextflow pipeline for Google’s DeepTrio, optimised for execution on NCI Gadi.

## Quickstart Guide

1. Edit the `pipeline_params.yml` file to include:
    - `samples`:  a list of family trio entries, where each entry includes a unique family id, child sample name, child BAM file path, parent1 sample name, parent1 BAM file path, parent2 sample name, parent2 BAM path, path to an optional regions-of-interest BED file (set to `''` if not required), and the model type. For BAM files, ensure corresponding .bai is in the same directory.
    - `ref`: path to the reference FASTA (ensure corresponding .fai is in the same directory).
    - `output_dir`: directory path to save output files.
    - `nci_project`, `nci_storage` : NCI project and storage.

2. Update `nextflow.config` to match the resource requirements for each stage of the pipeline. For NCI Gadi, you may need to adjust only `time` and `disk` (i.e. jobfs) parameters based on the size of the datasets used (the default values are tested to be suitable for a PacBio dataset with each family member's BAM being ~45GB in size).

3. Load the Nextflow module and run the pipeline using the following commands:
    ```bash
    module load nextflow/24.04.1
    nextflow run main.nf -params-file pipeline_params.yml
    ```

    Note: Additional Nextflow options can be included (e.g., `-resume` to resume from a previously paused/interrupted run)

4. For each family trio, output files will be stored in the directory `output_dir/family_id`.

## Notes  

1. It is assumed that the user has access to NCI's `if89` project (required for using DeepTrio via `module load`). If not, simply request access using this [link](https://my.nci.org.au/mancini/project/if89).

## Case Study

A case study was conducted using a PacBio dataset with each family member's BAM being ~45GB in size to evaluate the runtime and service unit (SU) efficiency of *deeptrio-nextflow* compared to the original DeepTrio running on a single node. The benchmarking results are summarised in the table below.

<table>
    <thead>
        <tr>
            <th>Version</th>
            <th>Gadi Resources</th>
            <th>Runtime (hh:mm:ss)</th>
            <th>SUs</th>
        </tr>
    </thead>
    <tbody>
        <tr>
            <td>Original DeepTrio</td>
            <td><code>gpuvolta</code> (48 CPUs, 4 GPUs, 384 GB memory)</td>
            <td>13:27:29</td>
            <td>1937.96</td>
        </tr>
        <tr>
            <td><i>deeptrio-nextflow</i></td>
            <td><code>normalsr</code> (104 CPUs, 500 GB memory) → <code>gpuvolta</code> (12 CPUs, 1 GPU, 96 GB memory) → <code>normalbw</code> (14 CPUs, 128 GB memory) </td>
            <td>05:22:16</td>
            <td>825.17</td>
        </tr>
    </tbody>
</table>

## Acknowledgments

The *deeptrio-nextflow* workflow was developed by Dr Kisaru Liyanage and Dr Matthew Downton (National Computational Infrastructure), with support from Australian BioCommons as part of the Workflow Commons project.

We thank Leah Kemp (Garvan Institute of Medical Research) for her collaboration in providing test datasets and assisting with pipeline testing.