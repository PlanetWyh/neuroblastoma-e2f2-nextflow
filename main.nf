nextflow.enable.dsl = 2

/*
 * Stage 2 bulk E2F2 / survival pipeline
 */

process BUILD_TARGET_MATRIX {

    publishDir "${params.outdir}/stage2/TARGET", mode: 'copy'

    input:
    path metadata
    val data_dir

    output:
    path "target_counts_matrix.csv", emit: counts
    path "target_logcpm_matrix.csv", emit: logcpm

    script:
    """
    python ${projectDir}/scripts/stage2/build_target_matrix.py \
        --metadata ${metadata} \
        --data_dir ${projectDir}/${data_dir} \
        --counts_out target_counts_matrix.csv \
        --logcpm_out target_logcpm_matrix.csv
    """
}


process CALCULATE_TARGET_E2F2_SSGSEA {

    publishDir "${params.outdir}/stage2/TARGET", mode: 'copy'

    input:
    path logcpm
    path geneset

    output:
    path "target_e2f2_ssgsea_scores.csv", emit: scores

    script:
    """
    python ${projectDir}/scripts/stage2/calculate_e2f2_ssgsea.py \
        --logcpm_file ${logcpm} \
        --geneset_file ${geneset} \
        --output_file target_e2f2_ssgsea_scores.csv
    """
}


process PREPARE_TARGET_SURVIVAL {

    publishDir "${params.outdir}/stage2/TARGET", mode: 'copy'

    input:
    path clinical
    path e2f2_scores
    path mycn_file

    output:
    path "target_survival_table.csv", emit: survival_table
    path "target_survival_summary.txt", emit: summary

    script:
    """
    python ${projectDir}/scripts/stage2/prepare_survival_data.py \
        --clinical_file ${clinical} \
        --e2f2_scores_file ${e2f2_scores} \
        --mycn_file ${mycn_file} \
        --output_file target_survival_table.csv \
        --summary_file target_survival_summary.txt
    """
}


process TARGET_SURVIVAL_ANALYSIS {

    publishDir "${params.outdir}/stage2/TARGET/figures", mode: 'copy'

    input:
    path survival_table

    output:
    path "km_e2f2_survival.png", emit: km_plot
    path "logrank_results.csv", emit: logrank

    script:
    """
    python ${projectDir}/scripts/stage2/survival_analysis.py \
        --input_file ${survival_table} \
        --output_dir .
    """
}


process VALIDATE_E2F2_GSE49711 {

    publishDir "${params.outdir}/stage2/GSE49711", mode: 'copy'

    input:
    path geneset

    output:
    path "GSE49711_clinical.csv", emit: clinical
    path "GSE49711_expression_log2.csv", emit: expression
    path "GSE49711_e2f2_ssgsea_scores.csv", emit: scores
    path "GSE49711_e2f2_clinical_merged.csv", emit: merged

    script:
    """
    python ${projectDir}/scripts/stage2/validate_e2f2_gse49711.py \
        --geneset_file ${geneset} \
        --output_dir .
    """
}


process PLOT_GSE49711 {

    publishDir "${params.outdir}/stage2/GSE49711/figures", mode: 'copy'

    input:
    path merged_csv

    output:
    path "*.png", emit: plots

    script:
    """
    python ${projectDir}/scripts/stage2/plot_gse49711_e2f2.py \
        --input_file ${merged_csv} \
        --output_dir .
    """
}

process LOAD_NB {

    publishDir "${params.outdir}/stage1", mode: 'copy'

    input:
    path input_dir

    output:
    path "nb_merged_raw.h5ad", emit: raw_h5ad

    script:
    """
    python ${projectDir}/scripts/stage1/load_nb.py \
        --input_dir ${input_dir} \
        --output_file nb_merged_raw.h5ad
    """
}


process PREPROCESS_NB {

    publishDir "${params.outdir}/stage1", mode: 'copy'

    input:
    path raw_h5ad

    output:
    path "nb_processed.h5ad", emit: processed_h5ad

    script:
    """
    python ${projectDir}/scripts/stage1/preprocess_nb.py \
        --input_file ${raw_h5ad} \
        --output_file nb_processed.h5ad
    """
}


process SUBSET_NE {

    publishDir "${params.outdir}/stage1", mode: 'copy'

    input:
    path processed_h5ad

    output:
    path "ne_cells.h5ad", emit: ne_h5ad

    script:
    """
    python ${projectDir}/scripts/stage1/subset_ne.py \
        --input_file ${processed_h5ad} \
        --output_file ne_cells.h5ad
    """
}


process CLEAN_NE {

    publishDir "${params.outdir}/stage1", mode: 'copy'

    input:
    path ne_h5ad

    output:
    path "ne_clean.h5ad", emit: ne_clean_h5ad

    script:
    """
    python ${projectDir}/scripts/stage1/clean_ne.py \
        --input_file ${ne_h5ad} \
        --output_file ne_clean.h5ad
    """
}


process EXPORT_FOR_PYSCENIC {

    publishDir "${params.outdir}/stage1", mode: 'copy'

    input:
    path ne_clean_h5ad

    output:
    path "ne_expr_matrix.csv", emit: expr_matrix
    path "ne_metadata.csv", emit: metadata

    script:
    """
    python ${projectDir}/scripts/stage1/export_for_pyscenic.py \
        --input_file ${ne_clean_h5ad} \
        --output_dir .
    """
}


process PYSCENIC_GRN {

    publishDir "${params.outdir}/stage1/pyscenic_output", mode: 'copy'

    input:
    path expr_matrix
    path tfs

    output:
    path "adjacencies.tsv", emit: adjacencies

    script:
    """
    arboreto_with_multiprocessing.py \
        ${expr_matrix} \
        ${tfs} \
        --method grnboost2 \
        --output adjacencies.tsv \
        --num_workers ${params.pyscenic_workers} \
        --seed 777
    """
}


process PYSCENIC_CTX {

    publishDir "${params.outdir}/stage1/pyscenic_output", mode: 'copy'

    input:
    path adjacencies
    path expr_matrix
    path ranking_10kb
    path ranking_500bp
    path annotations

    output:
    path "regulons.csv", emit: regulons

    script:
    """
    pyscenic ctx \
        ${adjacencies} \
        ${ranking_10kb} \
        ${ranking_500bp} \
        --annotations_fname ${annotations} \
        --expression_mtx_fname ${expr_matrix} \
        --output regulons.csv \
        --mask_dropouts \
        --num_workers ${params.pyscenic_workers}
    """
}


process PYSCENIC_AUCELL {

    publishDir "${params.outdir}/stage1/pyscenic_output", mode: 'copy'

    input:
    path expr_matrix
    path regulons

    output:
    path "auc_mtx.csv", emit: auc_matrix

    script:
    """
    pyscenic aucell \
        ${expr_matrix} \
        ${regulons} \
        --output auc_mtx.csv \
        --num_workers ${params.pyscenic_workers}
    """
}


process PLOT_PYSCENIC_E2F2 {

    publishDir "${params.outdir}/stage1/figures", mode: 'copy'

    input:
    path expr_matrix
    path tfs
    path ne_clean_h5ad
    path auc_matrix

    output:
    path "umap_pyscenic_e2f2.png", emit: umap_plot

    script:
    """
    python ${projectDir}/scripts/stage1/pyscenic_plots.py \
        --input_matrix ${expr_matrix} \
        --pyscenic_resources ${tfs} \
        --input_clean ${ne_clean_h5ad} \
        --auc ${auc_matrix} \
        --output_file umap_pyscenic_e2f2.png
    """
}

workflow {

    run_stage1 = params.run_stage1.toString() == 'true'
    run_stage2 = params.run_stage2.toString() == 'true'

    if (run_stage1) {

        raw_nb = LOAD_NB(
            file(params.stage1_input_dir)
        )

        processed_nb = PREPROCESS_NB(
            raw_nb.raw_h5ad
        )

        ne_cells = SUBSET_NE(
            processed_nb.processed_h5ad
        )

        ne_clean = CLEAN_NE(
            ne_cells.ne_h5ad
        )

        pyscenic_input = EXPORT_FOR_PYSCENIC(
            ne_clean.ne_clean_h5ad
        )

        grn = PYSCENIC_GRN(
            pyscenic_input.expr_matrix,
            file(params.pyscenic_tfs)
        )

        ctx = PYSCENIC_CTX(
            grn.adjacencies,
            pyscenic_input.expr_matrix,
            file(params.pyscenic_ranking_10kb),
            file(params.pyscenic_ranking_500bp),
            file(params.pyscenic_annotations)
        )

        aucell = PYSCENIC_AUCELL(
            pyscenic_input.expr_matrix,
            ctx.regulons
        )

        PLOT_PYSCENIC_E2F2(
            pyscenic_input.expr_matrix,
            file(params.pyscenic_tfs),
            ne_clean.ne_clean_h5ad,
            aucell.auc_matrix
        )
    }


    if (run_stage2) {

        geneset = file(params.e2f2_geneset)

        target_matrix = BUILD_TARGET_MATRIX(
            file(params.target_metadata),
            params.target_data_dir
        )

        target_scores = CALCULATE_TARGET_E2F2_SSGSEA(
            target_matrix.logcpm,
            geneset
        )

        target_survival = PREPARE_TARGET_SURVIVAL(
            file(params.target_clinical),
            target_scores.scores,
            file(params.target_mycn)
        )

        TARGET_SURVIVAL_ANALYSIS(
            target_survival.survival_table
        )

        gse49711 = VALIDATE_E2F2_GSE49711(
            geneset
        )

        PLOT_GSE49711(
            gse49711.merged
        )
    }
}