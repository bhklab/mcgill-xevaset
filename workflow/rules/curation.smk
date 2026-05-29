from pathlib import Path


rawdata = Path(config["directories"]["rawdata"])
procdata = Path(config["directories"]["procdata"])
results = Path(config["directories"]["results"])
logs = Path(config["directories"]["logs"])
scripts = Path("../scripts")

raw_inputs = config["raw_inputs"]


rule parse_tumor_workbook:
    input:
        workbook=rawdata / raw_inputs["tumor_workbook"],
    output:
        measurements=procdata / "treatmentResponse" / "tumor_measurements.csv",
        block_summary=procdata / "metadata" / "tumor_workbook_blocks.csv",
    log:
        logs / "parse_tumor_workbook.log",
    script:
        scripts / "parse_tumor_workbook.R"


rule build_metadata:
    input:
        measurements=rules.parse_tumor_workbook.output.measurements,
        clinical=rawdata / raw_inputs["clinical_info"],
        crosswalk=rawdata / raw_inputs["sample_crosswalk"],
    output:
        model=procdata / "metadata" / "model_raw.csv",
        sample=procdata / "metadata" / "sample_raw.csv",
        treatment=procdata / "metadata" / "treatment_raw.csv",
    log:
        logs / "build_metadata.log",
    script:
        scripts / "build_metadata.R"


rule annotate_metadata:
    input:
        model=rules.build_metadata.output.model,
        sample=rules.build_metadata.output.sample,
        treatment=rules.build_metadata.output.treatment,
    output:
        model=procdata / "metadata" / "annotations" / "McGill_modelMetadata_annotated.csv",
        sample=procdata / "metadata" / "annotations" / "McGill_sampleMetadata_annotated.csv",
        treatment=procdata / "metadata" / "annotations" / "McGill_treatmentMetadata_annotated.csv",
        unmatched_treatments=results / "unmatched_treatments.csv",
    log:
        logs / "annotate_metadata.log",
    script:
        scripts / "annotate_metadata.R"


rule make_RNASeq_SE:
    input:
        rnaseq=rawdata / raw_inputs["rnaseq"],
        clinical=rawdata / raw_inputs["clinical_info"],
        crosswalk=rawdata / raw_inputs["sample_crosswalk"],
    output:
        se=procdata / "rnaseq" / "RNAseq_SE.rds",
        matrix=procdata / "rnaseq" / "RNAseq_expression.csv",
    log:
        logs / "make_RNASeq_SE.log",
    script:
        scripts / "make_RNASeq_SE.R"


rule make_Mutation_SE:
    input:
        mutation=rawdata / raw_inputs["mutation"],
        crosswalk=rawdata / raw_inputs["sample_crosswalk"],
    output:
        se=procdata / "mutation" / "Mutation_SE.rds",
        matrix=procdata / "mutation" / "Mutation_expression.csv",
    log:
        logs / "make_Mutation_SE.log",
    script:
        scripts / "make_Mutation_SE.R"


rule make_CNV_SE:
    input:
        cnv=rawdata / raw_inputs["cnv"],
        crosswalk=rawdata / raw_inputs["sample_crosswalk"],
    output:
        se=procdata / "cnv" / "CNV_SE.rds",
        matrix=procdata / "cnv" / "CNV_expression.csv",
    log:
        logs / "make_CNV_SE.log",
    params:
        datatype="CNV",
    script:
        scripts / "make_CNV_SE.R"


rule make_Fusion_SE:
    input:
        fusion=rawdata / raw_inputs["fusion"],
        crosswalk=rawdata / raw_inputs["sample_crosswalk"],
    output:
        se=procdata / "fusion" / "Fusion_SE.rds",
        matrix=procdata / "fusion" / "Fusion_presence.csv",
    log:
        logs / "make_Fusion_SE.log",
    script:
        scripts / "make_Fusion_SE.R"


rule make_Segment_SE:
    input:
        segments=rawdata / raw_inputs["segments"],
        crosswalk=rawdata / raw_inputs["sample_crosswalk"],
    output:
        se=procdata / "segments" / "Segment_SE.rds",
        matrix=procdata / "segments" / "Segment_expression.csv",
    log:
        logs / "make_Segment_SE.log",
    script:
        scripts / "make_Segment_SE.R"


rule annotate_gene_metadata:
    input:
        rnaseq=rules.make_RNASeq_SE.output.se,
        mutation=rules.make_Mutation_SE.output.se,
        cnv=rules.make_CNV_SE.output.se,
        fusion=rules.make_Fusion_SE.output.se,
        segments=rules.make_Segment_SE.output.se,
    output:
        rnaseq=procdata / "gene_annotation" / "RNAseq_SE.rds",
        mutation=procdata / "gene_annotation" / "Mutation_SE.rds",
        cnv=procdata / "gene_annotation" / "CNV_SE.rds",
        fusion=procdata / "gene_annotation" / "Fusion_SE.rds",
        segments=procdata / "gene_annotation" / "Segment_SE.rds",
        unmapped=results / "unmapped_genes.csv",
    log:
        logs / "annotate_gene_metadata.log",
    script:
        scripts / "annotate_gene_metadata.R"


rule build_modToBiobaseMap:
    input:
        model=rules.annotate_metadata.output.model,
        rnaseq=rules.annotate_gene_metadata.output.rnaseq,
        mutation=rules.annotate_gene_metadata.output.mutation,
        cnv=rules.annotate_gene_metadata.output.cnv,
        fusion=rules.annotate_gene_metadata.output.fusion,
        segments=rules.annotate_gene_metadata.output.segments,
    output:
        modToBiobaseMap=procdata / "metadata" / "annotations" / "McGill_modToBiobaseMap.csv",
    log:
        logs / "build_modToBiobaseMap.log",
    script:
        scripts / "build_modToBiobaseMap.R"


rule build_experimentDesign:
    input:
        model=rules.annotate_metadata.output.model,
    output:
        expDesign=procdata / "metadata" / "annotations" / "McGill_expDesign.rds",
    log:
        logs / "build_experimentDesign.log",
    script:
        scripts / "build_experimentDesign.R"


rule build_MultiAssayExperiment:
    input:
        sample=rules.annotate_metadata.output.sample,
        assays=[
            rules.annotate_gene_metadata.output.rnaseq,
            rules.annotate_gene_metadata.output.mutation,
            rules.annotate_gene_metadata.output.cnv,
            rules.annotate_gene_metadata.output.fusion,
            rules.annotate_gene_metadata.output.segments,
        ],
    output:
        mae=results / "McGill_MultiAssayExperiment.rds",
    log:
        logs / "build_MultiAssayExperiment.log",
    script:
        scripts / "build_MultiAssayExperiment.R"


rule build_XevaSet:
    input:
        multiAssayExperiment=rules.build_MultiAssayExperiment.output.mae,
        model=rules.annotate_metadata.output.model,
        treatment=rules.annotate_metadata.output.treatment,
        experiment=rules.parse_tumor_workbook.output.measurements,
        expDesign=rules.build_experimentDesign.output.expDesign,
        modToBiobaseMap=rules.build_modToBiobaseMap.output.modToBiobaseMap,
    output:
        xeva=results / "Xeva_McGill.rds",
        csvDir=directory(results / "Xeva_McGill_csv"),
    log:
        logs / "build_XevaSet.log",
    script:
        scripts / "build_XevaSet.R"


rule archive_Xeva_csv:
    input:
        csvDir=rules.build_XevaSet.output.csvDir,
    output:
        archive=results / "Xeva_McGill_csv.tar.gz",
    log:
        logs / "archive_Xeva_csv.log",
    shell:
        """
        set -euo pipefail
        mkdir -p $(dirname {output.archive}) $(dirname {log})
        COPYFILE_DISABLE=1 tar --exclude="Xeva_McGill_csv/.snakemake_timestamp" -C data/results -cf - Xeva_McGill_csv | gzip -9 > {output.archive} 2> {log}
        """
