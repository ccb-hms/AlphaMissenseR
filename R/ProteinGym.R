#' @rdname ProteinGym
#'
#' @title Integrate ProteinGym Mave Scores with AlphaMissense Pathogenicity 
#' Scores
#'
#' @description `proteingym_data()` loads ProteinGym information from the
#'     supplemental table of the AlphaMissense
#'     [\[2023\]](https://www.science.org/doi/10.1126/science.adg7492)
#'     paper.
#'
#' @param record character(1) Zenodo record for the AlphaMissense data
#'     resources.
#'
#' @param bfc an object returned by `BiocFileCache()` representing the
#'     location of the AlphaMissenseR database. The default is the
#'     'global' BiocFileCache.
#'
#' @return
#'
#' `proteingym_data()` returns a tbl with 82872 rows and 5 variables:
#'
#' - `variant_id`: ClinVar variant identifier.
#' - `transcript_id`: Ensembl transcript identifier.
#' - `protein_variant`: UniProt accession:protein variant identifier.
#' - `AlphaMissense`: AlphaMissense pathogenicity score.
#' - `label`: Binary ClinVar class, either "benign" or "pathogenic".
#'
#' @examples
#' proteingym_data()
#' @export
proteingym_data <-
    function(record = ALPHAMISSENSE_RECORD, bfc = BiocFileCache())
    {
        db_rname <- paste0("AlphaMissense_", record)
        db_tbl_name <- "proteingym"
        if (!NROW(bfcquery(bfc, db_rname))) {
            spdl::info("creating AlphaMissense database for record '{}'", record)
            bfcnew(bfc, db_rname)
        }
        browser()
        fpath <- system.file(
            package = "AlphaMissenseR", "extdata", 
            "Supplementary_Data_S8_proteingym.csv.gz"
        )
        db_table(
            record, bfc, db_tbl_name, fpath = fpath,
            template = "import_csv", delim = ","
        )
    }
