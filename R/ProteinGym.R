#' @rdname ProteinGym
#'
#' @title Integrate ProteinGym DMS Scores with AlphaMissense Pathogenicity 
#' Scores
#'
#' @description `ProteinGym_data()` loads ProteinGym information from the
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
#' `ProteinGym_data()` returns a tbl with 82872 rows and 5 variables:
#'
#' - `variant_id`: ClinVar variant identifier.
#' - `transcript_id`: Ensembl transcript identifier.
#' - `protein_variant`: UniProt accession:protein variant identifier.
#' - `AlphaMissense`: AlphaMissense pathogenicity score.
#' - `label`: Binary ClinVar class, either "benign" or "pathogenic".
#'
#' @examples
#' ProteinGym_data()
#' @export
#' 
ProteinGym_data <- ## load in AlphaMissense PG supplement from EH
    function (metadata = FALSE)
{
    eh <- ExperimentHub::ExperimentHub()
    title <- "ProteinGym_Supplemental"
        
    eh <- AnnotationHub::query(eh, title)
    ehid <- eh$ah_id
        
    if (metadata == TRUE) {
        eh[ehid]
    }
    else eh[[ehid]]
}

#' @rdname ProteinGym
#'
#' @description `ProteinGym_correlate()` runs a Spearman correlation between 
#'    ProteinGym deep mutational scanning (DMS) assay scores with AlphaMissense 
#'    predicted scores. Returns a ggplot object for visualization.
#'
#' @param uniprotId `character()` a valid UniProt accession identifier.
#'
#' @param alphamissense_table a table containing AlphaMissense
#'    predictions for protein variants matching ProteinGym variants. The default 
#'    table is derived from the supplemental data of the AlphaMissense paper. 
#'    Alternatively, a user-defined [`tibble::tbl_df`] or [`data.frame`]
#'    can be supplied.
#'    
#' @param DMS_table a table containing ProteinGym deep mutational scanning (DMS) 
#'    assay scores for protein variants matching AlphaMissense variants. The
#'    default table loads in subsitutions from 
#'    [ProteinGym](https://proteingym.org/download).
#'    Alternatively, a user-defined [`tibble::tbl_df`] or [`data.frame`]
#'    can be supplied.
#'
#' @details
#'
#' For `ProteinGym_correlate()`, `alphamissense_table` columns must include:
#'
#' - `uniprot_id`: UniProt accession identifier.
#' - `DMS`
#' - `protein_variant`: variant identifier string, with protein
#'    position in the middle, and the reference and mutant amino acid
#'    residues to the left and right of the position, respectively.
#' - `DMS`: AlphaMissense classification of either "benign",
#'    "ambiguous", or "pathogenic".
#' - `am_pathogenicity`: AlphaMissense predicted score.
#'
#' `DMS_table` columns must include:
#'
#' - `uniprot_id`: UniProt accession identifier, matching `alphamissense_table`.
#' - `protein_variant`: variant identifier string, matching
#'    `alphamissense_table` format.
#' - `cv_class`: binary ClinVar classification of "benign" or "pathogenic".
#'
#' @return `clinvar_plot()` returns a `ggplot` object which overlays
#'    ClinVar classifications onto AlphaMissense predicted
#'    scores. Blue, gray, and red colors represent pathogenicity
#'    classifications for "likely benign", "ambiguous", or
#'    "likely pathogenic", respectively. Large, bolded points are
#'    ClinVar variants colored according to their clinical
#'    classification, while smaller points in the background are
#'    AlphaMissense predictions.
#'
#' @examples
#'
#' alphamissense_table <- am_data("aa_substitutions")
#'
#' clinvar_plot(
#'     uniprotId = "P37023",
#'     alphamissense_table = alphamissense_table
#' )
#'
#' @references Cheng et al.,
#' Accurate proteome-wide missense variant effect prediction with AlphaMissense.
#' \emph{Science} 381, eadg7492. DOI:10.1126/science.adg7492.
#'
#' @importFrom BiocBaseUtils isCharacter
#'
#' @importFrom utils data
#'
#' @importFrom dplyr is.tbl
#'
#' @export
#' 
ProteinGym_correlation <- function(uniprotId, alphamissense_table, clinvar_table){
    
    results <- cor.test(df$am_pathogenicity, df$DMS_score, method=c("spearman"), 
                        exact = FALSE)
    
    results <- abs(results$estimate)
    
    results <- unname(results)
    
    return(results)
}

uniprotId = "VKOR1"
### A function to extract ProteinGym studies matching UniProt specified in function

# Extract all ProteinGym studies matching UniProtId
# Extract only the UniProtId from the study names (first and second elements)
PG_names <- strsplit(names(RDS_ProGym), split = "_", fixed = TRUE)

# Define a function to extract the first two elements, apply to list elements
# Take the first element - SwissProtID
ProtId <- sapply(PG_names, function(input_string) {
    first_two_elements <- paste(input_string[1], input_string[2], sep = "_")
    return(first_two_elements)
})

# convert it to 

cPG <- all_ProGym[names(all_ProGym) %in% uniprotId]

indx <- paste("HUMAN")
ProGym_names <- names(all_ProGym)
human_studies <- ProGym_names[grep(indx, ProGym_names)] 
# 32 human studies

# Subset to ProteinGym list to only human studies
human_ProGym <- all_ProGym[names(all_ProGym) %in% human_studies]
human_ProGym