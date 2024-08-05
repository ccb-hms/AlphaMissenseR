#' @rdname ProteinGym
#'
#' @title Integrate ProteinGym DMS and AlphaMissense Pathogenicity Scores
#'
#' @description `ProteinGym_AlphaMissense_data()` loads in the AlphaMissense 
#'    pathogenicity scores for mutants in ProteinGym from the AlphaMissense 
#'    publication by Cheng et al. 
#'    ([2023](https://www.science.org/doi/10.1126/science.adg7492)).
#'
#' @return
#'
#' `ProteinGym_AlphaMissense_data()` returns a data.frame with 1622429 rows and
#'    4 variables:
#'
#' - `DMS_id`: ProteinGym assay identifier.
#' - `Uniprot_ID`: UniProt identifer. Note, these are not accession codes.
#' - `variant_id`: Mutant identifier string matching ProteinGym. 
#'    Protein position in the middle, and the reference and mutant 
#'    amino acid residues to the left and right of the position, respectively.
#' - `AlphaMissense`: AlphaMissense pathogenicity score.
#'
#' @examples
#' am_table <- ProteinGym_AlphaMissense_data()
#' 
#' @export
ProteinGym_AlphaMissense_data <-
    function ()
{
    eh <- ExperimentHub::ExperimentHub()
    data <- eh[['EH9554']]
    return(data)
}
#'
#'
#' @rdname ProteinGym
#'
#' @description `ProteinGym_DMS_data()` loads in 216 ProteinGym deep mutational 
#' scanning assays (DMS) for substitutions. The data is provided by Notin et. al
#' [(2023)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10723403/).
#'
#' @return
#'
#' `ProteinGym_DMS_data()` returns a list of 216 data.frames corresponding to
#'    individual DMS assays. Each table contains the following 6 columns:
#'
#' - `UniProt_id`: UniProt accession identifier.
#' - `DMS_id`: ProteinGym assay identifier.
#' - `mutant`: Mutant identifier string matching AlphaMissense variants. 
#'    Specifically, the set of substitutions to apply on the reference sequence 
#'    to obtain the mutated sequence (e.g., A1P:D2N implies the amino acid 'A' 
#'    at position 1 should be replaced by 'P', and 'D' at position 2 should be 
#'    replaced by 'N').
#' - `mutated_sequence`: Full amino acid sequence for the mutated protein.
#' - `DMS_score`: Experimental measurement in the DMS assay. 
#'    Higher values indicate higher fitness of the mutated protein.
#' - `DMS_score_bin`: Factor, indicating whether the DMS_score is 
#'    above the fitness cutoff (1 is fit, 0 is not fit).
#'
#' @examples
#' pg_data <- ProteinGym_DMS_data()
#' 
#' @export
ProteinGym_DMS_data <-
    function ()
{
    eh <- ExperimentHub::ExperimentHub()
    data <- eh[['EH9555']]
    return(data)
}
#'
#' Map Swiss-Prot entry name to UniProt Accession ID
#'
#' @noRd
#'
#' @import UniProt.ws
#' 
map_accessions <-
    function(entryNames)
{
    # Convert entryNames to UniProt Accession ID
    ws <- UniProt.ws::UniProt.ws()
    out <- UniProt.ws::select(ws, entryNames, columns = "UniProtKB", 
                              keytype = "UniProtKB")
    accessions <- out$Entry
    
    return(accessions)
}
#' 
#'
#' Filter the AlphaMissense table with UniprotID
#'
#' @noRd
#'
#' @importFrom dplyr filter as_tibble vars
#'
pg_filter_am_table <-
    function(am_table, uID)
{
    browser()
    ## Check if am_table is missing
    if (missing(am_table)) {
        spdl::info(paste(
            "'alphamissense_table' not provided, using default table from",
            "`ProteinGym_AlphaMissense_data()`"
        ))
        
        am_table <- ProteinGym_AlphaMissense_data()
        
        ## Add UniProt accessions to am_table
        swissprot_names <- am_table |> 
            select(.data$Uniprot_ID) |> 
            unique() |> pull()
        
        acc <- map_accessions(swissprot_names)
        
        accessions_lookup <-
            cbind(swissprot_names, acc) |> 
            as.data.frame()
        
        ## Default able uses SwissProt. Replace SwissProt with uID
        selected_swiss_protein <- 
            accessions_lookup |> 
            filter(.data$acc == uID) |> 
            pull(.data$swissprot_names)
        
        am_table <-
            am_table |>
            mutate(
                Uniprot_ID = case_when(
                    (.data$Uniprot_ID) == selected_swiss_protein ~ uID,
                    TRUE ~ as.character(Uniprot_ID)
                )
            )
        am_table
    }
        
    ## Take alphamissense_table and filter for the uniprotId
    alphamissense_table <-
        am_table |>
        filter(.data$Uniprot_ID == uID) |>
        as_tibble()

    ## Check if table is empty after filtering
    ## This will work for a tibble or a data.frame
    if (!NROW(alphamissense_table)) {
        stop(
            "no AlphaMissense information found for the protein ",
            "accession '", uID, "'; check that the UniProt ID is correct"
        )
    }
    
    alphamissense_table
}
#'
#'
#'#' Filter the DMS table with UniprotID
#'
#' @noRd
#'
#' @importFrom dplyr filter as_tibble vars
#'
pg_filter_DMS_table <-
    function(pg_table, uID)
    {
        ## Check if am_table is missing
        if (missing(pg_table)) {
            spdl::info(paste(
                "'DMS_table' not provided, using default table from",
                "`ProteinGym_DMS_data()`"
            ))
            
            pg_table <- ProteinGym_DMS_data()
        }
        
        ## Take alphamissense_table and filter for the uniprotId
        DMS_table <-
            pg_table |>
            filter(.data$UniProt_id == uID) |>
            as_tibble()
        
        ## Check if table is empty after filtering
        ## This will work for a tibble or a data.frame
        if (!NROW(DMS_table)) {
            stop(
                "no DMS information found for the protein ",
                "accession '", uID, "'"
            )
        }
        
        DMS_table
    }
#'
#'
#' Prepare data for the function ProteinGym_correlation_plot
#'
#' @noRd
#'
#' @importFrom dplyr left_join mutate case_when mutate_at group_by
#'     ungroup arrange
#'
proteingym_prepare_data_for_plot <-
    function(am_table, pg_table)
{
    ## grab amino acid positions
    am_table <- mutate(
        am_table,
        aa_pos = as.integer(
            gsub(".*?([0-9]+).*", "\\1", .data$protein_variant)
        )
    )
    
    ## join datasets
    combined_data <- left_join(
        am_table,
        cv_table,
        by = c('uniprot_id', 'protein_variant')
    )
    
    ## add color code matching AM and CV labels
    combined_data <-
        combined_data |>
        mutate(
            code_color = case_when(
                !is.na(.data$cv_class) & .data$cv_class == "benign" ~
                    "CV benign",
                !is.na(.data$cv_class) & .data$cv_class == "pathogenic" ~
                    "CV pathogenic",
                is.na(.data$cv_class) & .data$am_class == "pathogenic" ~
                    "AM pathogenic",
                is.na(.data$cv_class) & .data$am_class == "benign" ~
                    "AM benign",
                is.na(.data$cv_class) & .data$am_class == "ambiguous" ~
                    "AM ambiguous")
        ) |>
        mutate_at(vars(.data$code_color), factor) |>
        arrange(.data$code_color)
    
    ## Grab the thresholds for AM pathogenicity to plot
    combined_data <-
        combined_data |>
        group_by(.data$am_class) |>
        mutate(
            max = max(.data$am_pathogenicity, na.rm=TRUE),
            min = min(.data$am_pathogenicity, na.rm=TRUE)
        ) |>
        ungroup()
    
    combined_data
}
#'
#' @rdname ProteinGym
#' 
#' @description `ProteinGym_correlation_plot()` runs a Spearman correlation between 
#'    ProteinGym deep mutational scanning (DMS) assay scores against 
#'    AlphaMissense predicted scores. Returns a ggplot object for visualization.
#'
#' @param uniprotId `character()` a valid UniProt accession identifier.
#' 
#' @param alphamissense_table a table containing AlphaMissense predictions 
#'    for variants matching ProteinGym substitution mutants. The default 
#'    table is derived from the supplemental data of the AlphaMissense paper. 
#'    Alternatively, a user-defined [`tibble::tbl_df`] or [`data.frame`]
#'    can be supplied.
#'
#' @param DMS_table a table containing deep mutational scanning (DMS) 
#'    assay scores for protein variants matching AlphaMissense. The default 
#'    table loads in substitutions from 
#'    [ProteinGym](https://proteingym.org/download).
#'    Alternatively, a user-defined [`tibble::tbl_df`] or [`data.frame`]
#'    can be supplied.
#'
#' @details
#'
#' For `ProteinGym_correlation_plot()`, 
#'    `alphamissense_table` columns must include:
#'
#' - `Uniprot_ID`: UniProt accession identifier.
#' - `variant_id`: Mutant identifier string matching ProteinGym. 
#'    Protein position in the middle, and the reference and mutant 
#'    amino acid residues to the left and right of the position, respectively.
#' - `AlphaMissense`: AlphaMissense pathogenicity score.
#'
#' `DMS_table` columns must include:
#'
#' - `UniProt_id`: UniProt accession identifier.
#' - `mutant`: Mutant identifier string matching AlphaMissense variants. 
#'    Specifically, the set of substitutions to apply on the reference sequence 
#'    to obtain the mutated sequence (e.g., A1P:D2N implies the amino acid 'A' 
#'    at position 1 should be replaced by 'P', and 'D' at position 2 should be 
#'    replaced by 'N').
#' - `DMS_score`: Experimental measurement in the DMS assay. 
#'    Higher values indicate higher fitness of the mutated protein.
#'
#' @return `ProteinGym_correlation_plot()` returns a `ggplot` object visualizing the 
#'    Spearman correlation between experimental DMS scores and AlphaMissense 
#'    predicted scores. In our case, a stronger negative correlation correspond 
#'    to a stronger relationship between the two measures.
#'
#' @examples
#' 
#' ProteinGym_correlation_plot(uniprotId = "Q9NV35")
#' 
#' @export
ProteinGym_correlation_plot <-
    function(uniprotId, alphamissense_table, DMS_table)
{
    ## Validate uniprotId to start filtering alphamissense and clinvar tables
    stopifnot(isCharacter(uniprotId))
    
    ## Filter AM and PG tables with uniProtID
    alphamissense_table <-
        pg_filter_am_table(
            am_table = alphamissense_table,
            uID = uniprotId
        )
    
    DMS_table <-
        pg_filter_DMS_table(
            pg_table = DMS_table,
            uID = uniprotId
        )
    
 # results <- cor.test(df$am_pathogenicity, df$DMS_score, method=c("spearman")) 
                            #                         exact = FALSE)
                            #     
                            #     results <- abs(results$estimate)
                            #     
                            #     results <- unname(results)
                            #     
                            #     return(results)
}

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
# ProteinGym_correlation <- function(uniprotId, alphamissense_table, clinvar_table){
#     
#     results <- cor.test(df$am_pathogenicity, df$DMS_score, method=c("spearman"), 
#                         exact = FALSE)
#     
#     results <- abs(results$estimate)
#     
#     results <- unname(results)
#     
#     return(results)
# }
# 
# uniprotId = "VKOR1"
# ### A function to extract ProteinGym studies matching UniProt specified in function
# 
# # Extract all ProteinGym studies matching UniProtId
# # Extract only the UniProtId from the study names (first and second elements)
# PG_names <- strsplit(names(RDS_ProGym), split = "_", fixed = TRUE)
# 
# # Define a function to extract the first two elements, apply to list elements
# # Take the first element - SwissProtID
# ProtId <- sapply(PG_names, function(input_string) {
#     first_two_elements <- paste(input_string[1], input_string[2], sep = "_")
#     return(first_two_elements)
# })
# 
# # convert it to 
# 
# cPG <- all_ProGym[names(all_ProGym) %in% uniprotId]
# 
# indx <- paste("HUMAN")
# ProGym_names <- names(all_ProGym)
# human_studies <- ProGym_names[grep(indx, ProGym_names)] 
# # 32 human studies
# 
# # Subset to ProteinGym list to only human studies
# human_ProGym <- all_ProGym[names(all_ProGym) %in% human_studies]
# human_ProGym
# 
# 
# 
# ### ADDITIONAL CODE:
# all_ProGym <- readr::read_csv("~/AlphaMissense_data/Supplementary_Data_S8_proteingym.csv.gz")
# head(all_ProGym)
# 
# 
# RDS_ProGym <- readRDS("~/AlphaMissense_data/ProteinGym/all_ProGym.RDS")
# 
# # Spearman Correlation
# 
# #The absolute value of the Spearman correlation between predicted and observed
# #assay scores were calculated per MAVE experiment then averaged by UniProt ID.
# 
# #- Calculate Spearman Correlation between AM and MAVE assay scores for each MAVE
# #experiment. Make this an absolute value.
# 
# #- Then, average this value by UniprotID (per protein across all MAVE studies).
# 
# 
# # Spearman correlation function
# 
# ProteinGym_correlation <- function(uniprotId, alphamissense_table, clinvar_table){
#     
#     results <- cor.test(df$am_pathogenicity, df$DMS_score, method=c("spearman"), 
#                         exact = FALSE)
#     
#     results <- abs(results$estimate)
#     
#     results <- unname(results)
#     
#     return(results)
# }
# 
# 
# # Apply the function to each pair of data frames in the lists using lapply
# correlation_results <- lapply(1:length(filtered_dataframes), function(i) {
#     
#     cor_AM_ProGym(filtered_dataframes[[i]])
#     
# })
# 
# 
# #Average the absolute Spearman estimate by protein using the UniprotID.
# 
# # We can grab this information from the "UniprotID" column of each dataset.
# unique_vector <- unlist(lapply(filtered_dataframes, function(df) unique(df$uniprot_id)))
# 
# # Assign UniProt names to the correlation_results
# names(correlation_results) <- unique_vector
# 
# # Grab the unique vector (remove dup protein names)
# UniProt_IDs <- unique(unname(unique_vector))
# 
# # Average Spearman across UniProt
# avg_spearman <- lapply(1:length(UniProt_IDs), function(i) {
#     
#     # Grab each unique protein ID
#     c.id <- UniProt_IDs[i]
#     
#     # Grab the index/indices of that protein from the Spearman results list
#     idx <- which(names(correlation_results) == c.id)
#     
#     # Average the correlation
#     results <- mean(unlist(correlation_results[idx]))
# 
#     return(results)
#     
# })
# 
# # Name the results with the protein ID
# names(avg_spearman) <- UniProt_IDs