#' @rdname ProteinGym
#'
#' Load in AlphaMissense DMS supplemental table via ExperimentHub 
#' Through ProteinGymR later
#' 
#' @noRd
#'
#' @importFrom dplyr filter as_tibble vars
#'
pg_am_data <-
    function()
{
    # am_data <- ProteinGymR::AlphaMissense_scores()

    eh <- ExperimentHub::ExperimentHub()
    am_data <- eh[['EH9554']]
    return(am_data)
}

#' Load in ProteinGym DMS dataset via ExperimentHub (through ProteinGymR later)
#' 
#' @noRd
#'
#' @importFrom dplyr filter as_tibble vars
#' 
pg_dms_data <-
    function()
{
    # dms_data <- ProteinGymR::dms_substitutions()
        
    eh <- ExperimentHub::ExperimentHub()
    dms_data <- eh[['EH9555']]
    return(dms_data)
}

#' Map Swiss-Prot entry name to UniProt Accession ID
#' Later, use ProteinGymR function
#'
#' @noRd
#'
#' @import UniProt.ws
#' 
pg_map_accessions <-
    function(entryNames)
{
    # Convert SwissProt entries to UniProt accession ID
    ws <- UniProt.ws::UniProt.ws()
    out <- UniProt.ws::select(ws, entryNames, columns = "UniProtKB", 
                              keytype = "UniProtKB")
    accessions <- out$Entry
    
    return(accessions)
}

#' Filter the AlphaMissense table with UniprotID
#'
#' @noRd
#' 
#' @importFrom dplyr filter select pull as_tibble vars rename mutate
#' 
pg_filter_am_table <-
    function(am_table, uID)
{

    ## Check if am_table is missing
    if (missing(am_table)) {
        spdl::info(paste(
            "'alphamissense_table' not provided, using default table from",
            "`ProteinGymR::am_scores()`"
        ))
        
        am_table <- pg_am_data()
        
        ## Add UniProt accessions to am_table
        swissprot_names <- am_table |> 
            select(.data$Uniprot_ID) |> 
            unique() |> pull()
        
        acc <- pg_map_accessions(swissprot_names)
        
        accessions_lookup <-
            cbind(swissprot_names, acc) |> 
            as.data.frame()
        
        ## Default am_table uses SwissProt
        selected_swiss_protein <- 
            accessions_lookup |> 
            filter(.data$acc == uID) |> 
            pull(.data$swissprot_names)
        
        ## Replace SwissProt with corresponding UniProt ID in am_table
        am_table <-
            am_table |>
            mutate(
                Uniprot_ID = case_when(
                    (.data$Uniprot_ID) == selected_swiss_protein ~ uID,
                    TRUE ~ as.character(Uniprot_ID)
                )
            ) |> 
            as_tibble()
        
        ## Rename columns to match default dms_table
        new_cols <- c('UniProt_ID', 'mutant')

        am_table <- 
            am_table |> 
            rename_with(
                ~ new_cols, 
                .cols = c('Uniprot_ID', 'variant_id')
            )
        
        am_table
    }
        
    ## Filter for uID
    alphamissense_table <-
        am_table |>
        filter(.data$UniProt_ID == uID) |>
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
#' Filter the ProteinGym DMS table with UniprotID
#'
#' @noRd
#'
#' @importFrom dplyr filter as_tibble vars bind_rows
#' @importFrom purrr keep
#'
pg_filter_dms_table <-
    function(pg_table, uID)
{
    ## Check if pg_table is missing
    if (missing(pg_table)) {
        spdl::info(paste(
            "'dms_table' not provided, using default table from",
            "`ProteinGymR::dms_substitutions()`"
        ))
        
        pg_table <- pg_dms_data()
    }
    
    ## Take pg_table and filter for the uID, rbind into one data.frame
    filtered_dfs <- purrr::keep(pg_table, ~ any(.x$UniProt_id == uID))
    
    dms_table <-
        filtered_dfs |>
        bind_rows() |>
        as_tibble()
    
    ## Check if table is empty after filtering
    ## This will work for a tibble or a data.frame
    if (!NROW(dms_table)) {
        stop(
            "no DMS substitution information found for the protein ",
            "accession '", uID, "'"
        )
    }
    
    dms_table
}

#' Prepare dms and alphamissense data for correlation plotting
#'
#' @noRd
#'
#' @importFrom dplyr left_join mutate case_when mutate_at group_by
#'     ungroup arrange
#'
pg_prepare_data_for_plot <-
    function(am_table, pg_table)
{
    ## Merge am_table and pg_table by `DMS_id` and `mutant`
    ## Create a temporary identifier identical between datasets `tmp_id`
    pg_table <- 
        pg_table |> 
        mutate(tmp_id = paste(UniProt_id, mutant, sep = "_")) |> 
        select()
    
    am_table <- 
        am_table |> 
        mutate(tmp_id = paste(Uniprot_ID, variant_id, sep = "_"))
    
    
    
    both_table <- 
        left_join(
            am_table, pg_table, 
            by = "tmp_id",
            relationship = "many-to-many"
        )
        
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
#' @title Integrate ProteinGym DMS and AlphaMissense Pathogenicity Scores
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
