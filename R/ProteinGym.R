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
    # am_table <- ProteinGymR::AlphaMissense_scores()

    eh <- ExperimentHub::ExperimentHub()
    am_table <- eh[['EH9554']]
    return(am_table)
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
    # dms_table <- ProteinGymR::dms_substitutions()
        
    eh <- ExperimentHub::ExperimentHub()
    dms_table <- eh[['EH9555']]
    return(dms_table)
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
#' @importFrom queryup query_uniprot
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
        
        ## Load default AlphaMissense data from EH
        am_table <- pg_am_data()
        
        ## Rename columns to match default dms_table
        new_cols <- c('UniProt_id', 'mutant')

        am_table <- 
            am_table |> 
            rename_with(
                ~ new_cols, 
                .cols = c('Uniprot_ID', 'variant_id')
            )
        
        ## Default am_table IDs are in SwissProt. Convert to UniProt
        query <- list("accession_id" = uID)
        res <- query_uniprot(query = query, show_progress = TRUE)
        swissID <- res |> pull(`Entry Name`)
       
        ## Replace swissID observations with uID
        am_table <-
            am_table |>
            mutate(
                UniProt_id = case_when(
                    (.data$UniProt_id) == swissID ~ uID,
                    TRUE ~ as.character(UniProt_id)
                )
            )
        
        am_table
    }
        
    ## Filter for uID
    alphamissense_table <-
        am_table |>
        filter(.data$UniProt_id == uID) |> 
        as_tibble()

    ## Check if table is empty after filtering
    if (!NROW(alphamissense_table)) {
        stop(
            "no AlphaMissense information found for the protein ",
            "accession '", uID, "'; check that the UniProt ID is correct"
        )
    }
    
    alphamissense_table
}

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
    if (!NROW(dms_table)) {
        stop(
            "no DMS substitution information found for the protein ",
            "accession '", uID, "'"
        )
    }
    
    dms_table
}

#' Merge alphamissense and dms tables by UniProt and mutant IDs
#'
#' @noRd
#'
#' @importFrom dplyr left_join
#' 
pg_match_id <- 
    function(am_table, pg_table)
{
    ## Check that UniProt IDs are the same across tables
    stopifnot(
        unique(am_table$UniProt_id) == unique(pg_table$UniProt_id)
    )
        
    merged_table <- 
        left_join(
            am_table, pg_table, 
            by = c("UniProt_id", "mutant"),
            relationship = "many-to-many"
        ) |> 
        select(UniProt_id, mutant, AlphaMissense, DMS_score) |> 
        na.omit()
    
    ## Average am and dms scores across multiple studies
    merged_table <-
        merged_table |> 
        group_by(UniProt_id, mutant) |>    
        summarise(
            mean_am = mean(AlphaMissense),
            mean_dms= mean(DMS_score)
        )
    
    merged_table
}

#' Average Spearman correlation per protein
#'
#' @noRd
#'
#' @importFrom dplyr left_join mutate case_when mutate_at group_by
#'     ungroup arrange
#'
pg_correlate <- function(merged_table){
    
    cor_results <- 
        cor.test(
            merged_table$mean_am, merged_table$mean_dms, 
            method=c("spearman"), 
            exact = FALSE
        )
    
    cor_results
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
#' 
#' @references Cheng et al.,
#' Accurate proteome-wide missense variant effect prediction with AlphaMissense.
#' \emph{Science} 381, eadg7492. DOI:10.1126/science.adg7492.
#' 
#' @references Notin, P., Kollasch, A., Ritter, D., van Niekerk, L., Paul, S., 
#' Spinner, H., Rollins, N., Shaw, A., Orenbuch, R., Weitzman, R., Frazer, J., 
#' Dias, M., Franceschi, D., Gal, Y., & Marks, D. (2023). 
#' ProteinGym: Large-Scale 
#' Benchmarks for Protein Fitness Prediction and Design. In A. Oh, T. Neumann, 
#' A. Globerson, K. Saenko, M. Hardt, & S. Levine (Eds.), \emph{Advances in 
#' Neural Information Processing Systems} (Vol. 36, pp. 64331-64379). 
#' Curran Associates, Inc.
#' 
#' @importFrom ggplot2 ggplot geom_bin2d aes scale_colour_manual element_text
#'     scale_fill_manual scale_shape_manual scale_size_manual scale_fill_continuous
#'     element_blank scale_discrete_manual geom_hline labs xlab ylab 
#'     theme_classic annotate theme
#' 
#' @export
ProteinGym_correlation_plot <-
    function(uniprotId, alphamissense_table, dms_table)
{
    ## Validate uniprotId to start filtering alphamissense and clinvar tables
    stopifnot(isCharacter(uniprotId))
    
    ## Filter AM and PG tables with uniProtID
    alphamissense_table <-
        pg_filter_am_table(
            am_table = alphamissense_table,
            uID = uniprotId
        )
    
    dms_table <-
        pg_filter_dms_table(
            pg_table = dms_table,
            uID = uniprotId
        )
    
    ## Join tables by UniProt ID
    merged_table <-
        pg_match_id(am_table = alphamissense_table, pg_table = dms_table)
    
    ## Check if merged table is empty
    if (!NROW(merged_table)) {
        stop(
            "no common mutants between AlphaMissense and DMS scores for ",
            "accession '", uID, "'"
        )
    }
    
    cor_results <- pg_correlate(merged_table)
    
    ## Correlation density plot
    pg_density_plot <- 
        merged_table |> 
        ggplot(
            aes(y = .data$mean_am, x = .data$mean_dms)
        ) +
        geom_bin2d(bins = 60) +
        scale_fill_continuous(type = "viridis") +
        theme_classic() +
        labs(title = paste0("UniProt ID: ", uniprotId)) +
        xlab("DMS score") +
        ylab("AlphaMissense score") +
        theme_classic() +
        theme(
            axis.text.x = element_text(size = 16),
            axis.text.y = element_text(size = 16),
            axis.title.y = element_text(size = 16, vjust = 2),
            axis.title.x = element_text(size = 16, vjust = 0),
            legend.title = element_text(size = 16),
            legend.text = element_text(size = 16)
        ) +
        annotate("text", x = Inf, y = Inf, hjust = 1, vjust = 2,
            label = paste0("r = ", format(round(cor_results$estimate, 2)), 
                '\n',"Pval = ", cor_results$p.value),
            fontface="italic", size = 4
        )
    
    pg_density_plot
}