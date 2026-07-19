# Formats z-score into dataframe that can be used for plotting

format_score <- function(score_vec, signature_label) {
  if (!is.character(signature_label) ||
      length(signature_label) != 1 ||
      is.na(signature_label) ||
      signature_label == "") {
    stop("signature_label must be a single string.")
  }
  
  if (is.null(names(score_vec))) {
    stop("score_vec must be a named vector.")
  }
  
  score_col <- paste0(signature_label, "_z_score")
  sample_names <- names(score_vec)
  
  score_df <- data.frame(
    cell_line = sub("_.*", "", sample_names),
    parental = ifelse(grepl("_PARENTAL_", sample_names), "Parental", "MP1"),
    score_vec,
    check.names = FALSE
  )
  
  names(score_df)[names(score_df) == "score_vec"] <- score_col
  
  return(score_df)
}

# Create 2 complementary plots of signature activity
# 1. grouped_dot_plot
#    Shows Parental and MP1 samples separately for each cell line
# 2. interaction_plot
#    Connects the group means for each cell line, making it easier
#    to see whether the Parental-to-MP1 change differs by cell line

plot_signature_activity <- function(
    score_df,
    signature_label,
    x_axis_group = "parental",
    stratify_group = "cell_line",
    plot_title = NULL,
    y_axis_label = NULL,
    out_file = NULL,
    save_plot = TRUE
) {
  
  # ---------------------------------------------------------------------------
  # Validate the signature label
  # ---------------------------------------------------------------------------
  
  if (!is.character(signature_label) ||
      length(signature_label) != 1 ||
      is.na(signature_label) ||
      signature_label == "") {
    stop("signature_label must be a single non-empty string.")
  }
  
  # Construct the expected score-column name.
  #
  # Example:
  #   signature_label = "PPARG"
  #   score_col = "PPARG_z_score"
  score_col <- paste0(signature_label, "_z_score")
  
  if (!score_col %in% colnames(score_df)) {
    stop(
      paste0(
        "Column '",
        score_col,
        "' not found in score_df."
      )
    )
  }
  
  # ---------------------------------------------------------------------------
  # Validate the grouping variables
  # ---------------------------------------------------------------------------
  
  # These are the only two columns that can control the plot grouping.
  valid_groups <- c("cell_line", "parental")
  
  if (!x_axis_group %in% valid_groups) {
    stop(
      paste0(
        "x_axis_group must be either 'cell_line' or 'parental'. ",
        "Received: '",
        x_axis_group,
        "'."
      )
    )
  }
  
  if (!stratify_group %in% valid_groups) {
    stop(
      paste0(
        "stratify_group must be either 'cell_line' or 'parental'. ",
        "Received: '",
        stratify_group,
        "'."
      )
    )
  }
  
  # The same column cannot control both plot roles.
  #
  # Valid:
  #   x_axis_group = "parental"
  #   stratify_group = "cell_line"
  #
  # Also valid:
  #   x_axis_group = "cell_line"
  #   stratify_group = "parental"
  #
  # Invalid:
  #   x_axis_group = "parental"
  #   stratify_group = "parental"
  if (x_axis_group == stratify_group) {
    stop(
      paste0(
        "x_axis_group and stratify_group must be different. ",
        "Use one of the following combinations:\n",
        "  x_axis_group = 'parental', stratify_group = 'cell_line'\n",
        "  x_axis_group = 'cell_line', stratify_group = 'parental'"
      )
    )
  }
  
  # Confirm that both selected columns exist in the data frame.
  if (!x_axis_group %in% colnames(score_df)) {
    stop(
      paste0(
        "Column '",
        x_axis_group,
        "' not found in score_df."
      )
    )
  }
  
  if (!stratify_group %in% colnames(score_df)) {
    stop(
      paste0(
        "Column '",
        stratify_group,
        "' not found in score_df."
      )
    )
  }
  
  # ---------------------------------------------------------------------------
  # Create default labels
  # ---------------------------------------------------------------------------
  
  if (is.null(plot_title)) {
    plot_title <- paste(signature_label, "Activity")
  }
  
  if (is.null(y_axis_label)) {
    y_axis_label <- paste(
      signature_label,
      "Signature Z-Score"
    )
  }
  
  # Human-readable labels for plot legends.
  group_labels <- c(
    "cell_line" = "Cell line",
    "parental" = "Condition"
  )
  
  # ---------------------------------------------------------------------------
  # Define category order
  # ---------------------------------------------------------------------------
  
  # These levels determine the order displayed on the plot.
  group_levels <- list(
    parental = c(
      "Parental",
      "MP1"
    ),
    cell_line = c(
      "BBN963",
      "UPPL1541"
    )
  )
  
  # ---------------------------------------------------------------------------
  # Define colors
  # ---------------------------------------------------------------------------
  
  # Colors for each possible grouping variable.
  #
  # The function automatically selects the correct palette depending
  # on which column is used for each plot role.
  group_colors <- list(
    
    parental = c(
      "Parental" = "#556B2F",
      "MP1" = "#D8A3B0"
    ),
    
    cell_line = c(
      "BBN963" = "#0072B2",
      "UPPL1541" = "#D55E00"
    )
  )
  
  # Colors for the variable displayed on the x-axis.
  x_axis_colors <- group_colors[[x_axis_group]]
  
  # Colors for the variable represented by interaction lines.
  series_colors <- group_colors[[stratify_group]]
  
  # ---------------------------------------------------------------------------
  # Prepare the data
  # ---------------------------------------------------------------------------
  
  score_df <- score_df %>%
    mutate(
      
      # Convert the x-axis column to a factor and set its category order.
      "{x_axis_group}" := factor(
        .data[[x_axis_group]],
        levels = group_levels[[x_axis_group]]
      ),
      
      # Convert the series column to a factor and set its category order.
      "{stratify_group}" := factor(
        .data[[stratify_group]],
        levels = group_levels[[stratify_group]]
      )
    ) %>%
    
    # Remove rows missing information required for plotting.
    filter(
      !is.na(.data[[x_axis_group]]),
      !is.na(.data[[stratify_group]]),
      !is.na(.data[[score_col]])
    )
  
  # Stop if filtering removed all rows.
  if (nrow(score_df) == 0) {
    stop(
      "No complete rows remain after filtering the grouping and score columns."
    )
  }
  
  # ---------------------------------------------------------------------------
  # Create grouped dot plot
  # ---------------------------------------------------------------------------
  
  # The grouped dot plot:
  #
  #   - displays x_axis_group on the x-axis
  #   - colors points by x_axis_group
  #   - creates one panel for each stratify_group value
  #
  # Default orientation:
  #   x-axis: Parental and MP1
  #   panels: BBN963 and UPPL1541
  #
  # Swapped orientation:
  #   x-axis: BBN963 and UPPL1541
  #   panels: Parental and MP1
  grouped_dot_plot <- ggplot(
    score_df,
    aes(
      x = .data[[x_axis_group]],
      y = .data[[score_col]],
      color = .data[[x_axis_group]]
    )
  ) +
    
    # Add a reference line at a z-score of zero.
    geom_hline(
      yintercept = 0,
      linetype = "dotted",
      linewidth = 0.5
    ) +
    
    # Plot individual sample scores.
    geom_jitter(
      width = 0.08,
      height = 0,
      size = 2.8,
      alpha = 0.9
    ) +
    
    # Display the group mean as a thick horizontal marker.
    stat_summary(
      fun = mean,
      geom = "point",
      shape = 95,
      size = 10,
      linewidth = 1.2
    ) +
    
    # Create one panel for each series-group value.
    facet_wrap(
      vars(.data[[stratify_group]]),
      nrow = 1
    ) +
    
    # Apply the palette corresponding to the x-axis variable.
    scale_color_manual(
      values = x_axis_colors,
      drop = FALSE
    ) +
    ggplot2::coord_cartesian(
      ylim = c(-0.5, 0.5)
    ) +
    
    labs(
      title = plot_title,
      x = NULL,
      y = y_axis_label,
      color = group_labels[[x_axis_group]]
    ) +
    
    theme_classic(base_size = 12) +
    
    theme(
      plot.title = element_text(face = "bold"),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      legend.position = "right"
    )
  
  # ---------------------------------------------------------------------------
  # Calculate group means
  # ---------------------------------------------------------------------------
  
  # Calculate the mean for every combination of the two variables.
  #
  # These combinations are the same regardless of which variable is
  # displayed on the x-axis.
  mean_df <- score_df %>%
    group_by(
      .data[[stratify_group]],
      .data[[x_axis_group]]
    ) %>%
    summarise(
      mean_score = mean(
        .data[[score_col]],
        na.rm = TRUE
      ),
      .groups = "drop"
    )
  
  # ---------------------------------------------------------------------------
  # Create interaction plot
  # ---------------------------------------------------------------------------
  
  # The interaction plot:
  #
  #   - displays x_axis_group on the x-axis
  #   - creates one line for each stratify_group value
  #   - colors points and lines by stratify_group
  #
  # Default orientation:
  #   x-axis: Parental and MP1
  #   lines: BBN963 and UPPL1541
  #
  # Swapped orientation:
  #   x-axis: BBN963 and UPPL1541
  #   lines: Parental and MP1
  interaction_plot <- ggplot() +
    
    # Add a reference line at a z-score of zero.
    geom_hline(
      yintercept = 0,
      linetype = "dotted",
      linewidth = 0.5
    ) +
    
    # Display individual observations faintly in the background.
    geom_jitter(
      data = score_df,
      aes(
        x = .data[[x_axis_group]],
        y = .data[[score_col]],
        color = .data[[stratify_group]]
      ),
      width = 0.08,
      height = 0,
      size = 2.4,
      alpha = 0.45
    ) +
    
    # Connect the group means.
    geom_line(
      data = mean_df,
      aes(
        x = .data[[x_axis_group]],
        y = mean_score,
        group = .data[[stratify_group]],
        color = .data[[stratify_group]]
      ),
      linewidth = 1.1
    ) +
    
    # Add a point at each group mean.
    geom_point(
      data = mean_df,
      aes(
        x = .data[[x_axis_group]],
        y = mean_score,
        color = .data[[stratify_group]]
      ),
      size = 3.5
    ) +
    
    # Apply the palette corresponding to the series variable.
    scale_color_manual(
      values = series_colors,
      drop = FALSE
    ) +
    ggplot2::coord_cartesian(
      ylim = c(-0.5, 0.5)
    ) +
    
    labs(
      title = paste(
        plot_title,
        "\nInteraction View"
      ),
      x = NULL,
      y = y_axis_label,
      color = group_labels[[stratify_group]]
    ) +
    
    theme_classic(base_size = 12) +
    
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "right"
    )
  
  # ---------------------------------------------------------------------------
  # Save plots
  # ---------------------------------------------------------------------------
  
  if (save_plot) {
    
    if (is.null(out_file)) {
      
      # Convert the signature name into a filename-safe value.
      safe_label <- signature_label %>%
        str_replace_all("[^A-Za-z0-9]+", "_") %>%
        str_replace_all("^_|_$", "")
      
      # Include the orientation in the filename.
      #
      # This prevents the swapped version from overwriting plots
      # created using the default orientation.
      out_prefix <- paste0(
        safe_label,
        "_x_",
        x_axis_group,
        "_by_",
        stratify_group
      )
      
    } else {
      
      # Remove an extension supplied by the user.
      out_prefix <- tools::file_path_sans_ext(out_file)
    }
    
    out_dir <- here::here(
      "results",
      "plots"
    )
    
    dir.create(
      out_dir,
      recursive = TRUE,
      showWarnings = FALSE
    )
    
    ggsave(
      filename = file.path(
        out_dir,
        paste0(
          out_prefix,
          "_grouped_dot_plot.png"
        )
      ),
      plot = grouped_dot_plot,
      width = 6,
      height = 4,
      dpi = 300
    )
    
    ggsave(
      filename = file.path(
        out_dir,
        paste0(
          out_prefix,
          "_interaction_plot.png"
        )
      ),
      plot = interaction_plot,
      width = 5,
      height = 4,
      dpi = 300
    )
  }
  
  # Return the plots and calculated means.
  return(
    list(
      grouped_dot_plot = grouped_dot_plot,
      interaction_plot = interaction_plot,
      means = mean_df
    )
  )
}

# Compare BBN and UPPL models after selecting ONE condition,
# such as only Parental samples or only MP1 samples
# Ex. it'll show the signature scores for GATA3 between BBN and UPPL models for ONLY Parental samples

#ANSWERS: Among Parental samples, is signature activity different between the BBN and UPPL model families?

plot_signature_bbn_vs_uppl <- function(
    score_df, #df that contains signature score
    signature_label, #PPARG/GATA3
    parental_filter, #"Parental" or "MP1"
    plot_title = NULL, # optional custom title for plot
    out_file = NULL, #optional output file name
    save_plot = TRUE,
    model_cols = c( # assign colors for plots
      "BBN" = "#0072B2",
      "UPPL" = "#D55E00"
    ))
{
  
  # confirm signature label is valid
  if (!is.character(signature_label) ||
      length(signature_label) != 1 ||
      is.na(signature_label) ||
      signature_label == "") {
    stop("signature_label must be a single string.")
  }
  
  score_col <- paste0(signature_label, "_z_score")
  if (!score_col %in% colnames(score_df)) {
    stop(paste0("Column '", score_col, "' not found in score_df."))
  }
  
  if (!"cell_line" %in% colnames(score_df)) {
    stop("Column 'cell_line' not found in score_df.")
  }
  
  if (!"parental" %in% colnames(score_df)) {
    stop("Column 'parental' not found in score_df.")
  }
  
  if (is.null(plot_title)) {
    plot_title <- paste0(signature_label, " Activity: ", parental_filter, " Only")
  }
  
  # prepare data for plotting
  score_df <- score_df %>%
    #keep only selected condition
    filter(parental == parental_filter) %>% #ex. keep only parental samples
    mutate(
      # group cell lines into BBN or UPPL
      cell_model = case_when(
        grepl("^BBN", cell_line) ~ "BBN",
        grepl("^UPPL", cell_line) ~ "UPPL",
        TRUE ~ NA_character_
      ),
      cell_model = factor(cell_model, levels = c("BBN", "UPPL"))
    ) %>%
    filter(
      !is.na(cell_model),
      !is.na(.data[[score_col]])
    )
  
  # create BBN vs UPPL activity plot
  activity_plot <- ggplot(
    score_df,
    aes(
      x = cell_model,
      y = .data[[score_col]],
      color = cell_model
    )
  ) +
    # add horizontal reference line @z-score=0
    geom_hline(
      yintercept = 0,
      linetype = "dotted",
      linewidth = 0.5
    ) +
    #plot individual sample z scores
    geom_jitter(
      width = 0.08,
      height = 0,
      size = 2.8,
      alpha = 0.9
    ) +
    # add mean score as a thick horizontal marker
    stat_summary(
      fun = mean,
      geom = "point",
      shape = 95,
      size = 10,
      linewidth = 1.2
    ) +
    # apply colors for BBN and UPPL
    scale_color_manual(values = model_cols, drop = FALSE) +
    labs(
      title = plot_title,
      x = NULL,
      y = paste(signature_label, "Signature Z-Score"),
      color = NULL
    ) +
    ggplot2::coord_cartesian(
      ylim = c(-0.5, 0.5)
    ) +
    
    # plot title, axis labels, legend title
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "right"
    )
  # calculate mean signature score for BBN
  mean_df <- score_df %>%
    group_by(cell_model) %>%
    summarise(
      mean_score = mean(.data[[score_col]], na.rm = TRUE),
      .groups = "drop"
    )
  # save plot when save_plot=TRUE
  if (save_plot) {
    
    if (is.null(out_file)) {
      safe_label <- signature_label %>%
        str_replace_all("[^A-Za-z0-9]+", "_") %>%
        str_replace_all("^_|_$", "")
      
      safe_condition <- parental_filter %>%
        str_replace_all("[^A-Za-z0-9]+", "_") %>%
        str_replace_all("^_|_$", "")
      
      out_prefix <- paste0(safe_label, "_", safe_condition, "_bbn_vs_uppl")
    } else {
      out_prefix <- tools::file_path_sans_ext(out_file)
    }
    
    out_dir <- here::here("results", "plots")
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    
    ggsave(
      filename = file.path(out_dir, paste0(out_prefix, ".png")),
      plot = activity_plot,
      width = 4.5,
      height = 4,
      dpi = 400
    )
  }
  # return both plot and calculated group means
  return(list(
    activity_plot = activity_plot,
    means = mean_df
  ))
}


plot_signature_bbn_vs_uppl_combined <- function(
    score_df,                           # Data frame containing metadata and z-scores
    signature_label,                    # Signature name, such as "PPARG" or "GATA3"
    conditions = c("Parental", "MP1"),   # Conditions shown as separate subplots
    plot_title = NULL,                  # Optional custom plot title
    out_file = NULL,                    # Optional output filename
    save_plot = TRUE,                   # Whether to save the plot
    model_cols = c(
      "BBN" = "#0072B2",                # Color for BBN samples
      "UPPL" = "#D55E00"                # Color for UPPL samples
    )
) {
  
  # Confirm that signature_label is one non-empty string.
  if (!is.character(signature_label) ||
      length(signature_label) != 1 ||
      is.na(signature_label) ||
      signature_label == "") {
    stop("signature_label must be a single non-empty string.")
  }
  
  # Construct the expected score-column name.
  #
  # Example:
  #   signature_label = "PPARG"
  #   score_col = "PPARG_z_score"
  score_col <- paste0(signature_label, "_z_score")
  
  # Check that all required columns are present.
  required_cols <- c(
    "cell_line",
    "parental",
    score_col
  )
  
  missing_cols <- setdiff(
    required_cols,
    colnames(score_df)
  )
  
  if (length(missing_cols) > 0) {
    stop(
      paste0(
        "The following required columns are missing: ",
        paste(missing_cols, collapse = ", ")
      )
    )
  }
  
  # Confirm that only supported conditions were supplied.
  valid_conditions <- c(
    "Parental",
    "MP1"
  )
  
  if (!all(conditions %in% valid_conditions)) {
    stop(
      paste0(
        "conditions can only contain: ",
        paste(valid_conditions, collapse = ", ")
      )
    )
  }
  
  # Remove duplicated conditions while preserving their order.
  conditions <- unique(conditions)
  
  # Check that model_cols contains colors for both expected models.
  required_models <- c(
    "BBN",
    "UPPL"
  )
  
  missing_model_colors <- setdiff(
    required_models,
    names(model_cols)
  )
  
  if (length(missing_model_colors) > 0) {
    stop(
      paste0(
        "model_cols is missing colors for: ",
        paste(missing_model_colors, collapse = ", ")
      )
    )
  }
  
  # Create a default title when none is supplied.
  if (is.null(plot_title)) {
    plot_title <- paste(
      signature_label,
      "Activity: BBN vs UPPL"
    )
  }
  
  # Create the y-axis label.
  y_axis_label <- paste(
    signature_label,
    "Signature Z-Score"
  )
  
  # Prepare the data for plotting.
  plot_df <- score_df %>%
    dplyr::mutate(
      
      # Group individual cell lines into BBN or UPPL model families.
      #
      # Examples:
      #   BBN963   becomes BBN
      #   UPPL1541 becomes UPPL
      cell_model = dplyr::case_when(
        grepl("^BBN", cell_line) ~ "BBN",
        grepl("^UPPL", cell_line) ~ "UPPL",
        TRUE ~ NA_character_
      ),
      
      # Set the display order for BBN and UPPL.
      cell_model = factor(
        cell_model,
        levels = c("BBN", "UPPL")
      ),
      
      # Set the order of the Parental and MP1 subplot panels.
      parental = factor(
        parental,
        levels = conditions
      )
    ) %>%
    
    # Keep only the requested conditions and complete observations.
    dplyr::filter(
      parental %in% conditions,
      !is.na(parental),
      !is.na(cell_model),
      !is.na(.data[[score_col]])
    )
  
  # Stop if filtering removed all observations.
  if (nrow(plot_df) == 0) {
    stop("No valid observations remain after filtering.")
  }
  
  # Create a plot comparing BBN and UPPL within each condition.
  activity_plot <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(
      x = cell_model,
      y = .data[[score_col]],
      color = cell_model
    )
  ) +
    
    # Add a horizontal reference line at a z-score of zero.
    ggplot2::geom_hline(
      yintercept = 0,
      linetype = "dotted",
      linewidth = 0.5
    ) +
    
    # Plot individual sample scores.
    
    # Horizontal jitter prevents overlapping points.
    ggplot2::geom_jitter(
      width = 0.08,
      height = 0,
      size = 2.8,
      alpha = 0.9
    ) +
    
    # Add the mean as a thick horizontal marker.
    ggplot2::stat_summary(
      fun = mean,
      geom = "point",
      shape = 95,
      size = 10,
      linewidth = 1.2
    ) +
    
    # Create separate Parental and MP1 subplots.
    ggplot2::facet_wrap(
      ggplot2::vars(parental),
      nrow = 1
    ) +
    
    # Apply the BBN and UPPL colors.
    ggplot2::scale_color_manual(
      values = model_cols,
      drop = FALSE
    ) +
    ggplot2::coord_cartesian(
      ylim = c(-0.5, 0.5)
    ) +
    
    # Add plot labels.
    ggplot2::labs(
      title = plot_title,
      x = NULL,
      y = y_axis_label,
      color = NULL
    ) +
    
    ggplot2::theme_classic(base_size = 12) +
    
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold"
      ),
      
      # Add horizontal space between the Parental and MP1 subplots.
      panel.spacing.x = grid::unit(
        1.5,
        "lines"
      ),
      
      # Draw a border around each subplot.
      panel.border = ggplot2::element_rect(
        color = "grey60",
        fill = NA,
        linewidth = 0.7
      ),
      
      # Add a visible background and border to each subplot title.
      strip.background = ggplot2::element_rect(
        fill = "grey95",
        color = "grey60",
        linewidth = 0.7
      ),
      
      # Format the Parental and MP1 subplot labels.
      strip.text = ggplot2::element_text(
        face = "bold",
        size = 12,
        margin = ggplot2::margin(
          t = 5,
          r = 5,
          b = 5,
          l = 5
        )
      ),
      
      legend.position = "right"
    )
  
  # Calculate the mean score for every condition and model combination.
  #
  # This calculates means for:
  #   Parental BBN
  #   Parental UPPL
  #   MP1 BBN
  #   MP1 UPPL
  mean_df <- plot_df %>%
    dplyr::group_by(
      parental,
      cell_model
    ) %>%
    dplyr::summarise(
      mean_score = mean(
        .data[[score_col]],
        na.rm = TRUE
      ),
      .groups = "drop"
    )
  
  # Save the combined plot when save_plot is TRUE.
  if (save_plot) {
    
    if (is.null(out_file)) {
      
      # Convert the signature label into a filename-safe value.
      safe_label <- signature_label %>%
        stringr::str_replace_all(
          "[^A-Za-z0-9]+",
          "_"
        ) %>%
        stringr::str_replace_all(
          "^_|_$",
          ""
        )
      
      out_prefix <- paste0(
        safe_label,
        "_Parental_MP1_bbn_vs_uppl"
      )
      
    } else {
      
      # Remove any extension from the supplied filename.
      out_prefix <- tools::file_path_sans_ext(
        out_file
      )
    }
    
    # Define and create the output directory.
    out_dir <- here::here(
      "results",
      "plots"
    )
    
    dir.create(
      out_dir,
      recursive = TRUE,
      showWarnings = FALSE
    )
    
    # Save the plot as a high-resolution PNG.
    ggplot2::ggsave(
      filename = file.path(
        out_dir,
        paste0(out_prefix, ".png")
      ),
      plot = activity_plot,
      width = 7,
      height = 4,
      dpi = 400
    )
  }
  
  # Return the combined plot and calculated means.
  return(
    list(
      activity_plot = activity_plot,
      means = mean_df
    )
  )
}

# Compare Parental and MP1 models after selecting ONE condition,
# such as only BBN samples or only UPPL samples
# Ex. it'll show the signature scores for GATA3 between Parental and MP1 models for ONLY BBN samples

#ANSWERS: Within BBN subtype, is signature activity different between Parental and MP1 samples?

plot_signature_parental_vs_mp1 <- function(score_df, # input df from z-scores
                                           signature_label, #gata3/pparg
                                           cell_model_filter = "BBN", #subtype: BBN/UPPL 
                                           x_group = "parental",
                                           color_group = "parental",
                                           plot_title = NULL, #custom plot title
                                           out_file = NULL, # saves the plots
                                           save_plot = TRUE,
                                           condition_cols = c(
                                             "Parental" = "#556B2F", # color for parental
                                             "MP1" = "#D8A3B0" #color for MP1
                                           )) {
  
  if (!is.character(signature_label) ||
      length(signature_label) != 1 ||
      is.na(signature_label) ||
      signature_label == "") {
    stop("signature_label must be a single string.")
  }
  
  score_col <- paste0(signature_label, "_z_score")
  
  if (!score_col %in% colnames(score_df)) {
    stop(paste0("Column '", score_col, "' not found in score_df."))
  }
  
  if (!"cell_line" %in% colnames(score_df)) {
    stop("Column 'cell_line' not found in score_df.")
  }
  
  if (!"parental" %in% colnames(score_df)) {
    stop("Column 'parental' not found in score_df.")
  }
  
  if (is.null(plot_title)) {
    plot_title <- paste0(signature_label, " Activity: ", cell_model_filter, " Only")
  }
  y_axis_label <- paste(signature_label, "Signature Z-Score")
  
  score_df <- score_df %>%
    mutate(
      cell_model = case_when(
        grepl("^BBN", cell_line) ~ "BBN",
        grepl("^UPPL", cell_line) ~ "UPPL",
        TRUE ~ NA_character_
      ),
      parental = factor(parental, levels = c("Parental", "MP1"))
    ) %>%
    filter(
      cell_model == cell_model_filter,
      !is.na(parental),
      !is.na(.data[[score_col]])
    )
  
  activity_plot <- ggplot(
    score_df,
    aes(
      x = .data[[x_group]],
      y = .data[[score_col]],
      color = .data[[color_group]]
    )
  ) +
    geom_hline(
      yintercept = 0,
      linetype = "dotted",
      linewidth = 0.5
    ) +
    geom_jitter(
      width = 0.08,
      height = 0,
      size = 2.8,
      alpha = 0.9
    ) +
    stat_summary(
      fun = mean,
      geom = "point",
      shape = 95,
      size = 10,
      linewidth = 1.2
    ) +
    scale_color_manual(values = condition_cols, drop = FALSE) +
    labs(
      title = plot_title,
      x = NULL,
      y = y_axis_label,
      color = NULL
    ) +
    ggplot2::coord_cartesian(
      ylim = c(-0.5, 0.5)
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "right"
    )
  
  mean_df <- score_df %>%
    group_by(.data[[x_group]]) %>%
    summarise(
      mean_score = mean(.data[[score_col]], na.rm = TRUE),
      .groups = "drop"
    )
  
  if (save_plot) {
    
    if (is.null(out_file)) {
      safe_label <- signature_label %>%
        str_replace_all("[^A-Za-z0-9]+", "_") %>%
        str_replace_all("^_|_$", "")
      
      out_prefix <- paste0(safe_label, "_", cell_model_filter, "_parental_vs_mp1")
    } else {
      out_prefix <- tools::file_path_sans_ext(out_file)
    }
    
    out_dir <- here::here("results", "plots")
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    
    ggsave(
      filename = file.path(out_dir, paste0(out_prefix, ".png")),
      plot = activity_plot,
      width = 4.5,
      height = 4,
      dpi = 300
    )
  }
  
  return(list(
    activity_plot = activity_plot,
    means = mean_df
  ))
}

plot_signature_parental_vs_mp1_combined <- function(
    score_df,                          # Data frame containing metadata and z-scores
    signature_label,                   # Signature name, such as "GATA3" or "PPARG"
    cell_models = c("BBN", "UPPL"),     # Models shown as separate subplots
    plot_title = NULL,                 # Optional custom plot title
    out_file = NULL,                   # Optional output filename
    save_plot = TRUE,                  # Whether to save the plot
    condition_cols = c(
      "Parental" = "#556B2F",
      "MP1" = "#D8A3B0"
    )
) {
  
  # Check that signature_label is one non-empty string.
  if (!is.character(signature_label) ||
      length(signature_label) != 1 ||
      is.na(signature_label) ||
      signature_label == "") {
    stop("signature_label must be a single non-empty string.")
  }
  
  # Construct the expected score-column name.
  #
  # Example:
  #   signature_label = "PPARG"
  #   score_col = "PPARG_z_score"
  score_col <- paste0(signature_label, "_z_score")
  
  # Check that all required columns are present.
  required_cols <- c(
    "cell_line",
    "parental",
    score_col
  )
  
  missing_cols <- setdiff(
    required_cols,
    colnames(score_df)
  )
  
  if (length(missing_cols) > 0) {
    stop(
      paste0(
        "The following required columns are missing: ",
        paste(missing_cols, collapse = ", ")
      )
    )
  }
  
  # Check that cell_models contains only supported model names.
  valid_models <- c("BBN", "UPPL")
  
  if (!all(cell_models %in% valid_models)) {
    stop(
      paste0(
        "cell_models can only contain: ",
        paste(valid_models, collapse = ", ")
      )
    )
  }
  
  # Remove duplicated model names while preserving their order.
  cell_models <- unique(cell_models)
  
  # Create a default plot title.
  if (is.null(plot_title)) {
    plot_title <- paste(
      signature_label,
      "Activity: Parental vs MP1"
    )
  }
  
  # Create the y-axis label.
  y_axis_label <- paste(
    signature_label,
    "Signature Z-Score"
  )
  
  # Prepare the data for plotting.
  plot_df <- score_df %>%
    dplyr::mutate(
      
      # Assign each cell line to the BBN or UPPL model family.
      cell_model = dplyr::case_when(
        grepl("^BBN", cell_line) ~ "BBN",
        grepl("^UPPL", cell_line) ~ "UPPL",
        TRUE ~ NA_character_
      ),
      
      # Display Parental before MP1.
      parental = factor(
        parental,
        levels = c("Parental", "MP1")
      ),
      
      # Control the order of the subplot panels.
      cell_model = factor(
        cell_model,
        levels = cell_models
      )
    ) %>%
    
    # Keep only the requested models and complete observations.
    dplyr::filter(
      cell_model %in% cell_models,
      !is.na(cell_model),
      !is.na(parental),
      !is.na(.data[[score_col]])
    )
  
  # Stop if no valid observations remain.
  if (nrow(plot_df) == 0) {
    stop("No valid observations remain after filtering.")
  }
  
  # Create one subplot for each cell model.
  activity_plot <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(
      x = parental,
      y = .data[[score_col]],
      color = parental
    )
  ) +
    
    # Add a reference line at a z-score of zero.
    ggplot2::geom_hline(
      yintercept = 0,
      linetype = "dotted",
      linewidth = 0.5
    ) +
    
    # Plot individual sample scores.
    ggplot2::geom_jitter(
      width = 0.08,
      height = 0,
      size = 2.8,
      alpha = 0.9
    ) +
    
    # Add the mean as a thick horizontal marker.
    ggplot2::stat_summary(
      fun = mean,
      geom = "point",
      shape = 95,
      size = 10,
      linewidth = 1.2
    ) +
    
    # Create separate BBN and UPPL subplots.
    ggplot2::facet_wrap(
      ggplot2::vars(cell_model),
      nrow = 1
    ) +
    
    # Apply the Parental and MP1 colors.
    ggplot2::scale_color_manual(
      values = condition_cols,
      drop = FALSE
    ) +
    ggplot2::coord_cartesian(
      ylim = c(-0.5, 0.5)
    ) +
    
    # Add plot labels.
    ggplot2::labs(
      title = plot_title,
      x = NULL,
      y = y_axis_label,
      color = NULL
    ) +
    
    ggplot2::theme_classic(base_size = 12) +
    
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold"
      ),
      
      # Add more horizontal space between subplots.
      panel.spacing.x = grid::unit(
        1.5,
        "lines"
      ),
      
      # Draw a border around each subplot.
      panel.border = ggplot2::element_rect(
        color = "grey60",
        fill = NA,
        linewidth = 0.7
      ),
      
      # Add a visible background and border to each subplot title.
      strip.background = ggplot2::element_rect(
        fill = "grey95",
        color = "grey60",
        linewidth = 0.7
      ),
      
      # Format the BBN and UPPL subplot labels.
      strip.text = ggplot2::element_text(
        face = "bold",
        size = 12,
        margin = ggplot2::margin(
          t = 5,
          r = 5,
          b = 5,
          l = 5
        )
      ),
      
      legend.position = "right"
    )
  
  # Calculate the mean score for each model and condition.
  mean_df <- plot_df %>%
    dplyr::group_by(
      cell_model,
      parental
    ) %>%
    dplyr::summarise(
      mean_score = mean(
        .data[[score_col]],
        na.rm = TRUE
      ),
      .groups = "drop"
    )
  
  # Save the combined faceted plot.
  if (save_plot) {
    
    if (is.null(out_file)) {
      
      # Make the signature label safe for use in a filename.
      safe_label <- signature_label %>%
        stringr::str_replace_all(
          "[^A-Za-z0-9]+",
          "_"
        ) %>%
        stringr::str_replace_all(
          "^_|_$",
          ""
        )
      
      out_prefix <- paste0(
        safe_label,
        "_BBN_UPPL_parental_vs_mp1"
      )
      
    } else {
      
      # Remove any supplied file extension.
      out_prefix <- tools::file_path_sans_ext(
        out_file
      )
    }
    
    # Define and create the output directory.
    out_dir <- here::here(
      "results",
      "plots"
    )
    
    dir.create(
      out_dir,
      recursive = TRUE,
      showWarnings = FALSE
    )
    
    # Save the plot as a high-resolution PNG.
    ggplot2::ggsave(
      filename = file.path(
        out_dir,
        paste0(out_prefix, ".png")
      ),
      plot = activity_plot,
      width = 7,
      height = 4,
      dpi = 300
    )
  }
  
  # Return the plot and calculated means.
  return(
    list(
      activity_plot = activity_plot,
      means = mean_df
    )
  )
}

plot_signature_bbn_vs_uppl_pooled <- function(
    score_df,
    signature_label,
    plot_title = NULL,
    out_file = NULL,
    save_plot = TRUE,
    model_cols = c(
      "BBN" = "#1E3A8A",
      "UPPL" = "#F5A623"
    )
) {

  if (!is.character(signature_label) ||
      length(signature_label) != 1 ||
      is.na(signature_label) ||
      signature_label == "") {
    stop("signature_label must be a single non-empty string.")
  }

  score_col <- paste0(signature_label, "_z_score")

  required_cols <- c(
    "cell_line",
    score_col
  )

  missing_cols <- setdiff(
    required_cols,
    colnames(score_df)
  )

  if (length(missing_cols) > 0) {
    stop(
      paste0(
        "The following required columns are missing: ",
        paste(missing_cols, collapse = ", ")
      )
    )
  }

  required_models <- c("BBN", "UPPL")

  missing_model_colors <- setdiff(
    required_models,
    names(model_cols)
  )

  if (length(missing_model_colors) > 0) {
    stop(
      paste0(
        "model_cols is missing colors for: ",
        paste(missing_model_colors, collapse = ", ")
      )
    )
  }

  if (is.null(plot_title)) {
    plot_title <- paste(
      signature_label,
      "Activity: BBN vs UPPL"
    )
  }

  y_axis_label <- paste(
    signature_label,
    "Signature Z-Score"
  )

  plot_df <- score_df %>%
    dplyr::mutate(
      cell_model = dplyr::case_when(
        grepl("^BBN", cell_line) ~ "BBN",
        grepl("^UPPL", cell_line) ~ "UPPL",
        TRUE ~ NA_character_
      ),
      cell_model = factor(
        cell_model,
        levels = c("BBN", "UPPL")
      )
    ) %>%
    dplyr::filter(
      !is.na(cell_model),
      !is.na(.data[[score_col]])
    )

  if (nrow(plot_df) == 0) {
    stop("No valid BBN or UPPL observations remain after filtering.")
  }

  activity_plot <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(
      x = cell_model,
      y = .data[[score_col]],
      color = cell_model
    )
  ) +
    ggplot2::geom_hline(
      yintercept = 0,
      linetype = "dotted",
      linewidth = 0.5
    ) +
    ggplot2::geom_jitter(
      width = 0.08,
      height = 0,
      size = 2.8,
      alpha = 0.9
    ) +
    ggplot2::stat_summary(
      fun = mean,
      geom = "point",
      shape = 95,
      size = 10,
      linewidth = 1.2
    ) +
    ggplot2::scale_color_manual(
      values = model_cols,
      drop = FALSE
    ) +
    ggplot2::coord_cartesian(
      ylim = c(-0.5, 0.5)
    ) +
    ggplot2::labs(
      title = stringr::str_wrap(plot_title, width = 35),
      x = NULL,
      y = y_axis_label,
      color = NULL
    ) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold",
        hjust = 0.5
      ),
      panel.border = ggplot2::element_rect(
        color = "grey60",
        fill = NA,
        linewidth = 0.7
      ),
      legend.position = "right"
    )

  mean_df <- plot_df %>%
    dplyr::group_by(cell_model) %>%
    dplyr::summarise(
      mean_score = mean(
        .data[[score_col]],
        na.rm = TRUE
      ),
      n = dplyr::n(),
      .groups = "drop"
    )

  if (save_plot) {

    if (is.null(out_file)) {
      safe_label <- signature_label %>%
        stringr::str_replace_all("[^A-Za-z0-9]+", "_") %>%
        stringr::str_replace_all("^_|_$", "")

      out_prefix <- paste0(
        safe_label,
        "_pooled_bbn_vs_uppl"
      )
    } else {
      out_prefix <- tools::file_path_sans_ext(out_file)
    }

    out_dir <- here::here(
      "results",
      "plots"
    )

    dir.create(
      out_dir,
      recursive = TRUE,
      showWarnings = FALSE
    )

    ggplot2::ggsave(
      filename = file.path(
        out_dir,
        paste0(out_prefix, ".png")
      ),
      plot = activity_plot,
      width = 5,
      height = 4,
      dpi = 400
    )
  }

  list(
    activity_plot = activity_plot,
    means = mean_df
  )
}




plot_signature_parental_vs_mp1_pooled <- function(
    score_df,
    signature_label,
    plot_title = NULL,
    out_file = NULL,
    save_plot = TRUE,
    condition_cols = c(
      "Parental" = "#556B2F",
      "MP1" = "#D8A3B0"
    )
) {
  
  if (!is.character(signature_label) ||
      length(signature_label) != 1 ||
      is.na(signature_label) ||
      signature_label == "") {
    stop("signature_label must be a single non-empty string.")
  }
  
  score_col <- paste0(signature_label, "_z_score")
  
  required_cols <- c(
    "parental",
    score_col
  )
  
  missing_cols <- setdiff(
    required_cols,
    colnames(score_df)
  )
  
  if (length(missing_cols) > 0) {
    stop(
      paste0(
        "The following required columns are missing: ",
        paste(missing_cols, collapse = ", ")
      )
    )
  }
  
  required_conditions <- c("Parental", "MP1")
  
  missing_condition_colors <- setdiff(
    required_conditions,
    names(condition_cols)
  )
  
  if (length(missing_condition_colors) > 0) {
    stop(
      paste0(
        "condition_cols is missing colors for: ",
        paste(missing_condition_colors, collapse = ", ")
      )
    )
  }
  
  if (is.null(plot_title)) {
    plot_title <- paste(
      signature_label,
      "Signature Activity by Condition"
    )
  }
  
  y_axis_label <- paste(
    signature_label,
    "Signature Z-Score"
  )
  
  plot_df <- score_df %>%
    dplyr::mutate(
      parental = factor(
        parental,
        levels = c("Parental", "MP1")
      )
    ) %>%
    dplyr::filter(
      !is.na(parental),
      !is.na(.data[[score_col]])
    )
  
  if (nrow(plot_df) == 0) {
    stop("No valid Parental or MP1 observations remain after filtering.")
  }
  
  activity_plot <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(
      x = parental,
      y = .data[[score_col]],
      color = parental
    )
  ) +
    ggplot2::geom_hline(
      yintercept = 0,
      linetype = "dotted",
      linewidth = 0.5
    ) +
    ggplot2::geom_jitter(
      width = 0.05,
      height = 0,
      size = 2.8,
      alpha = 0.9
    ) +
    ggplot2::stat_summary(
      fun = mean,
      geom = "point",
      shape = 95,
      size = 10,
      linewidth = 1.2
    ) +
    ggplot2::scale_color_manual(
      values = condition_cols,
      drop = FALSE
    ) + 
    ggplot2::coord_cartesian(
      ylim = c(-0.5, 0.5)
    ) +
    ggplot2::labs(
      title = stringr::str_wrap(
        plot_title,
        width = 35
      ),
      x = NULL,
      y = y_axis_label,
      color = NULL
    ) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold",
        hjust = 0.5
      ),
      panel.border = ggplot2::element_rect(
        color = "grey60",
        fill = NA,
        linewidth = 0.7
      ),
      legend.position = "right"
    )
  
  mean_df <- plot_df %>%
    dplyr::group_by(parental) %>%
    dplyr::summarise(
      mean_score = mean(
        .data[[score_col]],
        na.rm = TRUE
      ),
      n = dplyr::n(),
      .groups = "drop"
    )
  
  if (save_plot) {
    
    if (is.null(out_file)) {
      safe_label <- signature_label %>%
        stringr::str_replace_all("[^A-Za-z0-9]+", "_") %>%
        stringr::str_replace_all("^_|_$", "")
      
      out_prefix <- paste0(
        safe_label,
        "_pooled_parental_vs_mp1"
      )
    } else {
      out_prefix <- tools::file_path_sans_ext(out_file)
    }
    
    out_dir <- here::here(
      "results",
      "plots"
    )
    
    dir.create(
      out_dir,
      recursive = TRUE,
      showWarnings = FALSE
    )
    
    ggplot2::ggsave(
      filename = file.path(
        out_dir,
        paste0(out_prefix, ".png")
      ),
      plot = activity_plot,
      width = 5,
      height = 6,
      dpi = 400
    )
  }
  
  list(
    activity_plot = activity_plot,
    means = mean_df
  )
}