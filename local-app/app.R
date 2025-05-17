# app.R

# --- 0. Load Libraries and Source Functions ---
library(shiny)
library(ggplot2)
library(dplyr)
library(patchwork) # For your plotting helper

# Source your DDM simulation functions and plotting helper
source("R/02_ddm_simulator_basic.R") # If you use simulate_diffusion_trial_with_path from here
source("R/03_ddm_simulator_variable.R")
source("R/utils/plot_basic_ddm_path_rt.R") # Or R/plot_ddm_paths_histograms.R

# --- 1. UI Definition ---
ui <- fluidPage(
  titlePanel("Interactive DDM Explorer with Parameter Variability"),
  
  sidebarLayout(
    sidebarPanel(
      h4("Simulation Settings"),
      numericInput("n_trials_shiny", "Number of Trials:", min = 1, max = 5000, value = 1000, step = 100),
      sliderInput("dt_shiny", "Time Step (dt):", value = 0.001, min = 0.0001, max = 0.01, step = 0.0001),
      helpText("Smaller dt is more accurate but slower."),
      actionButton("run_simulation_shiny", "Run Simulation", icon = icon("play-circle"), class="btn-success"),
      hr(),
      
      h4("Mean DDM Parameters"),
      numericInput("mean_v_shiny", "Mean Drift Rate (mean_v):", value = 0.15, step = 0.01),
      numericInput("a_shiny", "Threshold (a):", value = 1.0, min = 0.1, step = 0.05),
      sliderInput("mean_z_shiny_slider", "Mean Start Point (mean_z / a):", min = 0.01, max = 0.99, value = 0.5, step = 0.01), # Relative to 'a'
      # Hidden numeric input to store actual mean_z based on slider and 'a'
      # shiny::conditionalPanel(condition = "false", numericInput("mean_z_shiny", "", value=0.5)),
      numericInput("s_shiny", "Within-Trial Noise (s):", value = 0.1, min = 0.01, step = 0.01),
      numericInput("mean_ter_shiny", "Mean Non-Decision Time (mean_ter):", value = 0.2, min = 0, step = 0.01),
      hr(),
      
      h4("Across-Trial Variability Parameters"),
      numericInput("sv_shiny", "Drift Rate SD (sv):", value = 0.1, min = 0, step = 0.01),
      numericInput("sz_shiny", "Start Point Range (sz):", value = 0.2, min = 0, step = 0.01),
      numericInput("st0_shiny", "Non-Decision Time Range (st0):", value = 0.05, min = 0, step = 0.01),
      hr(),
      
      h4("Plotting Options"),
      numericInput("plot_max_time_shiny", "Max Time for Plots (s):", value = NA, min = 0.1, step = 0.1),
      helpText("Leave NA for automatic scaling."),
      numericInput("hist_binwidth_shiny", "Histogram Binwidth (s):", value = 0.04, min = 0.001, step = 0.001)
      
    ), # end sidebarPanel
    
    mainPanel(
      tabsetPanel(
        tabPanel("Paths & RT Histograms",
                 plotOutput("ddm_paths_plot_shiny", height = "700px") # Main combined plot
        ),
        tabPanel("Sampled Parameters",
                 plotOutput("sampled_v_hist_shiny", height = "250px"),
                 plotOutput("sampled_z_hist_shiny", height = "250px"),
                 plotOutput("sampled_ter_hist_shiny", height = "250px")
        ),
        tabPanel("Summary Statistics",
                 verbatimTextOutput("choice_prop_shiny"),
                 tableOutput("rt_summary_shiny")
        )
      ) # end tabsetPanel
    ) # end mainPanel
  ) # end sidebarLayout
) # end fluidPage

# --- 2. Server Logic ---
server <- function(input, output, session) {
  
  # Reactive value to store simulation results
  sim_results <- reactiveVal(NULL)
  
  # Calculate actual mean_z based on slider and 'a'
  # This is a common pattern for relative parameters
  actual_mean_z <- reactive({
    req(input$a_shiny) # Ensure a_shiny is available
    input$mean_z_shiny_slider * input$a_shiny
  })
  
  # Run simulation when button is clicked
  observeEvent(input$run_simulation_shiny, {
    # Show a notification that simulation is running
    showNotification("Running DDM simulation...", type = "message", duration = NULL, id="simNotify")
    
    # Use actual_mean_z() which is reactive
    current_mean_z <- actual_mean_z()
    
    # Validate sz relative to mean_z and a to prevent z_trial always being boundary
    max_sz_allowed <- min(current_mean_z * 2, (input$a_shiny - current_mean_z) * 2) * 0.98 # ensure it's within bounds
    current_sz <- min(input$sz_shiny, max_sz_allowed)
    if(input$sz_shiny > max_sz_allowed && input$sz_shiny > 0) {
      showNotification(paste0("sz reduced to ", round(current_sz,3) ," to keep z within (0,a)."), type="warning", duration=5)
    }
    
    
    # Simulate all trials (list of lists)
    n_trials <- input$n_trials_shiny
    trials_list_output <- vector("list", n_trials)
    
    # Using a progress bar for longer simulations
    withProgress(message = 'Simulating trials', value = 0, {
      for (i in 1:n_trials) {
        trials_list_output[[i]] <- simulate_diffusion_trial_variable_with_path( # <--- CALL NEW FUNCTION
          mean_v = input$mean_v_shiny,
          a = input$a_shiny,
          mean_z = current_mean_z,
          s = input$s_shiny,
          mean_ter = input$mean_ter_shiny,
          sv = input$sv_shiny,
          sz = current_sz,
          st0 = input$st0_shiny,
          dt = input$dt_shiny,
          max_decision_time = 7.0 
        )
        incProgress(1/n_trials, detail = paste("Trial", i))
      }
    })
    
    sim_results(list(
      trials_list = trials_list_output, 
      params = list( # Store params used for plotting
        v = input$mean_v_shiny, 
        a = input$a_shiny, 
        z = current_mean_z,
        max_time = if(is.na(input$plot_max_time_shiny)) NULL else input$plot_max_time_shiny,
        binwidth = input$hist_binwidth_shiny
      )
    ))
    removeNotification(id="simNotify") # Remove notification
  })
  
  # --- Output: Paths & RT Histograms Plot ---
  output$ddm_paths_plot_shiny <- renderPlot({
    res <- sim_results()
    req(res) # Require results to be available
    
    # NO aggressive filtering here for valid_trials_for_plot if we want to see NA paths.
    # The plotting function will handle what to do with NAs for paths vs. histograms.
    # We just need to ensure the list itself and the path_data within are valid.
    
    # Basic check: Ensure trials_list exists and has elements
    if (is.null(res$trials_list) || length(res$trials_list) == 0) {
      return(ggplot() + labs(title = "No simulation data available.") + theme_void())
    }
    
    # Optional: Filter out trials where path_data itself is NULL or empty,
    # as these cannot be plotted at all. This is different from an NA choice.
    trials_to_plot <- Filter(function(x) !is.null(x$path_data) && nrow(x$path_data) > 0, res$trials_list)
    
    if (length(trials_to_plot) < 1) { # Need at least one path to try plotting
      return(ggplot() + labs(title = "No trials with valid path data to plot.") + theme_void())
    }
    
    plot_ddm_paths_with_histograms(
      trials_data_list = trials_to_plot, # Pass the (potentially NA-choice-inclusive) list
      v_drift = res$params$v,
      z_start = res$params$z,
      a_threshold = res$params$a,
      max_time_to_plot = res$params$max_time,
      hist_binwidth = res$params$binwidth,
      main_plot_title = paste(length(trials_to_plot), "DDM Trials: Paths & RTs")
    )
  })
  
  # --- Output: Sampled Parameter Histograms ---
  # Helper function to create parameter histograms
  create_param_hist <- function(data_vector, param_name, mean_param_val) {
    if (length(na.omit(data_vector)) < 2) return(ggplot() + labs(title=paste("Not enough data for", param_name)) + theme_void())
    df <- data.frame(value = data_vector)
    ggplot(df, aes(x = value)) +
      geom_histogram(aes(y = ..density..), fill = "grey70", color = "black", bins = 30, alpha = 0.7) +
      geom_density(color = "blue", linewidth = 1) +
      geom_vline(xintercept = mean_param_val, color = "red", linetype = "dashed", linewidth = 1) +
      labs(title = paste("Sampled Trial", param_name), x = "Value", y = "Density") +
      theme_minimal()
  }
  
  output$sampled_v_hist_shiny <- renderPlot({
    res <- sim_results()
    req(res)
    v_trials <- sapply(res$trials_list, function(x) x$v_trial)
    create_param_hist(v_trials, "Drift Rates (v_trial)", input$mean_v_shiny)
  })
  
  output$sampled_z_hist_shiny <- renderPlot({
    res <- sim_results()
    req(res)
    z_trials <- sapply(res$trials_list, function(x) x$z_trial)
    create_param_hist(z_trials, "Start Points (z_trial)", actual_mean_z()) # Use reactive actual_mean_z
  })
  
  output$sampled_ter_hist_shiny <- renderPlot({
    res <- sim_results()
    req(res)
    ter_trials <- sapply(res$trials_list, function(x) x$ter_trial)
    create_param_hist(ter_trials, "Non-Decision Times (ter_trial)", input$mean_ter_shiny)
  })
  
  # --- Output: Summary Statistics ---
  # Reactive expression for full data frame from simulations for summary stats
  full_sim_df <- reactive({
    res <- sim_results()
    req(res) # Ensures res and res$trials_list are not NULL
    
    # Filter out any strictly NULL elements from the list of trials first
    # (e.g., if something went wrong during simulation for a particular iteration, though unlikely)
    valid_list_elements <- Filter(Negate(is.null), res$trials_list)
    
    if (length(valid_list_elements) == 0) {
      return(data.frame()) # Return empty df if no valid trial data
    }
    
    # Define expected scalar columns to extract
    expected_scalar_cols <- c("choice", "rt", "decision_time", "v_trial", "z_trial", "ter_trial")
    
    # Check if the first valid element has these expected columns
    # This assumes all elements have a similar structure
    if (length(valid_list_elements) > 0 && 
        all(expected_scalar_cols %in% names(valid_list_elements[[1]]))) {
      
      choice_vec        <- sapply(valid_list_elements, function(x) x$choice)
      rt_vec            <- sapply(valid_list_elements, function(x) x$rt)
      decision_time_vec <- sapply(valid_list_elements, function(x) x$decision_time)
      v_trial_vec       <- sapply(valid_list_elements, function(x) x$v_trial)
      z_trial_vec       <- sapply(valid_list_elements, function(x) x$z_trial)
      ter_trial_vec     <- sapply(valid_list_elements, function(x) x$ter_trial)
      
      temp_df <- data.frame(
        choice = choice_vec,
        rt = rt_vec,
        decision_time = decision_time_vec,
        v_trial = v_trial_vec,
        z_trial = z_trial_vec,
        ter_trial = ter_trial_vec,
        stringsAsFactors = FALSE # Good practice
      )
      
    } else {
      # If expected columns are missing, return an empty data frame to avoid errors downstream
      warning("Expected columns for summary statistics are missing from simulation output.")
      return(data.frame())
    }
    
    return(temp_df)
  })
  
  # Choice proportions: output$choice_prop_shiny
  output$choice_prop_shiny <- renderPrint({
    df <- full_sim_df()
    req(df, nrow(df) > 0)
    
    choice_table <- table(df$choice)
    choice_prop <- prop.table(choice_table)
    
    cat("Choice Proportions:\n")
    choices <- names(choice_table)
    for (ch in choices) {
      cat(paste0(ch, ": ", round(100*choice_prop[ch], 1), "%"), "\n")
    }
    cat("\nTotal Trials:", nrow(df), "\n")
    
    # Show proportion of timeouts if present
    if ("NA" %in% choices) {
      cat("Timeouts:", sum(is.na(df$choice)), "\n")
    }
  })
  
  # RT summary statistics: output$rt_summary_shiny
  output$rt_summary_shiny <- renderTable({
    df <- full_sim_df()
    req(df, nrow(df) > 0)
    
    # Filter out any NAs in rt (shouldn't happen if simulation is working properly)
    df_valid_rt <- df[!is.na(df$rt), ]
    
    if (nrow(df_valid_rt) == 0) {
      return(data.frame(Statistic = "No valid RT data"))
    }
    
    # Group by choice
    rt_by_choice <- split(df_valid_rt$rt, df_valid_rt$choice)
    
    # Create summary data frame
    summary_rows <- lapply(names(rt_by_choice), function(choice_name) {
      rts <- rt_by_choice[[choice_name]]
      data.frame(
        Choice = choice_name,
        Mean = round(mean(rts), 3),
        Median = round(median(rts), 3),
        SD = round(sd(rts), 3),
        Min = round(min(rts), 3),
        Max = round(max(rts), 3),
        N = length(rts)
      )
    })
    
    # Add All row - for all valid RTs
    all_row <- data.frame(
      Choice = "All",
      Mean = round(mean(df_valid_rt$rt), 3),
      Median = round(median(df_valid_rt$rt), 3),
      SD = round(sd(df_valid_rt$rt), 3),
      Min = round(min(df_valid_rt$rt), 3),
      Max = round(max(df_valid_rt$rt), 3),
      N = nrow(df_valid_rt)
    )
    
    # Combine all rows
    summary_df <- do.call(rbind, c(list(all_row), summary_rows))
    
    # Return formatted table
    summary_df
  })
}

# --- 3. Run the application ---
shinyApp(ui = ui, server = server) 