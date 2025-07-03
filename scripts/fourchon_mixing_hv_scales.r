#!/usr/bin/env Rscript

# Load required libraries
library(MixSIAR)
library(dplyr)
library(tidyverse)
library(hypervolume)

# Get input arguments: species code and TEF step
args <- commandArgs(trailingOnly = TRUE)
species_code <- tolower(args[1])  # Ensure species code is lowercase
tef_step <- args[2]  # TEF step value (e.g., "1.0", "2.5")
edge_val <- args[3]
buf_val <- args[4]
run_name <- "fa22_isotope_scales"


# Define output directory (include TEF step)
output_dir <- file.path("/work/jnlab/fourchon_isotope_runs", run_name,"outputs", species_code, paste0("tef_step_", tef_step), buf_val)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
cat("Final output directory:", output_dir, "\n")

# Construct paths for mix, source, and TEF files
mix_file <- file.path("/work/jnlab/fourchon_isotope_runs", run_name,"data/mix_files", 
                      paste0("edge_", edge_val, "_buf_",buf_val, "_", species_code, "_mix.csv"))

source_file <- file.path("/work/jnlab/fourchon_isotope_runs", run_name,"data/source_files", 
                         paste0(species_code, "_mix.csv"))

tef_file <- file.path("/work/jnlab/fourchon_isotope_runs", run_name,"data/tef_files", 
                      paste0("tef_", species_code, "_step_", tef_step, ".csv"))

# Ensure that the paired files exist
if (!file.exists(mix_file)) stop(paste("ERROR: Mix file not found:", mix_file))
if (!file.exists(source_file)) stop(paste("ERROR: Source file not found:", source_file))
if (!file.exists(tef_file)) stop(paste("ERROR: TEF file not found:", tef_file))

# Log paths for debugging
cat("Mix file path:", mix_file, "\n")
cat("Source file path:", source_file, "\n")
cat("TEF file path:", tef_file, "\n")

# Function to safely read CSV files
read_safe <- function(file) {
  tryCatch(read.csv(file), error = function(e) {
    cat("Error reading", file, ":", conditionMessage(e), "\n")
    NULL
  })
}

# Load data
animals <- read_safe(mix_file)
producers <- read_safe(source_file)
tef <- read_safe(tef_file)

# Check if any data is missing
if (is.null(animals) || is.null(producers) || is.null(tef)) {
  stop("Missing or invalid input files.")
}

cat("N total:", nrow(animals), "\n")
cat("N complete for d13C & d34S:", sum(complete.cases(animals[, c("d13C", "d34S")])), "\n")
cat("N complete with edge_l.mangrove:", sum(complete.cases(animals[, c("d13C", "d34S", "edge_l.mangrove")])), "\n")


# Prepare MixSIAR data
mix <- load_mix_data(
  filename = mix_file, 
  iso_names = c("d13C", "d34S"), 
  factors = NULL,
  fac_random = NULL,
  fac_nested = NULL,
  cont_effects = "edge_l.mangrove"
)

sources <- load_source_data(
  filename = source_file, 
  conc_dep = TRUE, 
  data_type = "raw", 
  mix = mix
)

discr <- load_discr_data(
  filename = tef_file, 
  mix = mix
)

######## N-Mixes ############
N_mix <- load_mix_data(
  filename = mix_file, 
  iso_names = c("d13C", "d15N"), 
  factors = NULL,
  fac_random = NULL,
  fac_nested = NULL,
  cont_effects = "edge_l.mangrove"
)

N_sources <- load_source_data(
  filename = source_file, 
  conc_dep = TRUE, 
  data_type = "raw", 
  mix = N_mix
)

N_discr <- load_discr_data(
  filename = tef_file, 
  mix = N_mix
)

##### NULL Mixes #########
null_mix <- load_mix_data(
  filename = mix_file, 
  iso_names = c("d13C", "d34S"), 
  factors = NULL,
  fac_random = NULL,
  fac_nested = NULL,
  cont_effects = NULL
)

null_sources <- load_source_data(
  filename = source_file, 
  conc_dep = TRUE, 
  data_type = "raw", 
  mix = null_mix
)

null_discr <- load_discr_data(
  filename = tef_file, 
  mix = null_mix
)

# Write and run the JAGS model
model_file <- file.path(output_dir, paste0("model_", species_code, "_tef_", tef_step, ".txt"))

# Write the JAGS model with corrected argument usage
write_JAGS_model(
  filename = model_file, 
  resid_err = FALSE,   # Residual error disabled
  process_err = TRUE,  # Process error enabled
  mix = mix, 
  source = sources
)

# Run the JAGS model with corrected argument usage
results <- run_model(
  run = "normal", 
  mix = mix, 
  source = sources, 
  discr = discr, 
  model_filename = model_file,
  alpha.prior = 1,
  resid_err = FALSE,   # Include resid_err
  process_err = TRUE   # Include process_err
)

# Run the JAGS model with corrected argument usage
null_results <- run_model(
  run = "normal", 
  mix = null_mix, 
  source = null_source, 
  discr = null_discr, 
  model_filename = model_file,
  alpha.prior = 1,
  resid_err = FALSE,   # Include resid_err
  process_err = TRUE   # Include process_err
)

saveRDS(null_results, file.path(output_dir, sprintf("null_results_%s_tef_%s_buf_%s.rds", species_code, tef_step, buf_val)))
saveRDS(results, file.path(output_dir, sprintf("results_%s_tef_%s_buf_%s.rds", species_code, tef_step, buf_val)))

# Save the current working directory
old_wd <- getwd()

# Set the working directory to the output directory
setwd(output_dir)
# Generate isospace plot and save to output directory
plot_data_filename <- file.path(output_dir, sprintf("isospace_plot_%s_tef_%s.png", species_code, tef_step))

png(filename = plot_data_filename, width = 800, height = 600)
plot_data(
  mix = N_mix,
  source = N_sources,
  discr = N_discr,
  plot_save_pdf = FALSE,
  plot_save_png = FALSE  # Don't let it save automatically
)
dev.off()

cat("Isospace plot saved to:", plot_data_filename, "\n")


# Output diagnostics and results
output_JAGS(
  results, 
  mix, 
  sources, 
  output_options = list(
    summary_save = TRUE,
    summary_name = sprintf("summary_statistics_%s_tef_%s", species_code, tef_step),
    sup_post = TRUE,
    plot_post_save_png = FALSE,
    plot_post_name = sprintf("posterior_density_%s_tef_%s", species_code, tef_step),
    sup_pairs = TRUE,
    plot_pairs_save_png = FALSE,
    plot_pairs_name = sprintf("pairs_plot_%s_tef_%s", species_code, tef_step),
    sup_xy = FALSE,
    plot_xy_save_png = FALSE,
    plot_xy_name = sprintf("xy_plot_%s_tef_%s", species_code, tef_step),
    gelman = TRUE,
    heidel = FALSE,
    geweke = TRUE,
    diag_save = TRUE,
    diag_name = sprintf("diagnostics_%s_tef_%s", species_code, tef_step),
    indiv_effect = FALSE,
    plot_post_save_pdf = FALSE,
    plot_pairs_save_pdf = FALSE,
    plot_xy_save_pdf = FALSE
  )
)

# Output diagnostics and results
output_JAGS(
  null_results, 
  null_mix, 
  null_sources, 
  output_options = list(
    summary_save = TRUE,
    summary_name = sprintf("null_summary_statistics_%s_tef_%s", species_code, tef_step),
    sup_post = TRUE,
    plot_post_save_png = FALSE,
    plot_post_name = sprintf("null_posterior_density_%s_tef_%s", species_code, tef_step),
    sup_pairs = TRUE,
    plot_pairs_save_png = FALSE,
    plot_pairs_name = sprintf("null_pairs_plot_%s_tef_%s", species_code, tef_step),
    sup_xy = FALSE,
    plot_xy_save_png = FALSE,
    plot_xy_name = sprintf("null_xy_plot_%s_tef_%s", species_code, tef_step),
    gelman = TRUE,
    heidel = FALSE,
    geweke = TRUE,
    diag_save = TRUE,
    diag_name = sprintf("null_diagnostics_%s_tef_%s", species_code, tef_step),
    indiv_effect = FALSE,
    plot_post_save_pdf = FALSE,
    plot_pairs_save_pdf = FALSE,
    plot_xy_save_pdf = FALSE
  )
)

# Restore the original working directory
setwd(old_wd)

# --- Start of new code integration ---
p_array <- results$BUGSoutput$sims.list$p.ind

# Dimensions
n_iter <- dim(p_array)[1]
n_ind <- dim(p_array)[2]
n_src <- dim(p_array)[3]

# Flatten to long format
p_df <- as.data.frame.table(p_array, responseName = "p") %>%
  rename(iter = Var1, ind = Var2, src = Var3) %>%
  mutate(
    ind = as.integer(ind),
    src = sources$source_names[as.integer(src)],
    covariate = mix$data$edge_l.mangrove[ind]  # replace with your actual covariate name
  )

# Extract posterior samples of source proportions
write.csv(p_df, file = file.path(output_dir, sprintf("posterior_pind_continuous_%s_tef_%s.csv", species_code, tef_step)), row.names = FALSE)

ilr_cont_summary <- results$BUGSoutput$summary[1:5, ]
write.csv(as.data.frame(ilr_cont_summary), file = file.path(output_dir, sprintf("ilr_cont_summary_%s_tef_%s.csv", species_code, tef_step)), row.names = FALSE)


# Process the posterior_samples data frame

posterior_samples <- as.data.frame(as.table(results$BUGSoutput$sims.list$p.ind))

posterior_samples <- posterior_samples %>%
  mutate(
    Source = sources$source_names[as.numeric(Var3)],
    species_code = species_code,
    ind = as.numeric(Var2),
    run = as.numeric(Var1)
  ) %>%
  select(species_code, Source, Contr = Freq,ind, run)

# Print a message indicating model completion
cat("Completed model for species:", species_code, "\n")

# Get c_d15N (consumer mean d15N)
if ("d15N" %in% names(animals)) {
  animal_d15N <- animals %>%
    mutate(ind = row_number()) %>%
    select(ind, c_d15N = d15N)
  
} else {
  stop("d15N data is not available in the mix data for trophic level calculation.")
}

# Print column names for debugging
cat("Column names in producers data frame:\n")
print(names(producers))

# Get s_d15N (source mean d15N)
if ("d15N" %in% names(producers)) {
  s_d15N_df <- producers %>%
    group_by(Source) %>%
    summarise(Source = first(Source), s_d15N = mean(d15N, na.rm = TRUE), .groups = "drop")
} else if ("Meand15N" %in% names(producers)) {
  s_d15N_df <- producers %>%
    select(Source = Source, s_d15N = Meand15N)
} else {
  stop("Neither d15N nor Meand15N found in the source data for TL calculation.")
}

write.csv(posterior_samples, file = file.path(output_dir, sprintf("posterior_samples_preTL_%s_tef_%s.csv", species_code, tef_step)))
# Merge data and calculate trophic levels
df_mm <- posterior_samples %>%
  left_join(s_d15N_df, by = "Source") %>%
  left_join(animal_d15N, by = "ind") %>%  # now c_d15N is individual-specific
  group_by(species_code, ind, run) %>%
  mutate(
    TL = (c_d15N - sum(Contr * s_d15N)) / 3.4 + 1
  ) %>%
  ungroup()

# Remove duplicates before pivot_wider
df_mm <- df_mm %>%
  distinct()

# Calculate summary statistics for source contributions
summary_contr <- posterior_samples %>%
  group_by(species_code, Source) %>%
  summarize(
    mean_Contr = mean(Contr),
    sd_Contr = sd(Contr),
    p5_Contr = quantile(Contr, 0.05),
    p95_Contr = quantile(Contr, 0.95),
    .groups = 'drop'
  )

# Calculate summary statistics for trophic levels
individual_summary <- df_mm %>%
  group_by(ind) %>%
  summarise(
    mean_tl = mean(TL),
    sd_tl = sd(TL),
    .groups = "drop"
  ) %>% mutate(species_code = species_code)

summary_TL <- df_mm %>%
  group_by(species_code) %>%
  summarize(
    mean_TL = mean(TL),
    sd_TL = sd(TL),
    p5_TL = quantile(TL, 0.05),
    p95_TL = quantile(TL, 0.95),
    .groups = 'drop'
  )


# Write the summary results to CSV files in the output directory
write.csv(summary_contr, file = file.path(output_dir, sprintf("source_contributions_summary_%s_tef_%s.csv", species_code, tef_step)), row.names = FALSE)
write.csv(summary_TL, file = file.path(output_dir, sprintf("trophic_levels_summary_%s_tef_%s.csv", species_code, tef_step)), row.names = FALSE)
write.csv(df_mm, file = file.path(output_dir, sprintf("mixing_model_df_%s_tef_%s.csv", species_code, tef_step)), row.names = FALSE)

