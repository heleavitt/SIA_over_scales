library(tidyverse)

# 1. Find all ilr_cont_summary files under outputs/  
ilr_files <- list.files(
  path       = "outputs",
  pattern    = "^ilr_cont_summary_.*\\.csv$",
  recursive  = TRUE,
  full.names = TRUE
)

# 2. Read & bind them, extracting metadata from the file path
ilr_df <- map_dfr(ilr_files, function(fp) {
  parts    <- str_split(fp, "/")[[1]]
  species  <- parts[2]
  tef_step <- parts[3] %>% str_remove("^tef_step_")
  scale    <- as.integer(parts[4])
  
  # read the ILR summary (one row per source)
  summ <- read_csv(fp, show_col_types = FALSE)
  
  # read your raw source‐data CSV
  raw_src_fp <- file.path(
    "data", "source_files",
    paste0(species, "_mix.csv")
  )
  raw_src <- read_csv(raw_src_fp, show_col_types = FALSE)
  
  # pull the source‐name column, dedupe, sort alphabetically
  source_names <- raw_src %>% 
    pull(Source) %>% 
    unique() %>% 
    sort()
  
  # sanity check
  if (length(source_names) != nrow(summ)) {
    stop(
      "Got ", length(source_names), 
      " unique sorted source names but ", nrow(summ), 
      " ILR rows for ", species, "/", tef_step
    )
  }
  
  # add the alphabetized names + metadata
  summ %>%
    mutate(
      Source   = source_names,
      species  = species,
      tef_step = tef_step,
      scale    = scale
    ) %>%
    select(Source, everything())
})

glimpse(ilr_df)
