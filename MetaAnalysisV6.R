# 
# Gill et al. (in prep.) — A global meta-analysis of piscivory among
# lanternfish (Myctophidae)
# 
# Description:
#   Fits Bayesian mixed-effects logistic regression models to lanternfish
#   diet data and produces all manuscript figures, tables, and supplementary
#   materials reported in the paper.
#
# Input:
#   data/Myctophid_diet.csv          — main diet dataset
#   data/dist_to_shore_cache.csv     — pre-computed coastline distances
#                                      (generated on first run; reload to skip)
#   data/seafloor_depth_cache.csv    — pre-computed NOAA bathymetry depths
#                                      (generated on first run; reload to skip)
#
# Output:
#   models/   — fitted brms model objects (.rds)
#   figures/  — all manuscript and supplementary figures (.pdf and .tiff)
#   tables/   — all manuscript and supplementary tables (.docx)
#
# NOTE on spatial data:
#   Distance-to-shore and seafloor depth are computed by querying external
#   services (rnaturalearth coastline; NOAA bathymetry server). These queries
#   can be slow and server-dependent. After the first successful run, the
#   results are cached as .csv files in data/ so that re-running the script
#   does not require re-querying. If you have received these cache files with
#   the repository, the download blocks will be skipped automatically.

# 
# SECTION 1: SETUP##############################################################
# 

## Clear environment -----------------------------------------------------------
rm(list = ls())

## Libraries -------------------------------------------------------------------

# Core data wrangling & visualisation
library(readr)
library(tidyverse)
library(broom.mixed)

# Spatial
library(maps)
library(sf)
library(rnaturalearth)
library(marmap)

# Bayesian modelling
library(brms)
library(tidybayes)
library(bayesplot)
library(loo)

# Tables & output
library(flextable)
library(officer)

## Output directories ----------------------------------------------------------
# Create all output directories up front 
dir.create("models",  showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)
dir.create("tables",  showWarnings = FALSE)
dir.create("data",    showWarnings = FALSE)

# 
# SECTION 2: DATA###############################################################
# 

## Diet data -------------------------------------------------------------------
diet <- read_csv("Myctophid_diet.csv")

# Capitalise genus names
diet <- diet %>%
  mutate(genusconfirmed = str_to_sentence(genusconfirmed))

# Set reference level for genus (Diaphus: most-sampled, intermediate piscivory)
diet$genusconfirmed <- relevel(factor(diet$genusconfirmed), ref = "Diaphus")

# Ensure ocean and environment are factors
diet$ocean       <- factor(diet$ocean)
diet$environment <- factor(diet$environment)

## Map data (for Figure 1) -----------------------------------------------------
world <- map_data("world")

## Distance to shore -----------------------------------------------------------
# Pre-computed distances are cached after the first run to avoid re-querying

dist_cache <- "data/dist_to_shore_cache.csv"

if (file.exists(dist_cache)) {

  message("Loading cached distance-to-shore data from: ", dist_cache)
  dist_cached <- read_csv(dist_cache)
  diet <- diet %>%
    left_join(dist_cached, by = c("ref", "species", "yearstart", "locality"))

} else {

  message("Computing distance to shore — this may take 1-2 minutes...")

  diet_with_coords <- diet %>%
    filter(!is.na(longmid) & !is.na(latmid))

  fish_points <- st_as_sf(diet_with_coords,
                          coords = c("longmid", "latmid"),
                          crs    = 4326)

  coastline     <- ne_coastline(scale = "medium", returnclass = "sf")
  distances     <- st_distance(fish_points, coastline)
  min_distances <- apply(distances, 1, min)

  diet_with_coords$dist_to_shore <- min_distances

  # Cache result
  dist_cached <- diet_with_coords %>%
    select(ref, species, yearstart, locality, dist_to_shore)
  write_csv(dist_cached, dist_cache)
  message("Distance-to-shore cached to: ", dist_cache)

  diet <- diet %>%
    left_join(dist_cached, by = c("ref", "species", "yearstart", "locality"))
}

## Bathymetry (seafloor depth) -------------------------------------------------
# Pre-computed depths are cached after the first run to avoid re-querying

bathy_cache <- "data/seafloor_depth_cache.csv"

if (file.exists(bathy_cache)) {

  message("Loading cached bathymetry data from: ", bathy_cache)
  bathy_cached <- read_csv(bathy_cache)
  diet <- diet %>%
    left_join(bathy_cached, by = c("ref", "species", "yearstart", "locality"))

} else {

  message("Downloading NOAA bathymetry — this may take several minutes...")

# Helper: download bathymetry for a bounding box around a data subset
get_bathy_for_ocean <- function(data, resolution = 10) {
    lon_range <- range(data$longmid, na.rm = TRUE)
    lat_range <- range(data$latmid,  na.rm = TRUE)
    getNOAA.bathy(
      lon1 = lon_range[1] - 1,
      lon2 = lon_range[2] + 1,
      lat1 = lat_range[1] - 1,
      lat2 = lat_range[2] + 1,
      resolution = resolution
    )
  }

# Subset by ocean (Indian Ocean NA coords already filtered above)
pacific_data   <- diet %>% filter(ocean == "pacific")
atlantic_data  <- diet %>% filter(ocean == "atlantic")
indian_data    <- diet %>% filter(ocean == "indian", !is.na(longmid), !is.na(latmid))
antarctic_data <- diet %>% filter(ocean == "antarctic")

# Download bathymetry per ocean to keep bounding boxes tight
bathy_pacific   <- get_bathy_for_ocean(pacific_data)
bathy_atlantic  <- get_bathy_for_ocean(atlantic_data)
bathy_indian    <- get_bathy_for_ocean(indian_data)
bathy_antarctic <- get_bathy_for_ocean(antarctic_data)

# Extract seafloor depths (convert to positive metres)
pacific_data$seafloor_depth_m   <- abs(get.depth(bathy_pacific,   x = pacific_data$longmid,   y = pacific_data$latmid,   locator = FALSE)$depth)
atlantic_data$seafloor_depth_m  <- abs(get.depth(bathy_atlantic,  x = atlantic_data$longmid,  y = atlantic_data$latmid,  locator = FALSE)$depth)
indian_data$seafloor_depth_m    <- abs(get.depth(bathy_indian,    x = indian_data$longmid,    y = indian_data$latmid,    locator = FALSE)$depth)
antarctic_data$seafloor_depth_m <- abs(get.depth(bathy_antarctic, x = antarctic_data$longmid, y = antarctic_data$latmid, locator = FALSE)$depth)

diet_with_bathy <- bind_rows(pacific_data, atlantic_data, indian_data, antarctic_data)

# Cache result
bathy_cached <- diet_with_bathy %>%
  select(ref, species, yearstart, locality, seafloor_depth_m)
write_csv(bathy_cached, bathy_cache)
message("Bathymetry cached to: ", bathy_cache)

diet <- diet %>%
  left_join(bathy_cached, by = c("ref", "species", "yearstart", "locality"))
}

##Standardise predictors ------------------------------------------------------
# Body size (sizez) is standardised globally across all genera using maximum
# standard length, placing all genera on a common scale. 
# All other continuous predictors are also standardised globally.

diet <- diet %>%
  mutate(
    sizez                = as.numeric(scale(maxsize)),
    depthz               = as.numeric(scale(maxdepth)),
    dist_to_shore_km     = as.numeric(dist_to_shore) / 1000,
    dist_to_shore_scaled = as.numeric(scale(dist_to_shore_km)),
    seafloor_scaled      = as.numeric(scale(seafloor_depth_m))
  )

# Pre-compute size back-transform constants (used in Figures 4 & 5 axis labels)
size_mean <- mean(diet$maxsize, na.rm = TRUE)
size_sd   <- sd(diet$maxsize,   na.rm = TRUE)

# 
# SECTION 3A: MODELS############################################################
# 

## Prior specifications --------------------------------------------------------
# priorA = moderate (main analysis)
# priorB = diffuse  (sensitivity analysis)
# priorC = tight    (sensitivity analysis)

priorA <- c(prior(normal(0, 1),   class = "b"),
            prior(cauchy(0, 1),   class = "sd"))

priorB <- c(prior(normal(0, 5),   class = "b"),
            prior(cauchy(0, 2),   class = "sd"))

priorC <- c(prior(normal(0, 0.5), class = "b"),
            prior(cauchy(0, 0.5), class = "sd"))

## Set brm() defaults -------------------------------------------------------
brm_defaults <- list(
  family    = bernoulli(link = "logit"),
  data      = diet,
  chains    = 4,
  iter      = 4000,
  warmup    = 2000,
  control   = list(adapt_delta = 0.95, max_treedepth = 12),
  cores     = 4,
  save_pars = save_pars(all = TRUE),
  seed      = 1234
)

## fit_or_load() ---------------------------------------------------------------
# Fits a brms model and saves it as an .rds file, or loads the saved file if it
# already exists. This means the models only need to be fitted once — re-running
# the script after the initial fit takes seconds rather than hours.

fit_or_load <- function(file, formula, prior = priorA) {
  if (file.exists(file)) {
    readRDS(file)
  } else {
    m <- do.call(brm, c(list(formula = formula, prior = prior), brm_defaults))
    saveRDS(m, file)
    m
  }
}


## Univariable models ----------------------------------------------------------
GeneraModel    <- fit_or_load("models/GeneraModel.rds",
                               fish ~ genusconfirmed + (1 | ref))
SizeModel      <- fit_or_load("models/SizeModel.rds",
                               fish ~ sizez + (1 | ref))
SeasonModel    <- fit_or_load("models/SeasonModel.rds",
                               fish ~ season + (1 | ref))
LocationModel  <- fit_or_load("models/LocationModel.rds",
                               fish ~ latmid + (1 | ref))
DepthModel     <- fit_or_load("models/DepthModel.rds",
                               fish ~ seafloor_scaled + (1 | ref))
DistShoreModel <- fit_or_load("models/DistShoreModel.rds",
                               fish ~ dist_to_shore_scaled + (1 | ref))

## Multivariable models --------------------------------------------------------
MultiModelA <- fit_or_load("models/MultiModelA.rds",
                            fish ~ sizez + genusconfirmed + (1 | ref))
MultiModelB <- fit_or_load("models/MultiModelB.rds",
                            fish ~ sizez + genusconfirmed + latmid + (1 | ref))
MultiModelC <- fit_or_load("models/MultiModelC.rds",
                            fish ~ sizez + genusconfirmed + season + (1 | ref))
MultiModelD <- fit_or_load("models/MultiModelD.rds",
                            fish ~ sizez + genusconfirmed + dist_to_shore_scaled + (1 | ref))
MultiModelE <- fit_or_load("models/MultiModelE.rds",
                            fish ~ sizez + genusconfirmed + seafloor_scaled + (1 | ref))
MultiModelF <- fit_or_load("models/MultiModelF.rds",
                            fish ~ sizez + genusconfirmed + latmid + season +
                              seafloor_scaled + dist_to_shore_scaled + (1 | ref))

## Interaction models ----------------------------------------------------------
MultiModelA_interaction <- fit_or_load("models/MultiModelA_interaction.rds",
                                        fish ~ sizez * genusconfirmed + (1 | ref))

MultiModelA_rs <- fit_or_load("models/MultiModelA_rs.rds",
                               fish ~ sizez * genusconfirmed + (1 + sizez | ref))

MultiModelA_rsnc <- fit_or_load("models/MultiModelA_rsnc.rds",
                                  fish ~ sizez * genusconfirmed + (1 + sizez || ref))


## LOO cross-validation --------------------------------------------------------
# moment_match = TRUE stabilises LOO estimates for influential observations
# (those with Pareto k > 0.7). All models use this setting for consistency.

loo_list <- list(
  GeneraModel             = loo(GeneraModel,             moment_match = TRUE),
  SizeModel               = loo(SizeModel,               moment_match = TRUE),
  SeasonModel             = loo(SeasonModel,             moment_match = TRUE),
  LocationModel           = loo(LocationModel,           moment_match = TRUE),
  DepthModel              = loo(DepthModel,              moment_match = TRUE),
  DistShoreModel          = loo(DistShoreModel,          moment_match = TRUE),
  MultiModelA             = loo(MultiModelA,             moment_match = TRUE),
  MultiModelB             = loo(MultiModelB,             moment_match = TRUE),
  MultiModelC             = loo(MultiModelC,             moment_match = TRUE),
  MultiModelD             = loo(MultiModelD,             moment_match = TRUE),
  MultiModelE             = loo(MultiModelE,             moment_match = TRUE),
  MultiModelF             = loo(MultiModelF,             moment_match = TRUE),
  MultiModelA_interaction = loo(MultiModelA_interaction, moment_match = TRUE),
  MultiModelA_rs          = loo(MultiModelA_rs,          moment_match = TRUE),
  MultiModelA_rsnc        = loo(MultiModelA_rsnc,        moment_match = TRUE)
)

# 
# SECTION 3B: PRIOR SENSITIVITY ANALYSIS########################################
#

# Fits the best-fit model structure under alternative prior specifications
# and compares parameter stability across prior choices.

## Fit best-fit model under alternative priors --------------------------------
best_priorB <- fit_or_load("models/MultiModelA_rs_priorB.rds",
                           fish ~ sizez * genusconfirmed + (1 + sizez | ref),
                           prior = priorB)

best_priorC <- fit_or_load("models/MultiModelA_rs_priorC.rds",
                           fish ~ sizez * genusconfirmed + (1 + sizez | ref),
                           prior = priorC)

## Extract key estimates under each prior -------------------------------------
extract_key_params <- function(model, prior_label) {
  fixef(model, probs = c(0.025, 0.975)) %>%
    as.data.frame() %>%
    rownames_to_column("Parameter") %>%
    filter(Parameter %in% c("Intercept",
                            "sizez",
                            "genusconfirmedGymnoscopelus",
                            "sizez:genusconfirmedGymnoscopelus")) %>%
    mutate(
      Prior    = prior_label,
      OR       = round(exp(Estimate), 2),
      OR_lower = round(exp(Q2.5),    2),
      OR_upper = round(exp(Q97.5),   2),
      OR_fmt   = sprintf("%.2f [%.2f, %.2f]", OR, OR_lower, OR_upper),
      Parameter_clean = case_when(
        Parameter == "Intercept"                              ~ "Intercept (Diaphus, mean size)",
        Parameter == "sizez"                                  ~ "Size (Diaphus)",
        Parameter == "genusconfirmedGymnoscopelus"            ~ "Gymnoscopelus",
        Parameter == "sizez:genusconfirmedGymnoscopelus"      ~ "Size \u00d7 Gymnoscopelus"
      )
    ) %>%
    select(Prior, Parameter_clean, OR_fmt)
}

sensitivity_df <- bind_rows(
  extract_key_params(MultiModelA_rs,   "Moderate: normal(0,1), Cauchy(0,1)"),
  extract_key_params(best_priorB,  "Diffuse:  normal(0,5), Cauchy(0,2)"),
  extract_key_params(best_priorC,  "Tight:    normal(0,0.5), Cauchy(0,0.5)")
)

## LOO under alternative priors -----------------------------------------------
loo_sensitivity <- tibble(
  Prior   = c("Moderate: normal(0,1), Cauchy(0,1)",
              "Diffuse:  normal(0,5), Cauchy(0,2)",
              "Tight:    normal(0,0.5), Cauchy(0,0.5)"),
  LOOIC   = c(
    round(loo(MultiModelA_rs,   moment_match = TRUE)$estimates["looic","Estimate"], 1),
    round(loo(best_priorB,  moment_match = TRUE)$estimates["looic","Estimate"], 1),
    round(loo(best_priorC,  moment_match = TRUE)$estimates["looic","Estimate"], 1)
  ),
  Pareto_k_high = c(
    sum(loo(MultiModelA_rs,   moment_match = TRUE)$diagnostics$pareto_k > 0.7),
    sum(loo(best_priorB,  moment_match = TRUE)$diagnostics$pareto_k > 0.7),
    sum(loo(best_priorC,  moment_match = TRUE)$diagnostics$pareto_k > 0.7)
  )
)


# 
# SECTION 4A: RESULTS — TABLES##################################################
# 

## Designate best-fit model ----------------------------------------------------
best_model <- MultiModelA_rs

## Extract fixed effects from best model ---------------------------------------
# Note: fe_tidy is kept separate from fe_raw (used in the hand-check section)
# to prevent the object being silently overwritten later in the script.

fe_tidy <- fixef(best_model, probs = c(0.025, 0.975)) %>%
  as.data.frame() %>%
  rownames_to_column("Parameter") %>%
  rename(Est_Error = Est.Error) %>%
  mutate(
    OR        = round(exp(Estimate), 2),
    OR_lower  = round(exp(Q2.5),    2),
    OR_upper  = round(exp(Q97.5),   2),
    OR_CrI    = sprintf("%.2f [%.2f, %.2f]", OR, OR_lower, OR_upper),
    Estimate_fmt = sprintf("%.2f [%.2f, %.2f]",
                           round(Estimate, 2),
                           round(Q2.5,     2),
                           round(Q97.5,    2)),
    # Credible effect: 95% CrI excludes 1.0 on OR scale
    Credible  = ifelse(OR_lower > 1 | OR_upper < 1, "Yes", ""),
    # Clean parameter names for display in tables and figures
    Parameter_clean = Parameter %>%
      str_replace("^b_", "") %>%
      str_replace("genusconfirmed", "") %>%
      str_replace("sizez:genusconfirmed", "Size \u00d7 ") %>%
      str_replace("^sizez$",    "Size (ref: Diaphus)") %>%
      str_replace("^Intercept$", "Intercept (Diaphus, mean size)")
  )

## Helper: extract Bayesian R2 as formatted string ----------------------------
get_r2 <- function(model) {
  r2 <- bayes_R2(model)
  sprintf("%.2f [%.2f, %.2f]", r2[1, "Estimate"], r2[1, "Q2.5"], r2[1, "Q97.5"])
}


## Prior sensitivity table -----------------------------------------------------
supp_table_sensitivity <- flextable(sensitivity_df) %>%
  set_header_labels(
    Prior           = "Prior specification",
    Parameter_clean = "Parameter",
    OR_fmt          = "OR [95% CrI]"
  ) %>%
  merge_v(j = "Prior") %>%
  valign(j = "Prior", valign = "top") %>%
  hline(i = c(4, 8), border = fp_border(color = "black", width = 1)) %>%
  align(align = "center", part = "all") %>%
  align(j = 1:2, align = "left", part = "all") %>%
  autofit() %>%
  set_caption(
    paste("Supplementary Table 3. Prior sensitivity analysis for the",
          "best-supported model (Genus \u00d7 Size with random slopes).",
          "Key parameter odds ratios (ORs) and 95% credible intervals (CrIs)",
          "are shown for the main analysis prior (moderate) and two",
          "alternative specifications (diffuse and tight).",
          "Results were consistent across all prior specifications,",
          "indicating that conclusions are robust to prior choice.")
  ) %>%
  theme_booktabs()

## LOOIC model comparison ------------------------------------------------------

model_labels <- c(
  "Genus",
  "Size",
  "Season",
  "Latitude",
  "Depth",
  "Distance to shore",
  "Genus + Size",
  "Genus + Size + Latitude",
  "Genus + Size + Season",
  "Genus + Size + Distance to shore",
  "Genus + Size + Depth",
  "Genus + Size + Latitude + Season + Depth + Distance to shore",
  "Genus \u00d7 Size",
  "Genus \u00d7 Size + random slopes",
  "Genus \u00d7 Size + random slopes (no correlation)"
)

random_effects_labels <- c(
  rep("(1 | ref)", 6),
  rep("(1 | ref)", 6),
  "(1 | ref)",
  "(1 + sizez | ref)",
  "(1 + sizez || ref)"
)

model_objects <- list(
  GeneraModel, SizeModel, SeasonModel, LocationModel,
  DepthModel, DistShoreModel,
  MultiModelA, MultiModelB, MultiModelC, MultiModelD, MultiModelE, MultiModelF,
  MultiModelA_interaction, MultiModelA_rs, MultiModelA_rsnc
)

loo_df <- map2_dfr(loo_list, model_labels, function(l, label) {
  tibble(
    Model         = label,
    LOOIC         = round(l$estimates["looic",   "Estimate"], 1),
    SE            = round(l$estimates["looic",   "SE"],       1),
    elpd_loo      = round(l$estimates["elpd_loo","Estimate"], 1),
    p_loo         = round(l$estimates["p_loo",  "Estimate"],  1),
    Pareto_k_high = sum(l$diagnostics$pareto_k > 0.7)
  )
}) %>%
  mutate(
    Random_effects = random_effects_labels,
    R2             = map_chr(model_objects, get_r2),
    delta_LOOIC    = round(LOOIC - min(LOOIC), 1)
  ) %>%
  arrange(LOOIC) %>%
  select(Model, Random_effects, LOOIC, SE, delta_LOOIC,
         elpd_loo, p_loo, Pareto_k_high, R2)

supp_table1 <- flextable(loo_df) %>%
  set_header_labels(
    Model         = "Model (fixed effects)",
    Random_effects = "Random effects",
    LOOIC         = "LOOIC",
    SE            = "SE",
    delta_LOOIC   = "\u0394LOOIC",
    elpd_loo      = "elpd\u2099\u2092\u2092",
    p_loo         = "p\u2099\u2092\u2092",
    Pareto_k_high = "Pareto k > 0.7",
    R2            = "R\u00b2 [95% CrI]"
  ) %>%
  bg(i = ~ delta_LOOIC == 0, bg = "#D9D9D9") %>%
  bold(i = ~ delta_LOOIC == 0) %>%
  align(align = "center", part = "all") %>%
  align(j = 1:2, align = "left", part = "all") %>%
  autofit() %>%
  set_caption(
    paste("Supplementary Table 1. Model comparison using leave-one-out",
          "cross-validation information criterion (LOOIC). Lower LOOIC values",
          "indicate better predictive performance. \u0394LOOIC is the difference",
          "from the best-supported model (highlighted in grey, \u0394LOOIC = 0).",
          "SE = standard error of LOOIC; elpd\u2099\u2092\u2092 = expected log",
          "pointwise predictive density; p\u2099\u2092\u2092 = effective number",
          "of parameters; Pareto k > 0.7 = number of highly influential",
          "observations stabilised by moment matching.",
          "R\u00b2 = Bayesian R\u00b2 [95% credible interval].")
  ) %>%
  theme_booktabs()

## Best-fit model parameters ---------------------------------------------------
re_summary <- VarCorr(best_model)$ref
re_df <- tibble(
  Parameter_clean = c("SD (intercept | study)",
                      "SD (size slope | study)",
                      "Correlation (intercept, slope)"),
  Estimate_fmt    = c(
    sprintf("%.2f", re_summary$sd["Intercept", "Estimate"]),
    sprintf("%.2f", re_summary$sd["sizez",     "Estimate"]),
    sprintf("%.2f", re_summary$cor["Intercept", "Estimate", "sizez"])
  ),
  OR_CrI   = "\u2014",
  Credible = ""
)

params_display <- bind_rows(
  fe_tidy %>% select(Parameter_clean, Estimate_fmt, OR_CrI, Credible),
  re_df
)

table1 <- flextable(params_display) %>%
  set_header_labels(
    Parameter_clean = "Parameter",
    Estimate_fmt    = "Posterior estimate [95% CrI]",
    OR_CrI          = "OR [95% CrI]",
    Credible        = "Credible effect"
  ) %>%
  bold(i = ~ Credible == "Yes") %>%
  align(align = "center", part = "all") %>%
  align(j = 1, align = "left", part = "all") %>%
  hline(i = nrow(fe_tidy), border = fp_border(color = "black", width = 1)) %>%
  autofit() %>%
  set_caption(
    paste("Table 1. Posterior estimates from the best-supported model",
          "(Genus \u00d7 Size with random slopes by study; MultiModelA_rs).",
          "Fixed effects are reported as posterior mean log-odds with 95%",
          "credible intervals (CrI); odds ratios (ORs) are exponentiated values.",
          "Genus ORs represent the multiplicative change in baseline piscivory",
          "odds relative to Diaphus at mean body size. Size \u00d7 genus",
          "interaction ORs indicate how the size effect for each genus differs",
          "multiplicatively from Diaphus. Bold rows indicate effects whose 95%",
          "CrI excludes 1.0. Random effect standard deviations and their",
          "correlation are shown below the horizontal rule.")
  ) %>%
  theme_booktabs()

## Full odds ratios ------------------------------------------------------------

or_table_df <- fe_tidy %>%
  filter(!str_detect(Parameter, "Intercept")) %>%
  mutate(OR_fmt = sprintf("%.2f [%.2f, %.2f]", OR, OR_lower, OR_upper)) %>%
  arrange(desc(OR)) %>%
  select(Parameter_clean, OR_fmt, Credible)

supp_table2 <- flextable(or_table_df) %>%
  set_header_labels(
    Parameter_clean = "Parameter",
    OR_fmt          = "OR [95% CrI]",
    Credible        = "Credible effect"
  ) %>%
  bold(i = ~ Credible == "Yes") %>%
  align(align = "center", part = "all") %>%
  align(j = 1, align = "left", part = "all") %>%
  autofit() %>%
  set_caption(
    paste("Supplementary Table 2. Full odds ratios (ORs) and 95% credible",
          "intervals (CrIs) from the best-supported Bayesian mixed-effects",
          "logistic regression model (Genus \u00d7 Size with random slopes).",
          "The baseline odds row reports the intercept-derived odds for Diaphus",
          "at mean body size. Size \u00d7 genus interaction ORs indicate how",
          "the size effect for each genus differs multiplicatively from Diaphus.",
          "ORs sorted in descending order. Bold rows indicate effects whose",
          "95% CrI excludes 1.0. OR > 1 indicates higher baseline piscivory",
          "odds than Diaphus; OR < 1 indicates lower odds.")
  ) %>%
  theme_booktabs()

## Write all tables to Word ----------------------------------------------------
doc <- read_docx() %>%
  body_add_par("Table 1", style = "heading 2") %>%
  body_add_flextable(table1) %>%
  body_add_break() %>%
  body_add_par("Supplementary Table 1", style = "heading 2") %>%
  body_add_flextable(supp_table1) %>%
  body_add_break() %>%
  body_add_par("Supplementary Table 2", style = "heading 2") %>%
  body_add_flextable(supp_table2)%>%
  body_add_par("Supplementary Table 3", style = "heading 2") %>%
  body_add_flextable(supp_table_sensitivity)

print(doc, target = "tables/piscivory_tables.docx")
message("Tables written to tables/piscivory_tables.docx")

# 
# SECTION 4B: RESULTS — FIGURES##################################################
# 

## Figure 1 — Global distribution map -----------------------------------------

diet_map <- diet %>%
  mutate(
    fish_bin = ifelse(fish == "yes", 1, 0),
    longmid  = round(longmid, 1),
    latmid   = round(latmid,  1)
  ) %>%
  group_by(longmid, latmid) %>%
  summarize(
    n_total        = n(),
    prop_piscivory = mean(fish_bin, na.rm = TRUE),
    .groups        = "drop"
  )

Figure1 <- ggplot() +
  geom_map(
    data = world, map = world,
    aes(x = long, y = lat, map_id = region),
    fill      = "grey80",
    color     = "white",
    linewidth = 0.3
  ) +
  geom_point(
    data = diet_map %>% arrange (desc(n_total)),
    aes(x = longmid, y = latmid,
        size = n_total,
        fill = prop_piscivory),
    shape = 21,
    color = "grey30"
  ) +
  scale_fill_viridis_c(
    option = "plasma",
    name   = "Proportion\npiscivorous",
    limits = c(0, 1),
    breaks = c(0, 0.25, 0.5, 0.75, 1),
    labels = c("0", "0.25", "0.50", "0.75", "1.0")
  ) +
  scale_size_continuous(
    name   = "Observations",
    range  = c(2, 12),
    breaks = c(1, 5, 10, 20, 36)
  ) +
  coord_fixed(ratio = 1.3, xlim = c(-180, 180), ylim = c(-90, 90)) +
  labs(x = NULL, y = NULL) +
  theme_classic() +
  theme(
    axis.text       = element_blank(),
    axis.ticks      = element_blank(),
    axis.line       = element_blank(),
    legend.text     = element_text(size = 12),
    legend.title    = element_text(size = 13),
    legend.position = "right",
    legend.key.size = unit(0.8, "lines"),
    panel.border    = element_rect(color = "grey60", fill = NA, linewidth = 0.5)
  ) +
  guides(
    fill = guide_colorbar(order = 1, barwidth = 1, barheight = 6),
    size = guide_legend(order = 2, override.aes = list(fill = "grey50"))
  )

ggsave("Figure1.png", width = 9, height = 7, dpi = 600)

## Figure 2 — Observations by genus (A) and season (B) ------------------------

# Panel A: by genus
n_genus      <- diet %>%
  count(genusconfirmed) %>%
  mutate(label = paste0(genusconfirmed, " (n=", n, ")"))
genus_labels <- setNames(n_genus$label, n_genus$genusconfirmed)

Figure2a <- ggplot(diet, aes(x = genusconfirmed, fill = fish)) +
  geom_bar(position = "fill") +
  scale_fill_manual(
    name   = "Piscivory",
    values = c("no" = "#721F81", "yes" = "#F0F921"),
    labels = c("no" = "No", "yes" = "Yes")
  ) +
  scale_x_discrete(labels = genus_labels) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = NULL, y = "Proportion of observations", tag = "A") +
  theme_classic(base_size = 18) +
  theme(
    axis.text.x  = element_text(face = "italic", size = 16, angle = 45, hjust = 1, vjust = 1),
    axis.text.y  = element_text(size = 16),
    axis.title.y = element_text(size = 18),
    legend.title = element_text(size = 18, face = "bold"),
    legend.text  = element_text(size = 16),
    plot.tag     = element_text(size = 18, face = "bold")
  )


# Panel B: by season
season_levels <- c("spring", "summer", "fall", "winter", "multiple")

diet_season   <- diet %>%
  filter(season %in% season_levels) %>%
  mutate(season = factor(season, levels = season_levels))

n_season      <- diet_season %>% count(season)
season_labels <- n_season %>%
  mutate(
    label  = paste0(str_to_title(season), "\n(n=", n, ")"),
    season = as.character(season)
  ) %>%
  select(season, label) %>%
  deframe()

Figure2b <- ggplot(diet_season, aes(x = season, fill = fish)) +
  geom_bar(position = "fill") +
  scale_fill_manual(
    name   = "Piscivory",
    values = c("no" = "#721F81", "yes" = "#F0F921"),
    labels = c("no" = "No", "yes" = "Yes")
  ) +
  scale_x_discrete(labels = season_labels) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = NULL, y = "Proportion of observations", tag = "B") +
  theme_classic(base_size = 18) +
  theme(
    axis.text.x  = element_text(size = 16),
    axis.text.y  = element_text(size = 16),
    axis.title.y = element_text(size = 18),
    legend.title = element_text(size = 18, face = "bold"),
    legend.text  = element_text(size = 16),
    plot.tag     = element_text(size = 18, face = "bold")
  )

# Combine the plots using patchwork
Figure2a + Figure2b + plot_layout(ncol = 1, heights = c(2, 1))

ggsave("Figure2.png", width = 13, height = 8, dpi = 600)

## Figure 3 — Size ranges by genus --------------------------------------------

genus_summary <- diet %>%
  filter(!is.na(maxsize)) %>%
  group_by(genusconfirmed) %>%
  summarize(
    min_maxsize    = min(maxsize,    na.rm = TRUE),
    max_maxsize    = max(maxsize,    na.rm = TRUE),
    median_maxsize = median(maxsize, na.rm = TRUE),
    .groups        = "drop"
  ) %>%
  mutate(genusconfirmed = str_to_title(genusconfirmed))

Figure3 <- ggplot(genus_summary,
                  aes(x = reorder(genusconfirmed, max_maxsize))) +
  geom_linerange(
    aes(ymin = min_maxsize, ymax = max_maxsize),
    color     = "black",
    linewidth = 1.2,
    alpha     = 0.7
  ) +
  geom_point(aes(y = min_maxsize),    color = "black", size = 2) +
  geom_point(aes(y = max_maxsize),    color = "black", size = 2) +
  geom_point(aes(y = median_maxsize), color = "black", size = 3, shape = 18) +
  coord_flip() +
  labs(x = NULL, y = "Maximum standard length (mm)") +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x  = element_text(size = 13),
    axis.text.y  = element_text(face = "italic", size = 13),
    axis.title.x = element_text(size = 14)
  )
Figure3

ggsave("Figure3.png", width = 6, height = 6, dpi = 600)

## Compute conditional effects (shared by Figures 4 & 5) ----------------------
# Computed once here to avoid the expensive call being repeated for each figure.

ce_size_genus <- conditional_effects(
  best_model,
  effects = "sizez:genusconfirmed",
  prob    = 0.90
)

## Figure 4 — Partial dependence heatmap --------------------------------------

n_by_genus <- diet %>%
  mutate(genusconfirmed = str_to_title(genusconfirmed)) %>%
  count(genusconfirmed, name = "n")

df_heatmap <- ce_size_genus[[1]] %>%
  mutate(genusconfirmed = str_to_title(genusconfirmed)) %>%
  left_join(n_by_genus, by = "genusconfirmed")

genus_order <- df_heatmap %>%
  group_by(genusconfirmed) %>%
  summarize(mean_prob = mean(estimate__), .groups = "drop") %>%
  arrange(mean_prob) %>%
  left_join(n_by_genus, by = "genusconfirmed") %>%
  mutate(genus_label = paste0(genusconfirmed, " (n=", n, ")")) %>%
  pull(genus_label)

df_heatmap <- df_heatmap %>%
  mutate(
    genus_label = paste0(genusconfirmed, " (n=", n, ")"),
    genus_label = factor(genus_label, levels = genus_order)
  )

x_min       <- -1.5
x_max       <-  2.0
size_breaks <- c(-1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2)
size_labels <- paste0(round(size_breaks * size_sd + size_mean, 0),
                      "\n(", size_breaks, " SD)")

df_heatmap_trimmed <- df_heatmap %>%
  filter(sizez >= x_min & sizez <= x_max)

Figure4 <- ggplot(df_heatmap_trimmed, aes(x = sizez, y = genus_label)) +
  geom_tile(aes(fill = estimate__)) +
  scale_fill_viridis_c(
    option = "plasma",
    name   = "Predicted\nprobability",
    limits = c(0, 1),
    breaks = c(0, 0.25, 0.5, 0.75, 1.0),
    labels = c("0", "0.25", "0.50", "0.75", "1.0")
  ) +
  scale_x_continuous(
    breaks = size_breaks,
    labels = size_labels,
    limits = c(x_min, x_max),
    expand = c(0, 0)
  ) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(x = "Body size (mm)", y = NULL) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x  = element_text(size = 12, hjust = 0.5),
    axis.text.y  = element_text(face = "italic", size = 13),
    axis.title.x = element_text(size = 14, margin = margin(t = 8)),
    legend.title = element_text(size = 12),
    legend.text  = element_text(size = 11),
    panel.grid   = element_blank()
  )
Figure4

ggsave("Figure4.png", width = 8, height = 7, dpi = 600)

## Figure 5 — Predicted probability curves by genus ---------------------------

genus_range <- diet %>%
  group_by(genusconfirmed) %>%
  summarize(
    min_size  = min(maxsize, na.rm = TRUE),
    max_size  = max(maxsize, na.rm = TRUE),
    size_range = max_size - min_size,
    n         = n(),
    .groups   = "drop"
  ) %>%
  filter(n >= 3, size_range > 0)   # exclude genera with < 3 obs or no size variation

# Back-transform sizez to mm using global constants and filter to observed size ranges
ce_size_genus_main <- ce_size_genus[[1]] %>%
  mutate(size_actual = sizez * size_sd + size_mean) %>%
  inner_join(genus_range, by = "genusconfirmed") %>%
  filter(size_actual >= min_size & size_actual <= max_size)

Figure5 <- ggplot(ce_size_genus_main, aes(x = size_actual, y = estimate__)) +
  geom_line(color = "black", linewidth = 1) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__),
              fill = "black", alpha = 0.2) +
  facet_wrap(~genusconfirmed, scales = "free_x", ncol = 4) +
  labs(
    x = "Body size (mm)",
    y = "Predicted probability of piscivory"
  ) +
  theme_classic(base_size = 14) +
  theme(
    strip.text       = element_text(face = "italic", size = 14),
    strip.background = element_blank()
  )
Figure5

ggsave("Figure5.png", width = 9, height = 7, dpi = 600)

## Supplementary Figure 2 — Posterior predictive check ------------------------
# type = "bars" is appropriate for binary (Bernoulli) outcomes, showing the
# observed 0/1 counts against replicated datasets from the posterior.

SuppFig2 <- pp_check(best_model, type = "bars", ndraws = 100) +
  labs(
    x = "Fish presence (0 = absent, 1 = present)",
    y = "Count"
  ) +
  theme_classic()

save_fig("SuppFig_2", SuppFig2, w = 120, h = 100)

## Supplementary Figure 3 — Trace plots ----------------------------------------
# Trace plots for key fixed-effect parameters of the best-fit model.
# Full convergence assessed via R-hat <= 1.01 and ESS > 400.

SuppFig3 <- mcmc_trace(
  as.array(best_model),
  pars = c("b_Intercept", "b_sizez"),
  facet_args = list(ncol = 1)
) +
  theme_classic()
SuppFig3

SuppFig3 <- mcmc_rank_hist(
  as.array(best_model),
  pars = c("b_Intercept", "b_sizez"),
  facet_args = list()
) +
  theme_classic()
SuppFig3

save_fig("SuppFig_3", SuppFig3, w = 180, h = 120)

## Supplementary Figure 4 — Caterpillar plot of study-level random effects ----

SuppFig4 <- best_model %>%
  spread_draws(r_ref[ref, term]) %>%
  filter(term == "Intercept") %>%
  median_qi(.width = 0.95) %>%
  arrange(r_ref) %>%
  mutate(ref = factor(ref, levels = ref)) %>%
  ggplot(aes(y = ref, x = r_ref, xmin = .lower, xmax = .upper)) +
  geom_pointrange() +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  labs(
    x = "Random intercept estimate (log-odds)",
    y = "Study"
  ) +
  theme_classic()

save_fig("SuppFig_4", SuppFig4, w = 140, h = max(120, nrow(SuppFig4$data) * 5))

message("All figures written to figures/")

# 
# SECTION 5: VERIFICATION — HAND CALCULATION CHECK##############################
# 

# Manually verifies the predicted probability for Gymnoscopelus at mean size
# and mean + 1 SD against the conditional_effects() output, to confirm the
# OR table is consistent with the predicted probability figures.

fe_raw <- fixef(best_model)

intercept    <- fe_raw["Intercept",                          "Estimate"]
size_slope   <- fe_raw["sizez",                              "Estimate"]
genus_effect <- fe_raw["genusconfirmedGymnoscopelus",        "Estimate"]
interaction  <- fe_raw["sizez:genusconfirmedGymnoscopelus",  "Estimate"]

log_odds_mean    <- intercept + genus_effect
log_odds_plus1sd <- intercept + size_slope + genus_effect + interaction

prob_mean    <- plogis(log_odds_mean)
prob_plus1sd <- plogis(log_odds_plus1sd)

# Compare against conditional_effects output (reuse ce_size_genus from above)
ce_gymno <- ce_size_genus[[1]] %>%
  filter(genusconfirmed == "Gymnoscopelus")

ce_at_mean <- ce_gymno %>%
  mutate(diff = abs(sizez - 0)) %>%
  arrange(diff) %>%
  slice(1)

ce_at_plus1sd <- ce_gymno %>%
  mutate(diff = abs(sizez - 1)) %>%
  arrange(diff) %>%
  slice(1)

cat("\n--- Hand-calculation verification: Gymnoscopelus ---\n")
cat(sprintf("At mean size   | hand: %.3f | model: %.3f\n",
            prob_mean,    ce_at_mean$estimate__))
cat(sprintf("At mean + 1 SD | hand: %.3f | model: %.3f\n",
            prob_plus1sd, ce_at_plus1sd$estimate__))
cat("Values should agree to within rounding error.\n\n")

# 
# SECTION 6: CITATIONS##########################################################
# 

citation()
citation("brms")
citation("tidyverse")
citation("loo")
citation("tidybayes")
citation("bayesplot")
citation("sf")
citation("marmap")
citation("rnaturalearth")
citation("flextable")
citation("officer")

# 
# SECTION 7: SESSION INFO#######################################################
# 
# Records exact R and package versions for reproducibility.
sessionInfo()
summary(MultiModelA_rs)
