# Widespread piscivory in lanternfishes (Myctophidae) revealed by meta-analysis
 
**Authors:** Alisha M. Gill, Stephen E. Swearer, Stephen R. Wing, Jeffrey S. Shima
**Journal:** Fish and Fisheries | **Year:** 2026
**DOI:** in prep
 
---
 
## Overview
 
This repository contains all data and code required to reproduce the analyses and
figures reported in Gill et al. (in prep.), which quantifies the prevalence and
predictors of piscivory across the family Myctophidae using a global synthesis of
published gut content studies.
 
This research was funded by the Royal Society of New Zealand (Marsden Fund Award
number MFP-VUW-2303).
 
---
 
## Repository contents
 
    Myctophid_diet.csv    — compiled diet dataset (321 records, 38 variables)
    MetaAnalysisV6.R      — R script to reproduce all models, tables, and figures
    READ_ME.txt           — this file
 
---
 
## Data
 
### Myctophid_diet.csv
 
Each row represents one species–study combination (321 records across 31 genera,
130 species, drawn from publications spanning 1970–2024). The dataset was compiled
from gut content literature and records the taxonomic identity, sampling location,
sampling period, body-size range, and dietary composition for each observation.
 
| Column          | Description |
|-----------------|-------------|
| ref             | Short citation key (FirstAuthorYEARletter; e.g. alwis1988a) |
| pubdate         | Year of publication (1970–2024) |
| months          | Calendar month(s) of sampling |
| yearstart       | First year of sampling |
| yearend         | Last year of sampling |
| season          | Season of sampling (spring / summer / fall / winter / multiple) |
| environment     | Broad environmental zone (tropical / temperate / polar / antarctic / arctic) |
| ocean           | Ocean basin (atlantic / pacific / indian / antarctic) |
| locality        | Named sampling locality or region |
| mindepth        | Minimum sampling depth (m) |
| maxdepth        | Maximum sampling depth (m) |
| latmin          | Minimum latitude of sampling area (decimal degrees) |
| latmid          | Midpoint latitude (decimal degrees) |
| latmax          | Maximum latitude (decimal degrees) |
| longmin         | Minimum longitude (decimal degrees) |
| longmid         | Midpoint longitude (decimal degrees) |
| longmax         | Maximum longitude (decimal degrees) |
| genus           | Genus as reported in source publication |
| genusconfirmed  | Accepted genus name (verified against FishBase) |
| species         | Species as reported in source publication |
| spconfirmed     | Accepted species name (verified against FishBase) |
| stage           | Life-history stage sampled (juvenile / mature / both), where reported |
| minsize         | Minimum standard length of sampled individuals (mm) |
| midsize         | Midpoint standard length (mm) |
| maxsize         | Maximum standard length (mm) |
| combmaxsize     | Combined maximum size: reported maxsize if available, otherwise genus-level FishBase maximum (mm) |
| genusfbms       | Genus-level maximum standard length from FishBase (mm) |
| fbmaxsize       | Species-level maximum standard length from FishBase (mm), where available |
| primarydiet     | Most abundant prey item (prey taxon or group) |
| primarygroup    | Broad taxonomic group of primary diet item |
| secondarydiet   | Second most abundant prey item, where reported |
| secondarygroup  | Broad taxonomic group of secondary diet item |
| tertiarydiet    | Third most abundant prey item, where reported |
| tertiarygroup   | Broad taxonomic group of tertiary diet item |
| copepods        | Whether copepods were recorded in the diet (yes / no) |
| fish            | Whether fish were recorded in the diet (yes / no) — primary response variable |
| fishid          | Identity of fish prey, where recorded to species or genus level |
| notes           | Miscellaneous notes |
 
---
 
## Code
 
### MetaAnalysisV6.R
 
A fully annotated R script that reproduces all analyses, tables, and figures in
the manuscript. The script is organised into seven sections:
 
1. **Setup** — loads required packages and creates output directories
   (`models/`, `figures/`, `tables/`, `data/`)
2. **Data** — imports `Myctophid_diet.csv`, computes distance-to-shore and
   seafloor depth for each record (queried from rnaturalearth and the NOAA
   bathymetry server; results are cached as .csv files in `data/` after the
   first run to avoid repeated downloads)
3. **Models** — fits Bayesian mixed-effects logistic regression models using
   `brms`, including univariable, multivariable, and interaction models, with
   LOO cross-validation for model comparison (Section 3A), and a prior
   sensitivity analysis (Section 3B)
4. **Results** — produces all manuscript tables (Section 4A) and figures
   (Section 4B), including the global distribution map (Figure 1), genus- and
   season-level summaries (Figure 2), size-range plots (Figure 3), partial
   dependence heatmaps (Figure 4), predicted probability curves (Figure 5), a
   forest plot of genus-level odds ratios (Figure 6), and supplementary
   diagnostic figures
5. **Verification** — hand-calculation checks on key model estimates
6. **Citations** — R package citations
7. **Session info** — records the R and package versions used
 
#### Key dependencies
 
    tidyverse, brms, tidybayes, bayesplot, loo,
    sf, rnaturalearth, marmap,
    flextable, officer
 
Run the script from the repository root directory. All file paths are relative.
Output files are written to `models/`, `figures/`, and `tables/`.
 
> **Note on spatial data:** distance-to-shore and seafloor depth are retrieved
> from external services on first run. If the cache files
> (`data/dist_to_shore_cache.csv`, `data/seafloor_depth_cache.csv`) are included
> in the repository, these download steps are skipped automatically.
 
---
 
## Contact
 
Alisha M. Gill — alisha.gill@vuw.ac.nz
