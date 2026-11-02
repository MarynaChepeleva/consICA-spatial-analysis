
# RUN CONFIG — one row per RUN (slide x K x perturbation fraction).
# Edit here, then source this file to (re)write config/runs.csv.
# bench_ica_perturbations.R and aggregate_results.R read config/runs.csv.
#
# Columns:
#   run_id    spells out the distinguishing params, like the output filenames
#   min_genes NA = auto spot-QC (slides with real background, e.g. Xenium);
#             1  = keep-all spots (full-tissue slides; no over-cut)
#   ncores    consICA cores on perturbation modes (pure_stability is always 1)
#   store_consica_objects  TRUE  = save every iteration's consica_object (getGO flexibility);
#                          FALSE = save only the GO representative (pure_stability iter 1) —
#                                  for big slides that would run out of memory (e.g. NCBI776)
# NB: actual n_genes/n_spots are logged per run in metadata AFTER filtering — not stored here.

df_params <- data.frame(
  run_id = character(), slide = character(), file = character(),
  platform = character(), tissue = character(),
  k = integer(), n_iter = integer(), gene_frac = numeric(), spot_frac = numeric(),
  min_genes = integer(), ncores = integer(), store_consica_objects = logical(),
  stringsAsFactors = FALSE)

# --- Xenium melanoma (targeted panels: auto spot-QC) -------------------------------
df_params <- rbind(df_params, data.frame(run_id = "TENX115__k20__gp05", slide = "TENX115",
  file = "data/HEST/skin/st/TENX115.h5ad", platform = "Xenium", tissue = "melanoma",
  k = 20L, n_iter = 10L, gene_frac = 0.05, spot_frac = 0.05,
  min_genes = NA_integer_, ncores = 4L, store_consica_objects = TRUE))

df_params <- rbind(df_params, data.frame(run_id = "TENX117__k20__gp05", slide = "TENX117",
  file = "data/HEST/skin/st/TENX117.h5ad", platform = "Xenium", tissue = "melanoma",
  k = 20L, n_iter = 10L, gene_frac = 0.05, spot_frac = 0.05,
  min_genes = NA_integer_, ncores = 4L, store_consica_objects = TRUE))

# --- Visium cutaneous SCC, Ji et al. 2020 (full-tissue: keep-all spots) ------------
df_params <- rbind(df_params, data.frame(run_id = "MEND40__k20__gp05", slide = "MEND40",
  file = "data/HEST/skin/st/MEND40.h5ad", platform = "Visium", tissue = "cSCC",
  k = 20L, n_iter = 10L, gene_frac = 0.05, spot_frac = 0.05,
  min_genes = 1L, ncores = 4L, store_consica_objects = TRUE))

df_params <- rbind(df_params, data.frame(run_id = "MEND40__k50__gp05", slide = "MEND40",
  file = "data/HEST/skin/st/MEND40.h5ad", platform = "Visium", tissue = "cSCC",
  k = 50L, n_iter = 10L, gene_frac = 0.05, spot_frac = 0.05,
  min_genes = 1L, ncores = 4L, store_consica_objects = TRUE))

df_params <- rbind(df_params, data.frame(run_id = "MEND40__k20__gp15", slide = "MEND40",
  file = "data/HEST/skin/st/MEND40.h5ad", platform = "Visium", tissue = "cSCC",
  k = 20L, n_iter = 10L, gene_frac = 0.15, spot_frac = 0.05,
  min_genes = 1L, ncores = 4L, store_consica_objects = TRUE))

# MEND38 ran but is not reported; download_hest.R does not fetch it.
# df_params <- rbind(df_params, data.frame(run_id = "MEND38__k50__gp05", slide = "MEND38",
#   file = "data/HEST/skin/st/MEND38.h5ad", platform = "Visium", tissue = "cSCC",
#   k = 50L, n_iter = 10L, gene_frac = 0.05, spot_frac = 0.05,
#   min_genes = 1L, ncores = 1L, store_consica_objects = TRUE))

# --- Visium breast (big -> keep only GO representative, avoids OOM) -----------------
df_params <- rbind(df_params, data.frame(run_id = "NCBI776__k20__gp05", slide = "NCBI776",
  file = "data/HEST/breast/st/NCBI776.h5ad", platform = "Visium", tissue = "breast",
  k = 20L, n_iter = 10L, gene_frac = 0.05, spot_frac = 0.05,
  min_genes = 1L, ncores = 1L, store_consica_objects = FALSE))

df_params <- rbind(df_params, data.frame(run_id = "NCBI776__k50__gp05", slide = "NCBI776",
  file = "data/HEST/breast/st/NCBI776.h5ad", platform = "Visium", tissue = "breast",
  k = 50L, n_iter = 10L, gene_frac = 0.05, spot_frac = 0.05,
  min_genes = 1L, ncores = 1L, store_consica_objects = FALSE))

# --- Visium bladder (low-power) ----------------------------------------------------
df_params <- rbind(df_params, data.frame(run_id = "NCBI600__k20__gp05", slide = "NCBI600",
  file = "data/HEST/bladder/st/NCBI600.h5ad", platform = "Visium", tissue = "bladder",
  k = 20L, n_iter = 10L, gene_frac = 0.05, spot_frac = 0.05,
  min_genes = 1L, ncores = 4L, store_consica_objects = TRUE))

df_params <- rbind(df_params, data.frame(run_id = "NCBI600__k20__gp15", slide = "NCBI600",
  file = "data/HEST/bladder/st/NCBI600.h5ad", platform = "Visium", tissue = "bladder",
  k = 20L, n_iter = 10L, gene_frac = 0.15, spot_frac = 0.05,
  min_genes = 1L, ncores = 4L, store_consica_objects = TRUE))

# --- write -------------------------------------------------------------------------
dir.create("config", showWarnings = FALSE)
# HEST_ROOT overrides the data folder.
HEST_DIR <- Sys.getenv("HEST_ROOT", "data/HEST")
df_params$file <- sub("^data/HEST", HEST_DIR, df_params$file)
write.csv(df_params, "config/runs.csv", row.names = FALSE)
cat("Wrote config/runs.csv with", nrow(df_params), "runs:\n")
print(df_params$run_id)
