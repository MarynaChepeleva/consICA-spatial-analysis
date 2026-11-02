# Download the five HEST-1k slides analysed in the paper, about 330 MB.

# To download from HEST-1k you need a HuggingFace token: accept the dataset
# terms at https://huggingface.co/datasets/MahmoodLab/hest, then create one
# under Settings, Access Tokens, and paste it below.
# REVISION is the dataset version used here.
# Files already downloaded are skipped.

token <- "" # paste your HuggingFace token here

# Set HEST_ROOT in .Renviron, or put the path in the line below.
# Sys.setenv(HEST_ROOT = "D:/HEST")
DATA_DIR <- Sys.getenv("HEST_ROOT", "data/HEST")
REVISION <- "7e8d5a0b0aace41d8c8ec0f6ecea80e4ad2a61ec"
BASE <- sprintf("https://huggingface.co/datasets/MahmoodLab/hest/resolve/%s", REVISION)

ORGAN <- c(TENX115 = "skin", TENX117 = "skin", MEND40 = "skin",
           NCBI776 = "breast", NCBI600 = "bladder")
SUBDIRS <- c(st = "h5ad", metadata = "json") # add wsis = "tif" (~7 GB) or patches = "h5"

handle <- curl::new_handle()
curl::handle_setheaders(handle, Authorization = paste("Bearer", token))

for (slide in names(ORGAN)) for (sub in names(SUBDIRS)) {
  f <- paste0(slide, ".", SUBDIRS[[sub]])
  dest <- file.path(DATA_DIR, ORGAN[[slide]], sub, f)
  if (file.exists(dest)) next

  dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
  curl::curl_download(file.path(BASE, sub, f), dest, handle = handle)
  message(dest)
}
