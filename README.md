# consICA for spatial transcriptomics <img src="logo_consICA.png" align="right" width="100" alt="consICA" />

Code for the manuscript *Consensus independent component analysis resolves
tumour spatial transcriptomics into interpretable transcriptional programmes*.
The citation and the link to the published version will be added here once it
appears.

consICA resolves a spatial expression matrix into transcriptional programmes
without being given the coordinates of the units. Each programme is a gene
signature with a weight over every unit, and its spatial organization emerges in
that map. On the tumour slides analysed here the components recover
cell-associated signals, melanocytic, immune, vascular, stromal and keratinocyte,
together with technical ones.

The analysis covers five HEST-1k slides, four tumour types on two platforms, and
a human melanoma growing in mouse brain. In the xenograft every gene carries the
genome it mapped to, a label the decomposition never sees, so the separation it
finds can be checked against the species the genes came from.

## Layout

```
scripts/
  download_hest.R  fetches the five slides, about 330 MB
  demo/            the demonstration: one slide end to end
  assembly/        the figures and tables of the manuscript
  benchmark/       the method comparison and its control arm
  xenograft/       the cross-species section
data/              download_hest.R writes here; two derived CSVs are included
output/            written here by the scripts; empty in the repository
sessionInfo.txt    the package versions the analysis ran under
```

Scripts use paths relative to the project root. Open the project, or set the
working directory there, before sourcing anything.

## Getting started

1. Install the packages. `sessionInfo.txt` records the versions the analysis ran
   under.

   ```r
   install.packages(c("Seurat", "Matrix", "ggplot2", "patchwork", "dplyr", "purrr",
                      "ica", "fastICA", "clue", "e1071", "matrixStats", "mgcv",
                      "spdep", "tiff", "jsonlite", "openxlsx", "curl", "ragg",
                      "remotes"))

   if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
   BiocManager::install(c("consICA", "SingleCellExperiment", "zellkonverter", "topGO",
                          "GOSemSim", "AnnotationDbi", "org.Hs.eg.db", "org.Mm.eg.db", "GO.db"))

   # CFS (CellsFromSpace), one of the configurations compared. Pinned
   # to the commit it was run on: the code calls an internal function of it.
   remotes::install_github("gustaveroussy/CFS@3926422e7b5264e92a75d665aa9ef447be88c273")
   ```
2. Accept the dataset terms at <https://huggingface.co/datasets/MahmoodLab/hest>
   and create a read token under Settings, Access Tokens. Paste it into
   `scripts/download_hest.R`.
3. `scripts/download_hest.R` downloads the slides. It pins the dataset to one
   revision, so the files are those the analysis used. The slides go to
   `data/HEST`. To read them from somewhere else, set `HEST_ROOT`.
4. `scripts/demo/consICA_Xenium.R` runs the demonstration on a Xenium slide,
   `scripts/demo/consICA_Visium.R` on a Visium one.
5. The xenograft runs on a separate deposit, see `data/README.md`.

## Description of scripts

**The demonstration** (`scripts/demo/`) runs one slide end to end: decomposed,
annotated, its components tested for spatial autocorrelation and mapped. Minutes
on the Xenium slide, longer on the transcriptome-wide Visium ones. It produces
sections 3.1 and 3.2 of the paper, where the vasculature separates into blood
endothelium, mural cells and lymphatic endothelium, and the keratinocytes into
differentiation stages.

**The comparison** (`scripts/benchmark/` is the record of how its numbers were
produced. `0_run_benchmark.R` holds the sixteen steps in order; its
header gives the file dependency that fixes the order. `DRY_RUN` is `TRUE` by
default, so the script prints the plan and runs nothing until it is set to
`FALSE`; `ONLY` restricts it to part of the list.

**The assembly scripts** (`scripts/assembly/`) build the figures and tables from the
outputs of the demonstration, the comparison and the xenograft.

**The xenograft** (`scripts/xenograft/`) works on reads aligned twice, against
the human and the mouse reference and stacked into one matrix. Those slides are
deposited separately from HEST-1k; `data/README.md` says where they come from and
where to unpack them. `xeno_data_build.R` builds the objects the rest of the
section reads, and the scripts after it run in the order their headers give.

## Data not redistributed

### The xenograft control slide

The public deposit carries only the mouse alignment of slide 1, the tumour-free
control. The cross-species decomposition cannot be run on it, and the scripts
leave it out of the xenograft tables and figure. `HAVE_SLIDE1_HM` is `FALSE` for
this reason; set it to `TRUE` if you obtain the human alignment from the authors
of that dataset.

### The H&E panels

The H&E panels of Figures 1 and 2 are not redistributed and come out as
placeholders. The whole-slide images are in HEST-1k, `wsis = "tif"` in
`download_hest.R`. Zoomed crops were obtained with QuPath.

## Citation

The manuscript is under review; its citation will be added here once it appears.
The code is archived at Zenodo, DOI zenodo.XXXXX.

consICA is published separately:

> Chepeleva M, Kaoma T, Zinovyev A, Toth R, Nazarov PV. consICA: an R package for
> robust reference-free deconvolution of multi-omics data. *Bioinformatics Advances*
> 2024;4(1):vbae102. doi:10.1093/bioadv/vbae102

## Licence

The code is released under the MIT licence, see `LICENSE`. One file is not:
`scripts/benchmark/cfs_ICASpatial.R` is copied from CellsFromSpace and carries
that package's GPL, whose text is in `LICENSE.GPL-2`. The data it reads are
third-party and carry their own terms: the HEST-1k expression matrices are
CC BY-NC-SA 4.0 and are not redistributed here, and `data/README.md` lists the
slides and the studies that contributed them.
