# Data

The five HEST-1k slides analysed in the paper are third-party data and are not
redistributed here. Run `scripts/download_hest.R` to download them from
HuggingFace into `data/HEST/`, about 330 MB in total.

| Slide   | Platform         | Tissue                            |
|---------|------------------|-----------------------------------|
| TENX117 | Xenium           | cutaneous melanoma                |
| TENX115 | Xenium           | cutaneous melanoma                |
| NCBI600 | Visium           | muscle-invasive bladder carcinoma |
| MEND40  | Visium           | cutaneous squamous cell carcinoma |
| NCBI776 | Visium CytAssist | breast invasive ductal carcinoma  |

## Xenograft

The xenograft slides are not part of HEST-1k and are not redistributed here
either. They are the Visium deposit of Zhang, Heiland and Thorsen,
<https://doi.org/10.5281/zenodo.22924632>. Unpack the archive into
`data/zenodo/`, one folder per sample, which is where `XENO_DIR` already points.

Reads of each sample were aligned twice, against GRCh38 and against mm10, so a
unit carries both genomes and a gene belongs to the reference it came from.
`scripts/xenograft/xeno_data_build.R` stacks the two matrices and writes the
objects the rest of that section reads.

Two derived tables are included here, because neither can be rebuilt from that
deposit:

| File               | Rows   | What it holds                                                |
|--------------------|--------|--------------------------------------------------------------|
| `xeno_genes.csv`   | 19,693 | the analysed gene set and the genome each gene mapped to, 8,472 human and 11,221 mouse |
| `tumor_spots.csv`  |    158 | the units annotated as tumour on the four slides that received cells, 47 of them on the reference slide |

The gene assignment was made on data that also carried the human alignment of
slide 1, which the public deposit does not hold; the tumour annotation comes
from a Space Ranger clustering annotated by hand, which is not part of it either.
The consequence for what can and cannot be recomputed is set out in the main
README.
