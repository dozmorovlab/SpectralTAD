# SpectralTAD

`SpectralTAD` is a TAD caller that uses a modified form of spectral clustering 
to quickly identify hierarchical topologically associating domains (TADs). 
Users input a contact matrix and receive a [BED](https://genome.ucsc.edu/FAQ/FAQformat.html#format1) file 
containing the coordinates of TADs and their levels in a hierarchy.
The Level 1 TADs are generally large, well-defined, while the subsequent levels
are less well-pronounced yet sufficiently distinct to be recognized as TADs.

The two main functions are `SpectralTAD()` and `SpectralTAD_Par()`. 
`SpectralTAD()` is a function for calling TADs. `SpectralTAD_Par()` 
is the parallelized version. The input data can be an $n \times n$, 
an $n \times (n+3)$, or a sparse 3-column matrix (see the [vignette](vignettes/SpectralTAD.Rmd)).

## Installation

The latest version of `SpectralTAD` can be directly installed from Github:

```
devtools::install_github('cresswellkg/SpectralTAD', build_vignettes = TRUE)
library(SpectralTAD)
```

Alternatively, the package can be installed from Bioconductor (to be submitted):

```
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install("SpectralTAD")
library(SpectralTAD)
```

## Input

There are three types of input accepted:

1. $n \times n$ contact matrices
2. $n \times (n+3)$ contact matrices
3. 3-column sparse contact matrices

These formats are explained in depth in the [vignette](vignettes/SpectralTAD.Rmd).

# Usage

## Multi-Level TADs

```
#Load example contact matrix
data("rao_chr20_25_rep")
#Find TADs
tads = SpectralTAD(rao_chr20_25_rep, chr = "chr20", levels = 2, qual_filter = FALSE)
```

The output is a list where each entry corresponds to a level of the TAD hierarchy.

First level sample output:

```
chr     start     end  Level
chr20   50000 1200000     1
chr20 1200000 2450000     1
chr20 2450000 3525000     1
chr20 3525000 4075000     1
```

Second level sample output:

```
chr     start     end  Level
chr20   50000  550000     2
chr20  550000  675000     2
chr20  675000 1200000     2
chr20 1200000 1750000     2
```

## Contributions & Support

Suggestions for new features and bug reports are welcome. Please, create a new 
issue for any of these or contact the author directly: 
@cresswellkg (cresswellkg@vcu.edu)

## Contributors

Authors: @cresswellkg (cresswellkg@vcu.edu) & @mdozmorov (mikhail.dozmorov@vcuhealth.org)
