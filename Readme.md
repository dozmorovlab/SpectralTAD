# SpectralTAD

`SpectralTAD` is a TAD caller that uses a modified form of spectral clustering to quickly identify hierarchical topologically associating domains (TADs). Users simply input a contact matrix and recieve a bed file containing the coordinates of TADs and their corresponding levels.

The two main functions are `SpectralTAD` and `SpectralTAD_Par`. `SpectralTAD` is a function for calling TADs and `SpectralTAD_Par` is the parallelized version of this function. The input data can be an n x n, an n x (n+3) or a sparse 3-column matrix.

## Installation

The latest version of SpectralTAD can be directly downloaded from Github using the following code

```{r, eval = FALSE}
library(devtools)
install_github('cresswellkg/SpectralTAD', build_vignettes = TRUE)
library(SpectralTAD)
```

Alternatively, one can download the package from Bioconductor

```{r, eval = FALSE}
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install("HiCcompare")
library(HiCcompare)
```

## Input

There are three types of input accepted:

1. n x n contact matrices
2. 3-column sparse contact matrices
3. n x (n+3) contact matrices

These formats are explained in depth in the vignette.

# Usage

## Multi-Level TADs

```{r}
#Load contact matrix
data("rao_chr20_25_rep")
#Find TADs
tads = SpectralTAD(rao_chr20_25_rep, chr = "chr20", levels = 2, qual_filter = FALSE)
```

The output is a list where each entry corresponds to a level of the TAD hierarchy

First level sample output:

```{r}
chr     start     end  Level
chr20   50000 1200000     1
chr20 1200000 2450000     1
chr20 2450000 3525000     1
chr20 3525000 4075000     1
```

Second level sample output:

```{r}
chr     start     end  Level
chr20   50000  550000     2
chr20  550000  675000     2
chr20  675000 1200000     2
chr20 1200000 1750000     2
```

## Contributions & Support
Suggestions for new features and bug reports are welcome. Please create a new issue for any of these or contact the author directly: @cresswellkg (cresswellkg@vcu.edu)

## Contributors
Authors: @cresswellkg (cresswellkg@vcu.edu) & @mdozmorov (mikhail.dozmorov@vcuhealth.org)




