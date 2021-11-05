# ruvIIInb
This package removes unwanted variation from raw sequencing count data using Negative Binomial and/or Zero Inflated Negative Binomial model. Thus far, the package has been applied to remove unwanted variation from single-cell RNA-seq and shotgun metagenomics data. The package takes raw sequencing count as input and thus in principle can be applied to data from any sequencing platforms.

# Installation 
```{r eval=FALSE}
devtools::install_github("limfuxing/ruvIIInb",build_vignettes=FALSE)
```
# Vignette
[RUV-III-NB Vignette for scRNA-seq](https://github.com/limfuxing/ruvIIInb/blob/master/inst/doc/ruvIIInbVignette.html)

