# WindowLickR

<!-- badges: start -->
[![R-CMD-check](https://github.com/kdm9/windowlickR/workflows/R-CMD-check/badge.svg)](https://github.com/kdm9/windowlickR/actions)
<!-- badges: end -->

A R/`foreach` interface to genome windows of genetic variation in BCF or VCF files.

**NB: please upgrade to at least version 0.3.0: previous versions produced incorrect output**.

# Usage

```R
library(windowlickr)
library(doMC)
registerDoMC()

# main function: a foreach loop over windows of SNPs across the genome
result =  windowlickr::windowlickr(
    bcf="path/to/samples.bcf",
    windowsize=1e6,
    func=function(window, genos) {
        # window is a list containing info on this window, e.g.:
        # window = list(region="Chr1:1-1000", contig="Chr1", start=1, stop=1000)

        # genos is the result of windowlickr::bcf_getGT, which has genotypes
        # (012-coded, missing is NA) and a bunch of per-snp stats & metadata.
        # genos = list(GT=matrix(...))

        # example code: return a list with a few stats
        list(
            n.snp = genos$nSNP,
            mean.maf = mean(genos$MAF),
            sample.dist = dist(genos$GT)
        )
    })
```

## Changes

- Version 0.2.0: implement VCF support
- Version 0.1.0: initial version


# Credits

[Dr. Kevin D Murray](https://kdmurray.id.au); of Weigel Group, MPI-DB, Tübingen, DE/Borevitz Group, ANU, Canberra, AU

With apology to Richard David James...
