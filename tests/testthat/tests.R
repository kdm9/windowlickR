context("vcf parse")

library(tidyverse)
library(testthat)
test_that("vcf_parse_correct", {

    v012 = read.table("data/test.raw", header=T) %>%
        select(-FID, -PAT, -MAT, -SEX, -PHENOTYPE) %>%
        column_to_rownames("IID") %>%
        as.matrix()
    class(v012) = "numeric"
    storage.mode(v012) = "numeric"
    colnames(v012) = rownames(v012) = NULL

    expect_equal(dim(v012), c(100, 1000))

    gt = windowlickr:::bcf_getGT("data/test.vcf.gz")
    gt0 = gt$GT
    colnames(gt0) = rownames(gt0) = NULL

    expect_equal(dim(gt0), c(100, 1000))


    expect_equal(gt0, v012)


})
