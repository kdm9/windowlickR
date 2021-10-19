test_that("vcf_parse_correct", {

    v012 = read.table("data/test.raw", header=T)
    vsamp = v012$IID
    v012 = as.matrix(v012[,-c(1:6)])
    class(v012) = "numeric"
    storage.mode(v012) = "numeric"
    colnames(v012) = rownames(v012) = NULL
    expect_equal(dim(v012), c(100, 1000))

    gt = windowlickr:::bcf_getGT("data/test.vcf.gz")
    gt0 = gt$GT
    colnames(gt0) = rownames(gt0) = NULL

    expect_equal(gt$Samples, vsamp)
    expect_equal(dim(gt0), c(100, 1000))
    expect_equal(gt0, v012)

})

test_that("vcf_parse_correct_window", {

    v012 = read.table("data/test.raw", header=T)
    vsamp = v012$IID
    v012 = as.matrix(v012[,-c(1:6)])
    class(v012) = "numeric"
    storage.mode(v012) = "numeric"
    colnames(v012) = rownames(v012) = NULL
    expect_equal(dim(v012), c(100, 1000))

    gt = windowlickr:::bcf_getGT("data/test.vcf.gz", region="1:100-199")
    gt0 = gt$GT
    colnames(gt0) = rownames(gt0) = NULL

    expect_equal(gt$Samples, vsamp)
    expect_equal(dim(gt0), c(100, 100))
    # these are offset by 1 as there is a snp at '0'th position of ref, and
    # therefore base 100 is the 101th snp.
    expect_equal(gt0, v012[,101:200])

})
