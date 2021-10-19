library(windowlickr)
library(doMC)
registerDoMC()

# main function: a foreach loop over windows of SNPs across the genome
result =  windowlickr::windowlickr(
    bcf="/data/kevin/work/genotyping-data/hpa/HaR.filtered_snps_final.PASS.bi.hardFiltered.indFiltered.noMit.reheader.bcf",
    windowsize=1e6,
    func=function(window, genos) {
        # window is a list containing info on this window, e.g.:
        # window = list(region="Chr1:1-1000", contig="Chr1", start=1, stop=1000)

        # genos is the result of windowlickr::bcf_getGT, which has genotypes
        # (012-coded, missing is NA) and a bunch of per-snp stats & metadata.
        # genos = list(GT=matrix(...))

        # example code: return a list with a few stats
        if (is.null(genos)) {
            return(NULL)
        }
        list(
            n.snp = genos$nSNP,
            mean.maf = mean(genos$MAF),
            sample.dist = dist(genos$GT)
        )
    })

