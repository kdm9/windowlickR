require(foreach)

#' windowlickr
#'
#' @param bcf Indexed BCF or VCF file
#' @param func Function to evaulate for each window. is called as func(window.metadata, window.genotypes) for each window. Func ideally returns a tibble or similar.
#' @param windowsize Size of each window (required unless `windows` provided)
#' @param slide Number of bases to skip between each window start (default `windowsize`)
#' @param samples Set of samples to calculate LD among
#' @param minMAF Minimum SNP minor allle freq
#' @param maxMissing Maximum SNP missing data rate
#' @param windows Restrict analyses to these windows (expected to be created with bcf_getWindows)
#' @param chrom Only generate windows on contig 'chrom'
#' @param from Only generate windows starting at `from` on `chrom`
#' @param to Only generate windows until `to` on `chrom`
#' @param errors One of "remove" or "stop", see ?foreach::foreach 's .errorhandling parameter
#'
#' @return list One list per window, with slots `regions`, `contig`, `start`, `stop`, `result`, where result is the return value of `func` applied to this window.
#' @export
windowlickr = function (bcf, func, windowsize=NULL, slide=windowsize, samples=NULL, minMAF=0.1, maxMissing=0.8,
                             windows=NULL, chrom=NULL, from=NULL, to=NULL, errors="stop", export=NULL, pkgs=NULL) {
  if (is.null(windows)) {
    if (is.null(windowsize) || is.null(slide)) stop("Invalid or missing windowsize/slide")
    cat("Making windows\n")
    windows = bcf_getWindows(bcf, windowsize=windowsize, slide = slide,
                             chrom=chrom, from=from, to=to)
    cat(paste("Using", nrow(windows), "generated windows\n"))
  } else {
    cat(paste("Using", nrow(windows), "user-supplied windows\n"))
  }
  window_i = seq_len(nrow(windows))
  export=c("bcf", "windows", "minMAF", "maxMissing", "samples", "func", export)
  pkgs = c("magrittr", pkgs)
  result = foreach(i=window_i, .export = export, .packages=pkgs, .errorhandling=errors) %dopar% {
      window = as.list(windows[i, ])
      geno =  bcf_getGT(bcf, window$region, minMAF=minMAF, maxMissing=maxMissing,
                          samples=samples)
      window$result = func(window, geno)
      window
  }
  result
}
