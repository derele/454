


find.primer.sequence <- function(primer){
  cat("on reverse complement",
      names(fullest.assembly)[grep (revcom(primer), fullest.assembly,)], "\n")
  cat("on contig sequence", 
      names(fullest.assembly)[grep (primer, fullest.assembly,)], "\n")
}
