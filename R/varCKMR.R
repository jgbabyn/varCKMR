#' Simulate from markers dataframe
#'
#' @param markers dataframe with the marker information geno,locus then freqs
#' @param n The number to simulate.
#' @export
simulate_markers <- function(markers,n){
    samp = tapply(1:nrow(markers),markers$locus,
                  FUN=function(x,n){
                      geno = markers$geno[x]
                      prob = markers$freqs[x]
                      y= sample(geno,2*n,TRUE,prob)
                      y},n=n)
    samp2 = lapply(seq_along(samp),function(x){
        larp = data.table(matrix(samp[[x]],ncol=2))
        names(larp) = paste0("m_",names(samp[x]),"_",c("A","B"))
        larp
    })
    do.call(cbind,samp2)    
}
