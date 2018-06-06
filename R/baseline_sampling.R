#' Gene Baseline Expression Sampling Function
#'
#' Function 'baseline_sampling' estimates a baseline expression for each gene,
#'according to a posterior Gaussian distribution.
#' @param Y gene expression data matrix
#' @param C sampled TF-gene binding network
#' @param A sampled regulatory strength matrix
#' @param X sampled transcription factor activity matrix
#' @param base_line_old prior gene expression baseline activity
#' @param C_prior prior TF-gene binding network
#' @param sigma_noise variance of gene expression fitting residuals
#' @param sigma_A variance of regulatory strength
#' @param sigma_baseline variance of gene expression baseline activity
#' @param sigma_X variance of transcription factor activity
#' @import stats
#' @keywords Baseline sampling



baseline_sampling<-function(Y, C, A, X, base_line_old, C_prior, sigma_noise, sigma_A, sigma_baseline, sigma_X){
  Num_genes=nrow(Y)
  Num_samples=ncol(Y)
  Num_TFs=ncol(C_prior)
  baseline_new=matrix(0, nrow=1, ncol=Num_genes)
  for (n in 1:Num_genes){
    temp=matrix(0, nrow=1, ncol=Num_samples)
    for (m in 1:Num_samples){
      for (t in 1:Num_TFs){
        temp[m]=temp[m]+A[n,t]*X[t,m]
      }
    }
    mean_baseline = sum(Y[n,]-temp)/Num_samples*sigma_baseline/(sigma_noise+sigma_baseline)
    variance_baseline = sigma_baseline*sigma_baseline/(sigma_noise+sigma_baseline)
    baseline_new[n]=rnorm(1, mean=mean_baseline, sd=sqrt(variance_baseline))
  }
  return(baseline_new)
}
