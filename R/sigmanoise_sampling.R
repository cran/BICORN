#' Fitting Residule Variance Sampling Function
#'
#' Function 'sigmanoise_sampling' estimates the variance of overal gene expression fitting residuals,
#'according to an inverse-gamma distribution.
#' @param Y gene expression data matrix
#' @param C sampled TF-gene binding network
#' @param A sampled regulatory strength matrix
#' @param X sampled transcription factor activity matrix
#' @param base_line sampled gene expression baseline activity
#' @param C_prior prior TF-gene binding network
#' @param sigma_A variance of regulatory strength
#' @param sigma_baseline variance of gene expression baseline activity
#' @param sigma_X variance of transcription factor activity
#' @param alpha hyper-parameter for inverse-gamma distribution
#' @param beta hyper-parameter for inverse-gamma distribution
#' @import stats
#' @keywords Noise sampling


sigmanoise_sampling<-function(Y, C, A, X, base_line, C_prior, sigma_A, sigma_baseline, sigma_X, alpha, beta){
  Num_genes=nrow(Y)# number of genes
  Num_samples=ncol(Y)# number of samples
  Num_TFs=ncol(C_prior)# number of TFs
  alpha_new=alpha+1/2
  beta_temp=0
  for (n in 1:Num_genes){
    temp=matrix(0, nrow=1,ncol=Num_samples)
    for (m in 1:Num_samples){
      for (t in 1:Num_TFs){
        temp[m]=temp[m]+A[n,t]*X[t,m]
      }
    }

    temp_y=Y[n,]-temp-base_line[n]
    temp_y=temp_y^2
    beta_temp=beta_temp+sum(temp_y)/2
  }

  beta_new=beta+beta_temp/(Num_samples*Num_genes)

  sigmanoise_new=1/rgamma(1, alpha_new,1/beta_new)
  #sigmanoise_new=beta_new/rchisq(1, 2*alpha_new)
  return(sigmanoise_new)
}
