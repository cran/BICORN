#' Transcription Factor Activity Sampling Function
#'
#' Function 'X_sampling' estimates the hidden activities of each transcription factor,
#'according to a posterior Gaussian random process.
#' @param Y gene expression data matrix
#' @param C sampled TF-gene binding network
#' @param A sampled regulatory strength matrix
#' @param X_old sampled transcription factor activity matrix from the previous round
#' @param base_line sampled gene expression baseline activity
#' @param C_prior prior TF-gene binding network
#' @param sigma_noise variance of gene expression fitting residuals
#' @param sigma_A variance of regulatory strength
#' @param sigma_baseline variance of gene expression baseline activity
#' @param sigma_X variance of transcription factor activity
#' @import stats
#' @keywords TFA sampling


X_sampling<-function(Y, C, A, X_old, base_line, C_prior, sigma_noise, sigma_A, sigma_baseline, sigma_X){
  Num_genes=nrow(Y)
  Num_samples=ncol(Y)
  Num_TFs=ncol(C_prior)
  X_new=matrix(0, nrow=Num_TFs, ncol=Num_samples)
  for (t2 in 1:Num_TFs){
    for (m in 1:Num_samples){
      temp_y=0
      for (n in 1:Num_genes){
        temp=0
        for (t1 in 1:Num_TFs){
          temp=temp+A[n,t1]*C[n,t1]*X_old[t1,m]
        }
        temp=temp-A[n,t2]*C[n,t2]*X_old[t2,m]
        temp_y=temp_y+(Y[n,m]-temp-base_line[n])*A[n,t2]*C[n,t2]
      }
      variance_X=sigma_X*sigma_noise*Num_genes/(sum(A[,t2]^2)*sigma_X+Num_genes*sigma_noise)
      mean_X=temp_y/sigma_noise/Num_genes*variance_X
      X_new[t2,m]=rnorm(1, mean=mean_X, sd=sqrt(variance_X))
    }
  }
  return(X_new)
}
