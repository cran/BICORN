#' Regulation Strength Sampling Function
#'
#' Function 'A_sampling' estimates a regulation strength for each sampled binding event in C,
#' according to a posterior Gaussian distribution.
#' @param Y gene expression data matrix
#' @param C sampled TF-gene binding network
#' @param A_old regulatory strength sampled from the previous round,
#' used as a prior in current function
#' @param X sampled transcription factor activity matrix
#' @param base_line sampled gene expression baseline activity
#' @param C_prior prior TF-gene binding network
#' @param sigma_noise variance of gene expression fitting residuals
#' @param sigma_A variance of regulatory strength
#' @param sigma_baseline variance of gene expression baseline activity
#' @param sigma_X variance of transcription factor activity
#' @import stats
#' @keywords Regulation_Strength


A_sampling<-function(Y, C, A_old, X, base_line, C_prior, sigma_noise, sigma_A, sigma_baseline, sigma_X){
  Num_genes=nrow(Y)
  Num_samples=ncol(Y)
  Num_TFs=ncol(C_prior)
  A_new=matrix(0, nrow=Num_genes, ncol=Num_TFs)

  for (n in 1:Num_genes){
    for (t2 in 1:Num_TFs){
      temp=matrix(0, nrow=1, ncol=Num_samples)
      for (m in 1:Num_samples){
        for (t1 in 1:Num_TFs){
          temp[m]=temp[m]+A_old[n,t1]*C[n,t1]*X[t1,m]
        }
      }
      if (C[n,t2]==1){
        temp_c=temp-A_old[n,t2]*X[t2,m]
        vairance_A=Num_samples*sigma_A*sigma_noise/(sum(X[t2,]^2)*sigma_A+Num_samples*sigma_noise)
        mean_A=sum((Y[n,]-temp_c-base_line[n])*X[t2,])/sigma_noise/Num_samples*vairance_A
        A_temp=rnorm(1, mean=mean_A, sd=sqrt(vairance_A))
        A_new[n,t2] = A_temp
      }else
        A_new[n,t2]=0
    }
  }
  return(A_new)
}
