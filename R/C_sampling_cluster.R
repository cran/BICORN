#' cis-Regulatory Module Sampling Function
#'
#' Function 'C_sampling_cluster' samples a candidate cis-regulatory module for each gene,
#' according to a discrete posterior probability distribution.
#' @param Y gene expression data matrix
#' @param C_old TF-gene binding network sampled from the previous round
#' @param A_old regulatory strength matrix sampled from the previous round
#' @param X_old transcription factor activity matrix sampled from the previous round
#' @param base_line_old gene expression baseline activity sampled from the previous round
#' @param C_prior prior TF-gene binding network
#' @param sigma_noise variance of gene expression fitting residuals
#' @param sigma_A variance of regulatory strength
#' @param sigma_baseline variance of gene expression baseline activity
#' @param sigma_X variance of transcription factor activity
#' @param BICORN_input this list structure contains TF symbols, Gene symbols and candidate modules
#' @import stats
#' @keywords Module sampling



C_sampling_cluster<-function(Y, C_old, A_old, X_old, base_line_old, C_prior, sigma_noise, sigma_A, sigma_baseline, sigma_X, BICORN_input){

  #Output
    #Sampled: list structure
      #Sampled$C: sampled TF-gene regulatory network
      #Sampled$A: sampled TF-gene regulatory strength (RS) matrix
      #Sampled$C_cluster: sampled module-gene regulatory network

  Num_genes=nrow(Y)# number of genes
  Num_samples=ncol(Y)# number of samples
  Num_TFs=ncol(C_prior)# number of TFs
  Num_clusters=nrow(BICORN_input$C_candidate_clusters)
  C=C_old
  A=A_old

  C_cluster=matrix(0, nrow=Num_genes, ncol=Num_clusters)

  for (n in 1:Num_genes){
    #gene
    post_probability=matrix(0, nrow=1,ncol=Num_clusters)
    A_candidate=matrix(0,nrow=Num_clusters,ncol=Num_TFs)

    for (s in 1:Num_clusters){#cluster values based on previous state
      for (t1 in 1:Num_TFs){
        if (C_old[n,t1]==1 && BICORN_input$C_candidate_clusters[s,t1]==1){
          A_candidate[s,t1]=A_old[n,t1]
        }else if (C_old[n,t1]==0 && BICORN_input$C_candidate_clusters[s,t1]==1){
          A_candidate[s,t1]=0
        }else if (C_old[n,t1]==1 && BICORN_input$C_candidate_clusters[s,t1]==0){
          A_candidate[s,t1]=0
        }else{
          A_candidate[s,t1]=0
        }
      }
    }


    for (s in 1:Num_clusters){# sample clusters
      for (t2 in 1:Num_TFs){
        temp_expression=matrix(0, nrow=1, ncol=Num_samples)
        for (m in 1:Num_samples){
          temp_expression[m]=sum(A_candidate[s,]*BICORN_input$C_candidate_clusters[s,]*X_old[,m])#based on previous binding state
        }

        if (C_old[n,t2]==1 && BICORN_input$C_candidate_clusters[s,t2]==1){
          temp_expression_c=temp_expression-A_candidate[s,t2]*BICORN_input$C_candidate_clusters[s,t2]*X_old[t2,]
          vairance_A=Num_samples*sigma_A*sigma_noise/(sum(X_old[t2,]^2)*sigma_A+Num_samples*sigma_noise)
          mean_A=sum((Y[n,]-temp_expression_c-base_line_old[n])*X_old[t2,])/sigma_noise/Num_samples*vairance_A
          probability_binding=-log(vairance_A)/2-(A_candidate[s,t2]-mean_A)^2/(2*vairance_A)+log(as.numeric(C_prior[n,t2]))

        }else if (C_old[n,t2]==0 && BICORN_input$C_candidate_clusters[s,t2]==1){
          temp_expression_c=temp_expression-A_candidate[s,t2]*BICORN_input$C_candidate_clusters[s,t2]*X_old[t2,]
          vairance_A=Num_samples*sigma_A*sigma_noise/(sum(X_old[t2,]^2)*sigma_A+Num_samples*sigma_noise)
          mean_A=sum((Y[n,]-temp_expression_c-base_line_old[n])*X_old[t2,])/sigma_noise/Num_samples*vairance_A
          A_candidate[s,t2]=rnorm(1, mean=mean_A, sd=sqrt(vairance_A))
          probability_binding=-log(vairance_A)/2-(A_candidate[s,t2]-mean_A)^2/(2*vairance_A)+log(as.numeric(C_prior[n,t2]))

        }else if (C_old[n,t2]==1 && BICORN_input$C_candidate_clusters[s,t2]==0){
          temp_expression_c=temp_expression
          vairance_A=Num_samples*sigma_A*sigma_noise/(sum(X_old[t2,]^2)*sigma_A+Num_samples*sigma_noise)
          mean_A=sum((Y[n,]-temp_expression_c-base_line_old[n])*X_old[t2,])/sigma_noise/Num_samples*vairance_A
          A_candidate[s,t2]=0
          probability_binding=-log(vairance_A)/2-(0-mean_A)^2/(2*vairance_A)+log(as.numeric(1-C_prior[n,t2]))

        }else {#C(n,t2)==1 && BICORN_input$C_candidate_clusters(s,t2)==0
          temp_expression_c=temp_expression
          vairance_A=Num_samples*sigma_A*sigma_noise/(sum(X_old[t2,]^2)*sigma_A+Num_samples*sigma_noise)
          mean_A=sum((Y[n,]-temp_expression_c-base_line_old[n])*X_old[t2,])/sigma_noise/Num_samples*vairance_A
          A_candidate[s,t2]=0
          probability_binding=-log(vairance_A)/2-(0-mean_A)^2/(2*vairance_A)+log(as.numeric(1-C_prior[n,t2]))
        }
        post_probability[s] = post_probability[s]+probability_binding#log format of probability
      }
      post_probability[s]=max(c(exp(post_probability[s]), 1e-6))
    }
    post_probability=post_probability/sum(post_probability)

    candidate_state=sample(Num_clusters, 1, prob=post_probability)#sample state according to probability
    C[n,]=BICORN_input$C_candidate_clusters[candidate_state,]
    A[n,]=A_candidate[candidate_state,]
    C_cluster[n,candidate_state]=1
  }
  Sampled=list('C'=C,
               'A'=A,
               'C_cluster'=C_cluster)
  return(Sampled)
}
