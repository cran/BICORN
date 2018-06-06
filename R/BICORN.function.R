#' BICORN Algorithm Function
#'
#' Function 'BICORN' infers a posterior module-gene regulatory network by iteratively sampling regulatory strength,
#' transcription factor activity and several key model parameters.
#' @param BICORN_input this list structure contains TF symbols, gene symbols and candidate modules
#' @param L total rounds of Gibbs Sampling.
#' @param output_threshold number of rounds after which we start to record results.
#' @keywords BICORN sampling
#' @export
#' @examples
#'
#' # load in the sample data input
#' data("sample.input")
#'
#' # Data initialization (Integerate prior binding network and gene expression data)
#' BICORN_input<-data_integration(Binding_matrix = Binding_matrix, Binding_TFs = Binding_TFs,
#' Binding_genes = Binding_genes, Exp_data = Exp_data, Exp_genes = Exp_genes,
#' Minimum_gene_per_module_regulate = 2)
#'
#'# Infer cis-regulatory modules (TF combinations) and their target genes
#'BICORN_output<-BICORN(BICORN_input, L = 2, output_threshold = 1)
#


BICORN<-function(BICORN_input = NULL, L = 100, output_threshold = 10){
  Y=BICORN_input$expression_data
  C_prior=BICORN_input$gene_cluster_mapping%*%BICORN_input$C_candidate_clusters
  C_prior[which(C_prior>0)]=0.9
  C_prior[which(C_prior==0)]=0.1

  #*******************************known value********************************
  alpha=1
  beta=1
  sigma_A=1
  sigma_baseline=1
  sigma_X=1

  Num_genes=nrow(Y)# number of genes
  Num_samples=ncol(Y)# number of samples
  Num_TFs=ncol(C_prior)# number of TFs
  Num_clusters=nrow(BICORN_input$C_candidate_clusters)

  #**********************initialize the baseline activity********************
  base_line_old=matrix(0, nrow=1,ncol=Num_genes)
  for (n in 1:Num_genes){
    base_line_old[1,n] = rnorm(1, mean=0, sd=sqrt(sigma_baseline))
  }

  #**********************initialize binding matrix and binding strength******
  C_old=matrix(0, nrow=Num_genes, ncol=Num_TFs)
  A_old=matrix(0, nrow=Num_genes, ncol=Num_TFs)
  threshold=0.5

  for (n in 1:Num_genes){
    for (t in 1:Num_TFs){
      if (C_prior[n,t]>threshold){
        C_old[n,t]=1
        A_old[n,t]=rnorm(1, mean=0, sd=sqrt(sigma_A))
      }else{
        C_old[n,t]=0
        A_old[n,t]=0
      }
    }
  }

  #**********************initialize TF activity******************************
  X_old=matrix(0, nrow=Num_TFs, ncol=Num_samples)
  for (t in 1:Num_TFs){
    for (m in 1:Num_samples){
      X_old[t,m]=rnorm(1, mean=0, sd=sqrt(sigma_X))
    }
  }

  #*************************Gibbs sampling***********************************
  sigma_noise_summition=matrix(0, nrow=L, ncol=1)
  C_summition=matrix(0, nrow=Num_genes,ncol=Num_TFs)
  C_cluster_summition=matrix(0, nrow=Num_genes,ncol=Num_clusters)
  A_summition=matrix(0, nrow=Num_genes,ncol=Num_TFs)
  X_summition=matrix(0, nrow=Num_TFs,ncol=Num_samples)

  for (i in 1:L){
    cat(paste("start", i, "round of sampling, in total", L, "rounds!\n"))
    #****sample noise variance *****************
    sigma_noise_new<-sigmanoise_sampling(Y, C_old, A_old, X_old, base_line_old, C_prior, sigma_A, sigma_baseline, sigma_X, alpha, beta)

    #****sample baseline expression ************
    base_line_new<-baseline_sampling(Y, C_old, A_old, X_old, base_line_old, C_prior, sigma_noise_new, sigma_A, sigma_baseline, sigma_X)

    #****sample binding matrix C and update A with C_old=0 if C_new=1******
    Sampled<-C_sampling_cluster(Y, C_old, A_old, X_old, base_line_new, C_prior, sigma_noise_new, sigma_A, sigma_baseline, sigma_X, BICORN_input)

    #****sample bidning strength matrix A ******
    A_new<-A_sampling(Y, Sampled$C, Sampled$A, X_old, base_line_new, C_prior, sigma_noise_new, sigma_A, sigma_baseline, sigma_X)

    #****sample TFA matrix X *******************
    X_new<-X_sampling(Y, Sampled$C, A_new, X_old, base_line_new, C_prior, sigma_noise_new, sigma_A, sigma_baseline, sigma_X)

    C_old=Sampled$C
    A_old=A_new
    X_old=X_new
    base_line_old=base_line_new

    if (i>output_threshold){
      C_summition=C_summition+Sampled$C
      C_cluster_summition=C_cluster_summition+Sampled$C_cluster
      A_summition=A_summition+A_new
      X_summition=X_summition+X_new
    }
  }

  C_cluster_summition=C_cluster_summition/(L-output_threshold)
  C_summition=C_summition/(L-output_threshold)

  BICORN_output=list('Posterior_module_gene_regulatory_network'=C_cluster_summition,
                     'Posterior_TF_gene_regulatory_network'=C_summition,
                     'Modules'=BICORN_input$C_candidate_clusters,
                     'TFs'=BICORN_input$TFs,
                     'Genes'=BICORN_input$genes)

  colnames(BICORN_output$Modules) <- BICORN_output$TFs
  rownames(BICORN_output$Modules) <- c(1:nrow(BICORN_output$Modules))

  colnames(BICORN_output$Posterior_module_gene_regulatory_network) <- c(1:nrow(BICORN_output$Modules))
  rownames(BICORN_output$Posterior_module_gene_regulatory_network) <- BICORN_output$Genes

  colnames(BICORN_output$Posterior_TF_gene_regulatory_network) <- BICORN_output$TFs
  rownames(BICORN_output$Posterior_TF_gene_regulatory_network) <- BICORN_output$Genes

  return(BICORN_output)
}
