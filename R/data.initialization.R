#' Data Initialization for BICORN
#'
#' Function 'data_integration' integrates the prior TF-gene binding network and gene expression data together. It will remove any genes missing either TF bindings or gene expression and identify a list of candidate cis-regulatory modules.

#' @param Binding_matrix loaded prior binding network
#' @param Binding_TFs loaded transcription factors
#' @param Binding_genes loaded genes in the prior binding network
#' @param Exp_data loaded properly normalized gene expression data
#' @param Exp_genes loaded genes in the gene expression data
#' @param Minimum_gene_per_module_regulate the minimum number of genes regulated by each module, used for candidate module filtering.
#' @keywords Initialization
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


data_integration<-function(Binding_matrix = NULL, Binding_TFs = NULL, Binding_genes = NULL,
                           Exp_data, Exp_genes = NULL, Minimum_gene_per_module_regulate = 2){

  match_list<-match(Binding_genes, Exp_genes)
  index=which(match_list>0)
  match_list=match_list[index]
  Sampling_window_flag<-Binding_matrix[index,c(1:length(Binding_TFs))]

  Cluster_pattern=unique(Sampling_window_flag)

  Candidate_cluster_num=matrix(0, nrow=nrow(Cluster_pattern), ncol=1)
  index_temp=which(rowSums(abs(sweep(Sampling_window_flag,2,Cluster_pattern[1,])))==0)
  if (sum(index_temp)>0){
    Candidate_cluster_num[1,1]=length(index_temp)
  }

  for (k in 2:nrow(Cluster_pattern)){
    index_temp=which(rowSums(abs(sweep(Sampling_window_flag,2,Cluster_pattern[k,])))==0)
    if (sum(index_temp)>0){
      Candidate_cluster_num[k,1]=length(index_temp)
    }
  }
  C_candidate_clusters=Cluster_pattern[which(Candidate_cluster_num>=Minimum_gene_per_module_regulate),]
  TF_index=which(colSums(C_candidate_clusters)>0)
  Candidate_TFs=Binding_TFs[TF_index]
  C_candidate_clusters=C_candidate_clusters[,TF_index]
  Sampling_window_flag<-Sampling_window_flag[,TF_index]
  C_candidate_clusters=rbind(C_candidate_clusters, 0)


  Num_regions=nrow(Sampling_window_flag)
  Num_clusters=nrow(C_candidate_clusters)

  Cluster_mapping=matrix(0, nrow=Num_regions, ncol=Num_clusters)
  for (k in 1:Num_regions){
    for (c in 1:Num_clusters){
      if (sum(abs(C_candidate_clusters[c,]-Sampling_window_flag[k,]*C_candidate_clusters[c,]))==0){
        Cluster_mapping[k,c]=1
      }
    }
  }

  sig_gene_index=unique(match_list)
  genes=Exp_genes[sig_gene_index]
  expression_data=Exp_data[sig_gene_index,]
  TFs=Candidate_TFs
  gene_cluster_mapping=matrix(0, nrow=length(sig_gene_index), ncol=Num_clusters)
  for (k in 1:length(sig_gene_index)){
    index_temp=which(match_list==sig_gene_index[k])
    if (length(index_temp)==1){
      gene_cluster_mapping[k,]=Cluster_mapping[index_temp,]
    }else{
      gene_cluster_mapping[k,]=colSums(Cluster_mapping[index_temp,])
    }
  }
  gene_cluster_mapping[which(gene_cluster_mapping>0)]=1
  BICORN_input=list('genes'=genes,
                    'expression_data'=expression_data,
                    'TFs'=TFs,
                    'gene_cluster_mapping'=gene_cluster_mapping,
                    'C_candidate_clusters'=C_candidate_clusters)
  return(BICORN_input)
}
