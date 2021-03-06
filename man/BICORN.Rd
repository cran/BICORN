% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BICORN.function.R
\name{BICORN}
\alias{BICORN}
\title{BICORN Algorithm Function}
\usage{
BICORN(BICORN_input = NULL, L = 100, output_threshold = 10)
}
\arguments{
\item{BICORN_input}{this list structure contains TF symbols, gene symbols and candidate modules}

\item{L}{total rounds of Gibbs Sampling.}

\item{output_threshold}{number of rounds after which we start to record results.}
}
\description{
Function 'BICORN' infers a posterior module-gene regulatory network by iteratively sampling regulatory strength,
transcription factor activity and several key model parameters.
}
\examples{

# load in the sample data input
data("sample.input")

# Data initialization (Integerate prior binding network and gene expression data)
BICORN_input<-data_integration(Binding_matrix = Binding_matrix, Binding_TFs = Binding_TFs,
Binding_genes = Binding_genes, Exp_data = Exp_data, Exp_genes = Exp_genes,
Minimum_gene_per_module_regulate = 2)

# Infer cis-regulatory modules (TF combinations) and their target genes
BICORN_output<-BICORN(BICORN_input, L = 2, output_threshold = 1)
}
\keyword{BICORN}
\keyword{sampling}
