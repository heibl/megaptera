\name{locus}
\alias{locus}
\alias{locusRef}
\title{Locus/Phylogenetic Marker Definition}
\description{
Creates a S4 classes defining a phylogenetic marker for a \code{\link{megapteraProj}}.
}
\usage{
locus(..., not, search.fields = c("gene", "title"),
      use.genomes = TRUE, align.method = "auto", 
      min.identity = 0.75, min.coverage = 0.5, 
      check = FALSE)
   
locusRef(..., not, search.fields = c("gene", "title"),
      use.genomes = TRUE, align.method = "auto",
      min.identity = 0.75, min.coverage = 0.5, 
      reference, adj.gene1 = NULL, adj.gene2 = NULL,
      check = FALSE)
}
\arguments{
  \item{\dots}{
a vector of mode \code{character} giving strings that should be seached for; the first element is taken to name the corresponding postgreSQL tables.
}
  \item{not}{
a vector of mode \code{character} giving strings that should be excluded from the search results; corresponds to the use of \code{NOT} in GenBank query.
  }
  \item{search.fields}{
a vector of mode \code{character} setting the search fields (or attributes) of the Nucleotide database to be searched for the strings specified via the \dots argument.
  }
  \item{use.genomes}{
logical, if \code{TRUE} sequences of loci will be extracted from annotated genomes of chloroplasts and mitochondria. Due to the existance of different annotation styles, which might (still) be incompatible with \bold{megaptera}, this option can be used to turn the extraction of sequences from whole genomes off.
  }
  \item{align.method}{
a character string giving the alignment method in MAFFT. Available accuracy-oriented methods for less than 200 sequences are \code{"localpair"}, \code{"globalpair"}, and \code{"genafpair"}; \code{"retree 1"} and  \code{"retree 2"} are for speed-oriented alignment. The default is \code{"auto"}, which lets MAFFT choose an opproriate alignment method.
  }
  \item{min.identity}{
numeric between 0 and 1, giving the minimum proportion of nucleotides required to be identical with the reference sequence in order to be included in an alignment (default: 0.75)
  }
  \item{min.coverage}{
numeric between 0 and 1, giving the minimum proportion of nucleotide positions a sequence must have in common with the reference sequence in order to be included in an alignment (default: 0.5)
  }
  \item{reference}{an object of class \code{\link{DNAbin}} containing reference sequences. Alternatively, one or more GI numbers (currently only of chloroplast or mitochondrial genomes) can be given and the appropriate sequences will be extracted automatically.
  }
  \item{adj.gene1}{a vector of mode \code{character} giving strings to identify the upstream coding region of an intergenic spacer (IGS); to be effective the vector of strings given via \dots has to identify an IGS.
  }
  \item{adj.gene2}{a vector of mode \code{character} giving strings to identify the downstream coding region of an intergenic spacer (IGS); to be effective the vector of strings given via \dots has to identify an IGS.
  }
  \item{check}{logical: if \code{TRUE}, the existance of a locus as specified by \code{...} and \code{not} will be checked. }
}
\details{The website \url{http://www.ncbi.nlm.nih.gov/genome/browse/?report=5} or \code{\link{ncbiGenome}} can be used to find complete organelle genomes that can be used as references.}
\value{
an object of class \code{locus} or \code{locusRef}
}
%%\references{}
\author{Christoph Heibl}
%%\note{}
\seealso{
\code{\link{ncbiGenome}} helps to select representative organelle genomes as references. 

Use \code{\link{dbPars}}, \code{\link{taxon}}, and \code{\link{megapteraPars}} for defining of database parameters, taxa, and the pipeline's parameters, respectively, and \code{\link{megapteraProj}} for the bundling of input data.
}
%\examples{}
