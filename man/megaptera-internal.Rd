\name{megaptera-internal}

% these function deserve proper documentation

\alias{giConsensus}
\alias{stepBOLD}
\alias{supermatrix.mrbayes}
\alias{conn}
\alias{cropToReference} % compareToRef
\alias{tax}
\alias{loci}
\alias{pars}
\alias{addNbAcc}
\alias{addRanks} % ncbiTaxonomy
\alias{alignGenus} % stepG
\alias{alignSisterClades} % stepG
\alias{alignSpecies} % stepC
\alias{as.Linnean} % stepE
\alias{boldLineage} % dbUpdateTaxonomy
\alias{boxplotMSA}
\alias{checkExcludedSpec}
\alias{checkStatus}
\alias{check.TaxAgainstPhy}
\alias{compareToRef} % stepE
\alias{compare.phylo.backbone}
\alias{coverage}
\alias{dbAutoFillTaxonomy}
\alias{dbconnect}
\alias{dbChooseRef}
\alias{dbExcludeIndet}
\alias{dbSummaryMSA}
\alias{dbTableNames}
\alias{dbUpdateSurrogate}
\alias{dbWriteMSA}
\alias{dbPReadDNA}
\alias{dbProgress} % stepB, stepBOLD, stepC, stepD, stepE, stepF, stepMAFFT, stepPASTA, stepPOLENTA, stepH, stepI
\alias{dbReadReference}
\alias{dbUpdateLocus}
\alias{dbUpdateReference} % stepD
\alias{dbUpdateReference_Clade} % stepD
\alias{dbWriteDNAbin}
\alias{dbWriteString}
\alias{DNAbin2pg} % dbWriteMSA
\alias{downloadSequences}
\alias{EFetchHigherTaxon}
\alias{EFetchLocus}
\alias{EFetchXML}
\alias{extractIGS}
\alias{extractLocus}
\alias{filter.alignment}
\alias{fixTaxonomy} % deserves explicit doc
\alias{force.mp}
\alias{formatSpecList} % stepA
\alias{fuzzy.match}
\alias{getGeneRange}
\alias{getLineage}
\alias{getMaxDist}
\alias{htmlTable} % megaptera2Rmarkdown
\alias{insert.tip.phylo}
\alias{is.alignable}
\alias{is.Linnean} % stepE
\alias{MAD}
\alias{md5} % stepF
\alias{markersGenbank}
\alias{maxDist}
\alias{maxDistMPI}
\alias{meanPairwiseIdentity}
\alias{megapteraHeatmap}
\alias{myDist}
\alias{ncbiLineage}
\alias{NCBI.wrap}
\alias{nucleotideOverlap}
\alias{pg2DNAbin} % dbReadMSA
\alias{proParte} % prune.phylo.rank
\alias{revCompTest}
\alias{robustXMLparse} % downloadSequences
\alias{select.seqs}
\alias{slog}
\alias{specCons} % stepD
\alias{speciesConsensus}
\alias{sortRanks}
\alias{splitAlignment} % stepH
\alias{splitGiTaxon}
\alias{splitNonoverlapping} % splitAlignment
\alias{sql.conform}
\alias{sqlTaxonomyHeader}
\alias{sql.wrap}
\alias{summary.sproj}
\alias{tab}
\alias{tax2tree}
\alias{taxonomy2html} % megaptera2Rmarkdown
\alias{term} % downloadSequences
\alias{valid.NCBI.taxon}
\alias{whereToInsert}
\alias{wrapSQL}
\alias{write.dna.spectable}
\alias{XML2acc}

\title{Internal MEGAPTERA Functions}
\description{Internal \pkg{megaptera} functions.}
\note{
These are internal functions in \code{megaptera} and are not intended to be called by the user.
}
\keyword{internal}
