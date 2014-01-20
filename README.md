recombination.analyzer
======================

Tools to identify recombinant sequences among up to 100.000 sequences based on the 3SEQ program.

To install this R package, please install first the packages 'roxygenize2' 'data.table', and 'big.phylo' from https://github.com/olli0601/big.phylo.
Then type on the command line 'R CMD build pkg' and then 'R CMD INSTALL recombination.analyzer_1.0-0.tar.gz'

Help files and documentation are available once the package is loaded in R with 'require(recombination.analyzer)': 
To see how candidate recombinant sequences can be identified, type '?pipeline.recom.run.3seq',
To process job output into a data.table format, type '?prog.recom.process.3SEQ.output'.