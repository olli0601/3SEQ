#
#	Neisseria example: 4 sequences
#
require(r3SEQ)
data(neisseria)
#Neisseria data is stored in 'seq' DNAbin matrix object
print( seq )							
#write sequences to directory for processing
indir		<- getwd()
infile		<- 'neisseria.R'
insignat	<- ''
save(seq, file=paste(indir, infile, sep='/'))	
#The data contains 4 sequences, so run 3SEQ in a single batch job
r3seq.pipe.run.3seq(indir, infile, batch.n=1, hpc.walltime=1, hpc.q=NA, hpc.mem="500mb", hpc.nproc=1)
