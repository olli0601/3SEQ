#
#	mtDNA example: 262 sequences
#
require(r3SEQ)
data(mtDNA)
#	mtDNA data is stored in 'seq' DNAbin matrix object
print( seq )							
#	write sequences to directory for processing
indir		<- getwd()
infile		<- 'mtDNA.R'
insignat	<- ''
save(seq, file=paste(indir, infile, sep='/'))	
#	Run 3SEQ in a single batch job
r3seq.pipe.run.3seq(indir, infile, batch.n=1, hpc.walltime=1, hpc.q=NA, hpc.mem="500mb", hpc.nproc=1)
#	Produce several shell files with a PBS header that can be submitted to an HPC system.
#	At the moment, only one PBS header is implemented. This PBS header is configured for the Imperial HPC system.
#	To interface with a different HPC system, a new PBS header function must be created. See 'cmd.hpcwrapper' and 'cmd.hpcwrapper.cx1.ic.ac.uk' for examples.
hpc.sys		<- "cx1.hpc.ic.ac.uk"
r3seq.pipe.run.3seq(indir, infile, batch.n=4, hpc.sys=hpc.sys, hpc.walltime=35, hpc.q=NA, hpc.mem="3850mb", hpc.nproc=1, verbose=1)
