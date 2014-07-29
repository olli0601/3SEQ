#
#	mtDNA example: 262 sequences
#
require(r3SEQ)
data(mtDNA)
#mtDNA data is stored in 'seq' DNAbin matrix object
print( seq )							
#write sequences to directory for processing
indir		<- getwd()
infile		<- 'mtDNA.R'
insignat	<- ''
save(seq, file=paste(indir, infile, sep='/'))	
#Run 3SEQ in a single batch job
r3seq.pipe.run.3seq(indir, infile, batch.n=1, hpc.walltime=1, hpc.q=NA, hpc.mem="500mb", hpc.nproc=1)
#Process 3SEQ output into R data.table
insignat	<- ''
argv		<<-	r3seq.cmd.process.3SEQ.output(indir, infile, insignat, resume=0, verbose=1) 
argv		<<- unlist(strsplit(argv,' '))
df.recomb	<- r3seq.prog.process.3SEQ.output()	
