#' This file contains parallel processing scripts that generate and submit PBS scripts to a high performance system. 

######################################################################################
pipeline.recom.run.3seq<- function()
{
	indir		<- paste(DATA,"tmp",sep='/')		
	infile		<- "ATHENA_2013_03_NoDRAll+LANL_Sequences"
	#infile		<- "ATHENA_2013_03_NoDRAll+LANL_Sequences100"
	insignat	<- "Thu_Aug_01_17/05/23_2013"
	
	batch.n			<- 100		#for 1e4 sequences, does about 100 in 25hrs, so request 100 batches for walltime 25 expected + 10hrs grace
	file			<- paste(indir,'/',infile,'_',gsub('/',':',insignat),".R",sep='')
	load(file)
	batch.seq		<- round(seq.int(0,nrow(seq.PROT.RT),len=batch.n),d=0)
	batch.seq		<- rbind(batch.seq[-length(batch.seq)], batch.seq[-1]-1)
	#batch.seq		<- batch.seq[,1:10]	#test run
	file			<- paste(indir,'/',infile,'_',gsub('/',':',insignat),".phylip",sep='')
	lapply(seq_len(ncol(batch.seq)),function(j)
			{					
				cmd			<- cmd.recombination.run.3seq(infile=file, outfile=paste(indir,'/',infile,'_',batch.seq[1,j],'-',batch.seq[2,j],'_',gsub('/',':',insignat),".3seq",sep=''), recomb.3seq.siglevel=0.1, nproc=1, recomb.3seq.testvsall.beginatseq=batch.seq[1,j], recomb.3seq.testvsall.endatseq=batch.seq[2,j], verbose=1)
				cmd			<- cmd.hpcwrapper(cmd, hpc.walltime=35, hpc.q="pqeph", hpc.mem="3850mb",  hpc.nproc=1)
				cat(cmd)
				outdir		<- indir
				outfile		<- paste("r3seq",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='.')									
				cmd.hpccaller(outdir, outfile, cmd)			
			})	
}
######################################################################################
pipeline.recom.get.phyloincongruence.for.candidates<- function()
{	
	indir		<- paste(DATA,"tmp",sep='/')		
	infile		<- "ATHENA_2013_03_NoDRAll+LANL_Sequences"		
	insignat	<- "Thu_Aug_01_17/05/23_2013"
	resume		<- 0
	verbose		<- 1
	
	argv				<<-	cmd.recombination.process.3SEQ.output(indir, infile, insignat, resume=1, verbose=1) 
	argv				<<- unlist(strsplit(argv,' '))
	df.recomb			<- prog.recom.process.3SEQ.output()	
	
	triplets			<- seq_len(nrow(df.recomb))
	triplets			<- 147:nrow(df.recomb)
	dummy	<- lapply(triplets, function(i)
			{					
				if(verbose)	cat(paste("\nprocess triplet number",i,"\n"))
				argv				<<- cmd.recombination.check.candidates(indir, infile, insignat, i, resume=resume, verbose=1)
				argv				<<- unlist(strsplit(argv,' '))
				prog.recom.get.incongruence()		#this starts ExaML for the ith triplet			
			})	
}
######################################################################################
pipeline.recom.plot.phyloincongruence.for.candidates<- function()
{	
	indir		<- paste(DATA,"tmp",sep='/')		
	infile		<- "ATHENA_2013_03_NoDRAll+LANL_Sequences"		
	insignat	<- "Thu_Aug_01_17/05/23_2013"
	resume		<- 0
	verbose		<- 1
	
	argv				<<- cmd.recombination.plot.incongruence(indir, infile, insignat, triplet.id=NA, opt.select="ng2", verbose=1)
	argv				<<- unlist(strsplit(argv,' '))
	prog.recom.plot.incongruence()		
	
	argv				<<- cmd.recombination.plot.incongruence(indir, infile, insignat, triplet.id=NA, opt.select="g2", verbose=1)
	argv				<<- unlist(strsplit(argv,' '))
	prog.recom.plot.incongruence()	
}
######################################################################################
pipeline.recom<- function()
{
	if(0)	pipeline.recom.run.3seq()
	if(0)	pipeline.recom.get.phyloincongruence.for.candidates()
	if(0)	pipeline.recom.plot.phyloincongruence.for.candidates()
}



