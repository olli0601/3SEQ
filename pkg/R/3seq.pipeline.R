######################################################################################
hivc.pipeline.recombination<- function()
{
	if(0)	#generate candidate recombinants	- 	this is computationally expensive
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
		stop()
	}
	if(1)	#validate candidate recombinants
	{
		indir		<- paste(DATA,"tmp",sep='/')		
		infile		<- "ATHENA_2013_03_NoDRAll+LANL_Sequences"		
		insignat	<- "Thu_Aug_01_17/05/23_2013"
		resume		<- 0
		verbose		<- 1
		
		argv				<<-	cmd.recombination.process.3SEQ.output(indir, infile, insignat, resume=1, verbose=1) 
		argv				<<- unlist(strsplit(argv,' '))
		df.recomb			<- hivc.prog.recombination.process.3SEQ.output()	
		
		triplets			<- seq_len(nrow(df.recomb))
		triplets			<- 147:nrow(df.recomb)
		dummy	<- lapply(triplets, function(i)
				{					
					if(verbose)	cat(paste("\nprocess triplet number",i,"\n"))
					argv				<<- cmd.recombination.check.candidates(indir, infile, insignat, i, resume=resume, verbose=1)
					argv				<<- unlist(strsplit(argv,' '))
					hivc.prog.recombination.check.candidates()		#this starts ExaML for the ith triplet			
				})
		stop()
	}
	if(1)	#process validation step to produce plots
	{
		indir		<- paste(DATA,"tmp",sep='/')		
		infile		<- "ATHENA_2013_03_NoDRAll+LANL_Sequences"		
		insignat	<- "Thu_Aug_01_17/05/23_2013"
		resume		<- 0
		verbose		<- 1
		
		argv				<<- cmd.recombination.plot.incongruence(indir, infile, insignat, triplet.id=NA, opt.select="ng2", verbose=1)
		argv				<<- unlist(strsplit(argv,' '))
		hivc.prog.recombination.plot.incongruence()		
		
		argv				<<- cmd.recombination.plot.incongruence(indir, infile, insignat, triplet.id=NA, opt.select="g2", verbose=1)
		argv				<<- unlist(strsplit(argv,' '))
		hivc.prog.recombination.plot.incongruence()
		
	}
	if(0)	#	collect likely recombinants or those likely confounding the phylogeny		-- 	identified by eye
	{
		verbose				<- 1
		
		argv				<<-	cmd.recombination.process.3SEQ.output(indir, infile, insignat, resume=1, verbose=1) 
		argv				<<- unlist(strsplit(argv,' '))
		df.recomb			<- hivc.prog.recombination.process.3SEQ.output()
		setnames(df.recomb, "dummy", "triplet.id")
		setkey(df.recomb, triplet.id)
		#	collect likely recombinants or those likely confounding the phylogeny
		recombinants.ng2	<-	c(	subset( df.recomb, triplet.id==60 )[,child],		#likely confounding
									subset( df.recomb, triplet.id==52 )[,child],		
									subset( df.recomb, triplet.id==55 )[,parent2],"R09-26706","R11-12152","R12-15108",		#others cluster in addition
									subset( df.recomb, triplet.id==48 )[,parent2],"2007G319",								#others cluster in addition
									subset( df.recomb, triplet.id==112 )[,child],
									subset( df.recomb, triplet.id==129 )[,child],
									subset( df.recomb, triplet.id==148 )[,child],		#length only 50
									subset( df.recomb, triplet.id==135 )[,child],
									subset( df.recomb, triplet.id==102 )[,child],		#likely confounding
									subset( df.recomb, triplet.id==85 )[,child],
									subset( df.recomb, triplet.id==81 )[,child],		#likely confounding
									subset( df.recomb, triplet.id==125 )[,child],		#likely confounding but might be recombinant
									subset( df.recomb, triplet.id==120 )[,child]	)
		recombinants.g2		<-	c(	"2006G052", "2007G263", "M3621708072010", "M4048713072011", 
									"M4203226082011", "R03-14636", "TN_B.HT.2005.05HT_129336.EU439719", "TN_B.HT.2004.04HT_129732.EU439728", "TN_B.PH.2008.08R_01_361.AB587101", "TN_B.ZA.2011.patient_1720_seq_1746.KC423374", "TN_B.ZA.2012.patient_1720_seq_2734.KC423805", "TN_B.DO.2008.HIV_PRRT_PJ01967_48.JN713614",					#these cluster with  M4048713072011									
									"R11-15440", "R12-00343", "R10-09812",				#these are children of R08-20970
									"R12-07939" )
		recombinants		<- c( recombinants.ng2, recombinants.g2 )
		if(verbose)	cat(paste("\ncollected recombinants or likely confounding sequences, n=",length(recombinants)))
		#
		#	save non-recombinant sequence dataset
		#
		indir		<- paste(DATA,"tmp",sep='/')		
		infile		<- "ATHENA_2013_03_NoDRAll+LANL_Sequences"	
		insignat	<- "Thu_Aug_01_17/05/23_2013"
		outfile		<- "ATHENA_2013_03_NoRCDRAll+LANL_Sequences"	
		outsignat	<- "Fri_Nov_01_16/07/23_2013"
		
		file		<- paste(indir,'/',infile,'_',gsub('/',':',insignat),".R",sep='')
		load(file)
		tmp			<- setdiff(rownames(seq.PROT.RT), recombinants)
		seq.PROT.RT	<- seq.PROT.RT[tmp,]
		if(verbose)	cat(paste("\nnumber of sequences without (likely) recombinants, n=",nrow(seq.PROT.RT)))
		file		<- paste(indir,'/',outfile,'_',gsub('/',':',outsignat),".R",sep='')
		if(verbose)	cat(paste("\nsave new file to",file))
		save(seq.PROT.RT, file=file)
	}
}



