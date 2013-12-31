if(!exists("HIVC.CODE.HOME"))	
{	
	HIVC.CODE.HOME	<- getwd()
	INST			<- paste(HIVC.CODE.HOME,"inst",sep='/')
}

#' @export
PR.BLASTMASK	<- "windowmasker"

#' @export
PR.BLASTMAKEDB	<- "makeblastdb"

#' @export
PR.BLASTN		<- "blastn"

#' @export
PR.CLUSTALO		<- "clustalo"

#' @export
PR.CLUSTALO.HMM	<- paste(INST,"align_HIV-1_pol_DNA.hmm",sep='/')

#' @export
PR.FIRSTSEQ		<- paste(HIVC.CODE.HOME,"pkg/misc/hivclu.startme.R -exeFIRSTSEQ",sep='/')

#' @export
PR.GENDISTMAT	<- paste(HIVC.CODE.HOME,"pkg/misc/hivclu.startme.R -exeGENDISTMAT",sep='/')

#' @export
PR.PRECLUST		<- paste(HIVC.CODE.HOME,"pkg/misc/hivclu.startme.R -exePRECLUST",sep='/')

#' @export
PR.CLUST		<- paste(HIVC.CODE.HOME,"pkg/misc/hivclu.startme.R -exeCLUST",sep='/')

#' @export
PR.CLUSTTPTN	<- paste(HIVC.CODE.HOME,"pkg/misc/hivclu.startme.R -exeCLUSTTPTN",sep='/')

#' @export
PR.CLUSTMSM		<- paste(HIVC.CODE.HOME,"pkg/misc/hivclu.startme.R -exeCLUSTMSM",sep='/')

#' @export
PR.EXAML.BSCREATE	<- paste(HIVC.CODE.HOME,"pkg/misc/hivclu.startme.R -exeBOOTSTRAPSEQ",sep='/')

#' @export
PR.RECOMB.3SEQ	<- system.file(package="hivclust", "ext", "3seq") 

#' @export
PR.RECOMB.PROCESS3SEQOUTPUT	<- paste(HIVC.CODE.HOME,"pkg/misc/hivclu.startme.R -exeRECOMB.PROCESS3SEQOUT",sep='/')

#' @export
PR.RECOMB.CHECKCANDIDATES	<- paste(HIVC.CODE.HOME,"pkg/misc/hivclu.startme.R -exeRECOMB.CHECKCANDIDATES",sep='/')

#' @export
PR.RECOMB.PLOTINCONGRUENCE	<- paste(HIVC.CODE.HOME,"pkg/misc/hivclu.startme.R -exeRECOMB.PLOTINCONGRUENCE",sep='/')

#' @export
PR.EXAML.PARSER	<- system.file(package="hivclust", "ext", "ExaML-parser") 

#' @export
PR.EXAML.STARTTREE	<- system.file(package="hivclust", "ext", "ExaML-parsimonator")

#' @export
PR.EXAML.EXAML	<- system.file(package="hivclust", "ext", "examl")

#' @export
PR.EXAML.BS		<- system.file(package="hivclust", "ext", "ExaML-raxml")

#' @export
PR.BEAST		<- {tmp<- c("/Applications/BEAST_1.7.5/bin/beast","beast"); names(tmp)<- c("debug","cx1.hpc.ic.ac.uk"); tmp } 

#' @export
PR.BEASTMCC		<- {tmp<- c("/Applications/BEAST_1.7.5/bin/treeannotator","treeannotator"); names(tmp)<- c("debug","cx1.hpc.ic.ac.uk"); tmp }

#' @export
PR.BEASTEVALRUN	<- paste(HIVC.CODE.HOME,"pkg/misc/hivclu.startme.R -exeBEASTEVALRUN",sep='/')

#' @export
PR.BEASTPOOLRUN	<- paste(HIVC.CODE.HOME,"pkg/misc/hivclu.startme.R -exeBEASTPOOLRUN",sep='/')

#' @export
PR.BEAST2		<- system.file(package="hivclust", "ext", "beast2.jar") 

#' @export
PR.BEAST2SA		<- system.file(package="hivclust", "ext", "beast2-SA.jar")

#' @export
HPC.NPROC		<- {tmp<- c(1,4); names(tmp)<- c("debug","cx1.hpc.ic.ac.uk"); tmp}

#' @export
HPC.MPIRUN		<- {tmp<- c("mpirun","mpiexec"); names(tmp)<- c("debug","cx1.hpc.ic.ac.uk"); tmp}

#' @export
HPC.CX1.IMPERIAL<- "cx1.hpc.ic.ac.uk"		#this is set to system('domainname',intern=T) for the hpc cluster of choice

#' @export
HPC.MEM			<- "1750mb"

#' @export
HPC.LOAD		<- "module load intel-suite mpi R/2.15 raxml examl/2013-05-09 beast/1.7.4"


#generate clustalo command
#' @export
hivc.cmd.blast.makedb<- function(indir, infile, signat=paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''), outdir=indir, with.mask=1, prog.mask= PR.BLASTMASK, prog.makedb=PR.BLASTMAKEDB, nproc=1, verbose=1)
{
	if(verbose) cat(paste("\nprocess",infile,"\n"))
	file	<- paste(indir, infile,sep='/')
	#verbose stuff
	cmd<- "#######################################################
# run blast.makedb
#######################################################"
	cmd<- paste(cmd,paste("\necho \'run ",prog.makedb,"\'",sep=''))
	if(with.mask)
	{
		cmd<- paste(cmd,'\n',prog.mask," -in ",paste(indir,'/',infile,'_',signat,".fasta",sep='')," -infmt fasta -mk_counts  -parse_seqids -out ",paste(outdir,'/',infile,'_',signat,".count",sep=''), sep='')
		cmd<- paste(cmd,'\n',prog.mask,  " -in ",paste(indir,'/',infile,'_',signat,".fasta",sep='')," -infmt fasta -ustat ",paste(outdir,'/',infile,'_',signat,".count",sep='')," -outfmt maskinfo_asn1_bin -parse_seqids -out ",paste(outdir,'/',infile,'_',signat,".asnb",sep=''), sep='')
		cmd<- paste(cmd,'\n',prog.makedb," -in ",paste(indir,'/',infile,'_',signat,".fasta",sep='')," -input_type fasta -dbtype nucl -parse_seqids -mask_data ",paste(outdir,'/',infile,'_',signat,".asnb",sep='')," -out ",paste(outdir,'/',infile,'_',signat,".blastdb",sep=''),' -title \"',infile,'\"',sep='')
	}
	else
		cmd<- paste(cmd,'\n',prog.makedb," -in ",paste(indir,'/',infile,'_',signat,".fasta",sep='')," -input_type fasta -dbtype nucl -parse_seqids -out ",paste(outdir,'/',infile,'_',signat,".blastdb",sep=''),' -title \"',infile,'\"',sep='')
	#clean up
	
	#verbose stuff
	cmd<- paste(cmd,paste("\necho \'end ",prog.makedb,"\'\n\n",sep=''))
	cmd
}

#' default values taken from http://indra.mullins.microbiol.washington.edu/viroblast/viroblast.php
#' @export
hivc.cmd.blast<- function(indir, infile, insignat, dbdir, dbfile, dbsignat, outdir=indir, outfile=infile, outsignat=insignat, prog.blastn= PR.BLASTN, blast.task="blastn", blast.max_target_seqs=10, blast.evalue=10, blast.wordsize=11, blast.gapopen=5, blast.gapextend=2, blast.penalty=-3, blast.reward= 2, blast.dust= "no", nproc=1, verbose=1)
{
	if(verbose) cat(paste("\nprocess",infile,"\n"))
	file	<- paste(indir, infile,sep='/')
	#verbose stuff
	cmd<- "#######################################################
# run blastn
#######################################################"
	cmd<- paste(cmd,paste("\necho \'run ",prog.blastn,"\'",sep=''))	
	cmd<- paste(cmd,'\n',prog.blastn," -query ",paste(indir,'/',infile,'_',insignat,".fasta",sep='')," -db ",paste(dbdir,'/',dbfile,'_',dbsignat,".blastdb",sep='')," -out ",paste(outdir,'/',outfile,'_',outsignat,".blast",sep=''),sep='')
	cmd<- paste(cmd," -task ",blast.task," -max_target_seqs ",blast.max_target_seqs," -evalue ",blast.evalue," -word_size ",blast.wordsize," -gapopen ",blast.gapopen," -gapextend ",blast.gapextend," -penalty ",blast.penalty," -reward ",blast.reward," -dust ",blast.dust," -outfmt 6",sep='')	
	#verbose stuff
	cmd<- paste(cmd,paste("\necho \'end ",prog.blastn,"\'\n\n",sep=''))
	cmd
}

#generate 3seq command
#' @export
hivc.cmd.recombination.run.3seq<- function(infile, outfile=paste(infile,".3s.rec",sep=''), recomb.3seq.siglevel=0.1, recomb.3seq.testvsall.beginatseq=NA, recomb.3seq.testvsall.endatseq=NA, prog= PR.RECOMB.3SEQ, nproc=1, verbose=1)
{
	cmd<- "#######################################################
# start: run 3Seq
#######################################################"
	cmd<- paste(cmd,paste("\necho \'run ",prog,"\'\n",sep=''))				
	#default commands
	cmd<- paste(cmd,prog,sep=" ")	
	cmd<- paste(cmd, " -x ",infile," -id ",outfile, sep='' )	
	if(!is.na(recomb.3seq.testvsall.beginatseq))
		cmd<- paste(cmd, " -b",recomb.3seq.testvsall.beginatseq, sep='' )
	if(!is.na(recomb.3seq.testvsall.endatseq))
		cmd<- paste(cmd, " -e",recomb.3seq.testvsall.endatseq, sep='' )	
	cmd<- paste(cmd, " -hs -t",recomb.3seq.siglevel, sep='')					
	#verbose stuff
	cmd<- paste(cmd,paste("\necho \'end ",prog,"\'\n",sep=''))
	cmd<- paste(cmd,"#######################################################
# end: run 3Seq
#######################################################\n",sep='')
	cmd			
}

#process 3seq output
#' @export
hivc.cmd.recombination.process.3SEQ.output<- function(indir, infile, insignat, prog= PR.RECOMB.PROCESS3SEQOUTPUT, resume=1, verbose=1)
{
	cmd<- "#######################################################
# start: run hivc.prog.recombination.process.3SEQ.output
#######################################################"
	cmd<- paste(cmd,paste("\necho \'run ",prog,"\'\n",sep=''))
	#default commands
	cmd<- paste(cmd,prog," -v=",verbose," -resume=",resume,sep='')
	cmd<- paste(cmd," -indir=",indir," -infile=",infile," -insignat=",insignat,sep='')
	#verbose stuff
	cmd<- paste(cmd,paste("\necho \'end ",prog,"\'\n",sep=''))
	cmd<- paste(cmd,"#######################################################
# end: run hivc.prog.recombination.process.3SEQ.output
#######################################################\n",sep='')
	cmd	
}

#check 3seq candidate recombinants
#' @export
hivc.cmd.recombination.check.candidates<- function(indir, infile, insignat, triplet.id, prog= PR.RECOMB.CHECKCANDIDATES, resume=1, verbose=1)
{
	cmd<- "#######################################################
# start: hivc.prog.recombination.check.candidates
#######################################################"
	cmd<- paste(cmd,paste("\necho \'run ",prog,"\'\n",sep=''))
	#default commands
	cmd<- paste(cmd,prog," -v=",verbose," -resume=",resume,sep='')
	cmd<- paste(cmd," -indir=",indir," -infile=",infile," -insignat=",insignat," -tripletid=",triplet.id,sep='')
	#verbose stuff
	cmd<- paste(cmd,paste("\necho \'end ",prog,"\'\n",sep=''))
	cmd<- paste(cmd,"#######################################################
# end: hivc.prog.recombination.check.candidates
#######################################################\n",sep='')
	cmd	
}

hivc.cmd.recombination.plot.incongruence<- function(indir, infile, insignat, triplet.id=NA, prog= PR.RECOMB.PLOTINCONGRUENCE, opt.select=NA,verbose=1)
{
	cmd<- "#######################################################
# start: hivc.prog.recombination.plot.incongruence
#######################################################"
	cmd<- paste(cmd,paste("\necho \'run ",prog,"\'\n",sep=''))
	#default commands
	cmd<- paste(cmd,prog," -v=",verbose,sep='')
	cmd<- paste(cmd," -indir=",indir," -infile=",infile," -insignat=",insignat, sep='')
	if(!is.na(opt.select))
		cmd<- paste(cmd," -select=",opt.select,sep='')
	if(!is.na(triplet.id))
		cmd<- paste(cmd," -tripletid=",triplet.id,sep='')
	#verbose stuff
	cmd<- paste(cmd,paste("\necho \'end ",prog,"\'\n",sep=''))
	cmd<- paste(cmd,"#######################################################
# end: hivc.prog.recombination.plot.incongruence
#######################################################\n",sep='')
	cmd	
}

#generate clustalo command
#' @export
hivc.cmd.clustalo<- function(indir, infiles, signat=paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''), outdir=indir, prog= PR.CLUSTALO, hmm=PR.CLUSTALO.HMM, nproc=1, verbose=1)
{
	ans<- lapply(infiles,function(x)
			{
				if(verbose) cat(paste("\nprocess",x,"\n"))				
				#verbose stuff
				cmd<- "#######################################################
# run clustalo
#######################################################"
				cmd<- paste(cmd,paste("\necho \'run ",PR.CLUSTALO,"\'\n",sep=''))
				
				#default commands
				cmd<- paste(cmd,PR.CLUSTALO,sep=" ")								
				cmd<- paste(cmd, " --infmt=fa --outfmt=fa --force", sep='')
				cmd<- paste(cmd, " --threads=",nproc," ", sep='' )
				#cmd<- paste(cmd, " --hmm-in=",hmm," ", sep='' )		not supported in clustalo v1.1.0
				
				#file in/out
				tmp<- paste(indir,x,sep='/')
				cmd<- paste(cmd, paste("--in",tmp,sep=' ') )				
				tmp<- paste(outdir,paste(x,ifelse(nchar(signat),'.',''),signat,".clustalo",sep=''),sep='/')
				cmd<- paste(cmd, paste("--out",tmp,sep=' ') )
				
				#verbose stuff
				cmd<- paste(cmd,paste("\necho \'end ",PR.CLUSTALO,"\'\n\n",sep=''))
				cmd
			})
	if(length(ans)==1)
		ans<- unlist(ans)
	ans
}

hivc.cmd.beast.poolrunxml<- function(indir, infile, insignat, indircov, infilecov, infiletree, infilexml, outsignat, pool.ntip, infilexml.opt="txs4clu", infilexml.template="standard", opt.brl="dist.brl.casc", thresh.brl=0.096, thresh.bs=0.8, prog= PR.BEASTPOOLRUN, resume=1, verbose=1)
{
	cmd<- "#######################################################
# start: run BEAST for pooled clusters capturing MSM transmission
#######################################################"
	cmd<- paste(cmd,paste("\necho \'run ",prog,"\'\n",sep=''))
	#default commands
	cmd<- paste(cmd,prog," -v=",verbose," -resume=",resume,sep='')
	cmd<- paste(cmd," -indir=",indir," -infile=",infile," -insignat=",insignat," -indircov=",indircov," -infilecov=",infilecov,sep='')
	cmd<- paste(cmd," -infilexml=",infilexml," -infilexml.opt=",infilexml.opt," -infilexml.template=",infilexml.template," -infiletree=",infiletree," -outdir=",indir," -outsignat=",outsignat,sep='')
	cmd<- paste(cmd," -pool.ntip=",pool.ntip," -thresh.brl=",thresh.brl," -thresh.bs=",thresh.bs," -opt.brl=",opt.brl,sep='')
	#verbose stuff
	cmd<- paste(cmd,paste("\necho \'end ",prog,"\'\n\n",sep=''))
	cmd	
}

hivc.cmd.beast.evalrun<- function(indir, infile, insignat, infilexml.opt, infilexml.template, pool.n, prog= PR.BEASTEVALRUN, verbose=1)
{
	cmd<- "#######################################################
# start: evaluate BEAST MCC trees
#######################################################"
	cmd<- paste(cmd,paste("\necho \'run ",prog,"\'\n",sep=''))
	#default commands
	cmd<- paste(cmd,prog," -v=",verbose," -indir=",indir," -infile=",infile," -insignat=",insignat,sep='')
	cmd<- paste(cmd," -infilexml.opt=",infilexml.opt," -infilexml.template=",infilexml.template," -pool.n=",pool.n,sep='')
	#verbose stuff
	cmd<- paste(cmd,paste("\necho \'end ",prog,"\'\n",sep=''))
	cmd<- paste(cmd,"#######################################################
# end: evaluate BEAST MCC trees
#######################################################\n",sep='')
	cmd	
}

hivc.cmd.clustering.msm<- function(indir, infile, insignat, indircov, infilecov, opt.brl, thresh.brl, thresh.bs, prog= PR.CLUSTMSM, resume=1, verbose=1)
{
	cmd<- "#######################################################
# extract clusters capturing MSM transmission
#######################################################"
	cmd<- paste(cmd,paste("\necho \'run ",prog,"\'\n",sep=''))
	#default commands
	cmd<- paste(cmd,prog," -v=",verbose," -resume=",resume,sep='')
	cmd<- paste(cmd," -indir=",indir," -infile=",infile," -insignat=",insignat," -indircov=",indircov," -infilecov=",infilecov," -thresh.brl=",thresh.brl," -thresh.bs=",thresh.bs," -opt.brl=",opt.brl,sep='')
	#verbose stuff
	cmd<- paste(cmd,paste("\necho \'end ",prog,"\'\n\n",sep=''))
	cmd
}

hivc.cmd.clustering.tptn<- function(indir, infile, insignat, indircov, infilecov, opt.brl="dist.brl.casc", patient.n=15700, prog= PR.CLUSTTPTN, resume=1, verbose=1)
{
	cmd<- "#######################################################
# compute clustering true pos and true neg statistics
#######################################################"
	cmd<- paste(cmd,paste("\necho \'run ",prog,"\'\n",sep=''))
	#default commands
	cmd<- paste(cmd,prog," -v=",verbose," -resume=",resume,sep='')
	cmd<- paste(cmd," -indir=",indir," -infile=",infile," -insignat=",insignat," -indircov=",indircov," -infilecov=",infilecov," -patient.n=",patient.n," -opt.brl=",opt.brl,sep='')
	#verbose stuff
	cmd<- paste(cmd,paste("\necho \'end ",prog,"\'\n\n",sep=''))
	cmd
}

hivc.cmd.clustering<- function(indir, infile, insignat, opt.brl, thresh.brl, thresh.bs, prog= PR.CLUST, resume=1, verbose=1)
{
	cmd<- "#######################################################
# compute clustering 
#######################################################"
	cmd<- paste(cmd,paste("\necho \'run ",prog,"\'\n",sep=''))
	#default commands
	cmd<- paste(cmd,prog," -v=",verbose," -resume=",resume,sep='')
	cmd<- paste(cmd," -indir=",indir," -infile=",infile," -insignat=",insignat," -thresh.brl=",thresh.brl," -thresh.bs=",thresh.bs," -opt.brl=",opt.brl,sep='')
	#verbose stuff
	cmd<- paste(cmd,paste("\necho \'end ",prog,"\'\n\n",sep=''))
	cmd
}

#' @export
hivc.cmd.preclustering<- function(indir, infile, insignat, indircov, infilecov, prog= PR.PRECLUST, resume=0, verbose=1)
{
	cmd<- "#######################################################
# precompute clustering helper objects that are memory intensive
#######################################################"
	cmd<- paste(cmd,paste("\necho \'run ",prog,"\'\n",sep=''))
	#default commands
	cmd<- paste(cmd,prog," -v=",verbose," -resume=",resume,sep='')
	cmd<- paste(cmd," -indir=",indir," -infile=",infile," -insignat=",insignat," -indircov=",indircov," -infilecov=",infilecov,sep='')
	#verbose stuff
	cmd<- paste(cmd,paste("\necho \'end ",prog,"\'\n\n",sep=''))
	cmd
}

#' @export
hivc.cmd.get.geneticdist<- function(indir, infile, signat, gd.max, outdir=indir, prog= PR.GENDISTMAT, resume=1, verbose=1)
{
	cmd<- "#######################################################
# start: compute genetic distance matrix of sequence alignment
#######################################################"
	cmd<- paste(cmd,paste("\necho \'run ",prog,"\'\n",sep=''))
	#default commands
	cmd<- paste(cmd,prog," -v=",verbose," -resume=",resume,sep='')
	cmd<- paste(cmd," -indir=",indir," -infile=",infile," -outdir=",outdir," -signat=",signat," -maxgd=",gd.max,sep='')
	#verbose stuff
	cmd<- paste(cmd,paste("\necho \'end ",prog,"\'\n",sep=''))
	cmd<- paste(cmd,"#######################################################
# end: compute genetic distance matrix of sequence alignment
#######################################################\n",sep='')
	cmd
}

#' @export
hivc.cmd.get.firstseq<- function(indir, infile, signat.in, signat.out, outdir=indir, prog= PR.FIRSTSEQ, resume=1, verbose=1)
{
	cmd<- "#######################################################
# extract firstseq
#######################################################"
	cmd<- paste(cmd,paste("\necho \'run ",prog,"\'\n",sep=''))
	#default commands
	cmd<- paste(cmd,prog," -v=",verbose," -resume=",resume,sep='')
	cmd<- paste(cmd," -indir=",indir," -infile=",infile," -outdir=",outdir," -insignat=",signat.in," -outsignat=",signat.out,sep='')
	#verbose stuff
	cmd<- paste(cmd,paste("\necho \'end ",prog,"\'\n\n",sep=''))
	cmd
}

#' @export
hivc.cmd.examl<- function(indir, infile, signat.in, signat.out, outdir=indir, prog.parser= PR.EXAML.PARSER, prog.starttree= PR.EXAML.STARTTREE, args.starttree.seed=12345, args.starttree.bsid= NA, prog.examl= PR.EXAML.EXAML, args.examl="-m GAMMA -D", resume=1, verbose=1)
{
	if(is.na(args.starttree.bsid))
		args.starttree.bsid	<- "000"
	else
		args.starttree.bsid	<-	sprintf("%03d",args.starttree.bsid)
	cmd<- "#######################################################
# start: compute ExaML tree
#######################################################"
	cmd<- paste(cmd,paste("\necho \'run ",prog.parser,"\'\n",sep=''))
	#if output files are found and resume, don t do anything
	if(resume)
	{
		cmd		<- paste(cmd,"[ -s ",outdir,'/ExaML_result.',infile,'_',signat.in,".finaltree.",args.starttree.bsid," ] && ", sep='')	#if file non-zero
		cmd		<- paste(cmd,"[ -s ",outdir,'/ExaML_info.',infile,'_',signat.in,".finaltree.",args.starttree.bsid," ] ", sep='')		#and if file non-zero
		cmd		<- paste(cmd,"&& exit 1\n",sep='')		#then exit
	}
	#default commands for parser					
	cmd			<- paste(cmd,"CWDEXAML=$(pwd)\n",sep='')
	cmd			<- paste(cmd,"cd ",outdir,'\n',sep='')
	tmp			<- paste(indir,paste(infile,'_',signat.in,".phylip.",args.starttree.bsid,sep=''),sep='/')
	cmd			<- paste(cmd,prog.parser," -m DNA -s ",tmp,sep='')
	tmp			<- paste(infile,'_',signat.out,".phylip.examl.",args.starttree.bsid,sep='')
	cmd			<- paste(cmd," -n ",tmp,sep='')
	#verbose stuff for parser	
	cmd			<- paste(cmd,paste("\necho \'end ",prog.parser,"\'",sep=''))
		
	cmd			<- paste(cmd,paste("\necho \'run ",prog.starttree,"\'\n",sep=''))
	#default commands for start tree
	tmp			<- paste(indir,paste(infile,'_',signat.in,".phylip.",args.starttree.bsid,sep=''),sep='/')	
	cmd			<- paste(cmd,prog.starttree," -p",args.starttree.seed," -s ",tmp,sep='')	
	tmp			<- paste(infile,'_',signat.out,".starttree.",args.starttree.bsid,sep='')
	cmd			<- paste(cmd," -n ",tmp,sep='')
	#verbose stuff
	cmd			<- paste(cmd,paste("\necho \'end ",prog.starttree,"\'",sep=''))
	
	cmd			<- paste(cmd,paste("\necho \'run ",prog.examl,"\'\n",sep=''))
	#default commands for final tree
	tmp			<- hivc.get.hpcsys()
	if(tmp=="debug")
		cmd		<- paste(cmd,HPC.MPIRUN[tmp]," -np ",HPC.NPROC[tmp],' ',prog.examl,' ',args.examl,sep='')
	else if(tmp==HPC.CX1.IMPERIAL)
		cmd		<- paste(cmd,HPC.MPIRUN[tmp],' ',prog.examl,' ',args.examl,sep='')
	else	
		stop("unknown hpc sys")
	tmp			<- paste(infile,'_',signat.out,".phylip.examl.",args.starttree.bsid,".binary",sep='')
	cmd			<- paste(cmd," -s ",tmp,sep='')
	tmp			<- paste("RAxML_parsimonyTree.",infile,'_',signat.out,".starttree.",args.starttree.bsid,".0",sep='')
	cmd			<- paste(cmd," -t ",tmp,sep='')
	tmp			<- paste(infile,'_',signat.out,".finaltree.",args.starttree.bsid,sep='')
	cmd			<- paste(cmd," -n ",tmp,sep='')		
	cmd			<- paste(cmd,paste("\necho \'end ",prog.examl,"\'",sep=''))
	
	#delete ExaML output that is not further needed 
	cmd			<- paste(cmd,paste("\necho \'start cleanup\'",sep=''))	
	cmd			<- paste(cmd,"\nfind . -name \'*phylip.",args.starttree.bsid,"*\' -delete",sep='')
	cmd			<- paste(cmd,"\nfind . -name \'*phylip.examl.",args.starttree.bsid,"*\' -delete",sep='')
	cmd			<- paste(cmd,"\nfind . -name \'*starttree.",args.starttree.bsid,"*\' -delete",sep='')
	cmd 		<- paste(cmd,"\nfind . -name \'ExaML_binaryCheckpoint.*?finaltree.",args.starttree.bsid,"*\' -delete", sep='' )
	cmd 		<- paste(cmd,"\nfind . -name \'ExaML_log.*?finaltree.",args.starttree.bsid,"*\' -delete", sep='' )
	cmd			<- paste(cmd,paste("\necho \'end cleanup\'",sep=''))		
	cmd			<- paste(cmd,"\ncd $CWDEXAML",sep='')
	cmd			<- paste(cmd,"\n#######################################################
# end: compute ExaML tree
#######################################################\n",sep='')
	cmd
}

#' @export
#' 	creates a shell command to create a new bootstrap alignment over codon positions of an input alignment 
hivc.cmd.examl.bsalignment<- function(indir, infile, signat.in, signat.out, bs.id, outdir=indir, prog.bscreate= PR.EXAML.BSCREATE, opt.bootstrap.by="codon",resume=0, verbose=1)
{
	cmd			<- paste("#######################################################
# start: create and check bootstrap alignment
#######################################################\n",sep='')
	cmd			<- paste(cmd,prog.bscreate," -resume=",resume," -bootstrap=",bs.id," -by=",opt.bootstrap.by,sep='')
	cmd			<- paste(cmd," -indir=",indir," -infile=",infile," -outdir=",outdir,sep='')
	cmd			<- paste(cmd," -insignat=",signat.in," -outsignat=",signat.out,sep='')
	cmd			<- paste(cmd,"\n#######################################################
# end: create and check bootstrap alignment
#######################################################",sep='')			
	cmd
}

#' @export
hivc.cmd.examl.bootstrap.on.one.machine<- function(indir, infile, signat.in, signat.out, bs.from=0, bs.to=99, bs.n=bs.to-bs.from+ifelse(bs.from==0,1,0), outdir=indir, prog.parser= PR.EXAML.PARSER, prog.starttree= PR.EXAML.STARTTREE, prog.examl=PR.EXAML.EXAML, opt.bootstrap.by="codon", args.examl="-m GAMMA -D", prog.supportadder=PR.EXAML.BS, tmpdir.prefix="examl", resume=1, verbose=1)
{
	hpcsys			<- hivc.get.hpcsys()
	hpcsys			<- "cx1.hpc.ic.ac.uk"
	#create number of seeds for the number of runs being processed, which could be less than bs.n
	bs.id			<- seq.int(bs.from,bs.to)
	bs.seeds		<- floor( runif(length(bs.id), 1e4, 1e5-1) )
	tmpdir			<- paste("$CWD/",tmpdir.prefix,'_',format(Sys.time(),"%y-%m-%d-%H-%M-%S"),sep='')
	cmd				<- sapply(seq_along(bs.seeds), function(i)
			{
				cmd			<- ''
				if(hpcsys=="debug")						#my MAC - don t use scratch
				{
					cmd		<- paste(cmd,hivc.cmd.examl.bsalignment(indir, infile, signat.in, signat.out, bs.id[i], opt.bootstrap.by=opt.bootstrap.by, outdir=indir, verbose=verbose),sep='\n')
					cmd		<- paste(cmd,hivc.cmd.examl(indir, infile, signat.in, signat.out, outdir=outdir, prog.parser= prog.parser, prog.starttree= prog.starttree, args.starttree.seed=bs.seeds[i], args.starttree.bsid= bs.id[i], prog.examl=prog.examl, args.examl=args.examl, resume=resume, verbose=verbose),sep='\n')
				}
				else if(hpcsys=="cx1.hpc.ic.ac.uk")		#imperial - use scratch directory
				{										
					if(i==1)
					{
						cmd	<- paste(cmd,"\nCWD=$(pwd)",sep='')
						cmd	<- paste(cmd,"\necho $CWD",sep='')
						cmd	<- paste(cmd,"\nmkdir -p ",tmpdir,sep='')
						tmp	<- paste(indir,'/',infile,'_',gsub('/',':',signat.in),".R",sep='')
						cmd	<- paste(cmd,"\ncp ",tmp," ",tmpdir,sep='')
					}
					cmd		<- paste(cmd,hivc.cmd.examl.bsalignment(tmpdir, infile, signat.in, signat.out, bs.id[i], opt.bootstrap.by=opt.bootstrap.by, outdir=tmpdir, verbose=verbose),sep='\n')
					cmd		<- paste(cmd,hivc.cmd.examl(tmpdir, infile, signat.in, signat.out, outdir=tmpdir, prog.parser= prog.parser, prog.starttree= prog.starttree, args.starttree.seed=bs.seeds[i], args.starttree.bsid= bs.id[i], prog.examl=prog.examl, args.examl=args.examl, resume=resume, verbose=verbose),sep='\n')
				}
				cmd
			})
	cmd				<- paste(cmd, collapse='')
	cmd				<- paste(cmd,"#######################################################
# start: check if all ExaML boostrap trees have been computed and if yes create ExaML bootstrap tree
#######################################################",sep='')	
	cmd			<- paste(cmd,"\nCWD=$(pwd)\n",sep='')
	cmd			<- paste(cmd,"cd ",tmpdir,sep='') 
	cmd			<- paste(cmd,paste("\necho \'check if all bootstrap samples have been computed\'",sep=''))
	#compute bs tree even when some errors	
	tmp			<- paste("\nif [ $(find . -name 'ExaML_result.",infile,'_',gsub('/',':',signat.in),".finaltree*' | wc -l) -gt ",round(bs.n*0.95)," ]; then",sep='')
	cmd			<- paste(cmd,tmp,sep='')
	cmd			<- paste(cmd,paste("\n\techo \'all bootstrap samples computed -- find best tree and add bootstrap support values\'",sep=''))				
	tmp			<- c(	paste("ExaML_result.",infile,'_',signat.out,".finaltree",sep=''), 	paste("ExaML_result.",infile,'_',signat.out,".bstree",sep=''))
	#add all bootstrap final trees into bs tree file
	cmd			<- paste(cmd,"\n\tfor i in $(seq 0 ",bs.n-1,"); do cat ",tmp[1],".$(printf %03d $i) >> ",tmp[2],"; done",sep='')
	#identify suffix finaltree.XXX of final tree with highest likelihood
	cmd			<- paste(cmd,"\n\tBSBEST=$(grep 'Likelihood of best tree' ExaML_info.* | awk '{print $5,$1;}' | sort -n | tail -1 | grep -o 'finaltree.*' | cut -d':' -f 1)",sep='')
	cmd			<- paste(cmd,paste("\n\techo \"best tree is $BSBEST\"",sep=''))		
	#create final tree with bootstrap support values
	tmp			<- c(	paste(infile,'_',signat.out,".phylip.examl.binary",sep=''),	paste("ExaML_result.",infile,'_',signat.out,".$BSBEST",sep=''),	paste("ExaML_result.",infile,'_',signat.out,".bstree",sep=''), paste(infile,'_',signat.out,".supporttree",sep=''))
	cmd			<- paste(cmd,"\n\t",prog.supportadder," -f b -m GTRCAT -s ",tmp[1]," -t ",tmp[2]," -z ",tmp[3]," -n ",tmp[4],sep='' )
	cmd			<- paste(cmd,paste("\n\techo \'all bootstrap samples computed -- found best tree and added bootstrap support values\'",sep=''))										
	#delete ExaML output that is not further needed
	cmd			<- paste(cmd,paste("\n\techo \'start cleanup\'",sep=''))
	cmd			<- paste(cmd,"\n\trm RAxML_info*",sep=' ')								
	cmd			<- paste(cmd,paste("\n\tmv RAxML_bipartitions.",infile,'_',signat.out,".supporttree ",infile,"_examlbs",bs.n,'_',signat.out,".newick",sep=''),sep='')				
	cmd			<- paste(cmd,paste("\n\trm RAxML_bipartitionsBranchLabels.",infile,'_',signat.out,".supporttree",sep=''),sep='')
	cmd			<- paste(cmd,paste("\n\trm ExaML_result.",infile,'_',signat.out,".bstree",sep=''),sep='')									
	cmd			<- paste(cmd,paste("\n\techo \'end cleanup\'",sep=''))
	#copy bstree to outdir
	cmd			<- paste(cmd,paste("\n\tcp -f ",infile,"_examlbs",bs.n,'_',signat.out,".newick",' ',outdir,sep=''),sep='')
	cmd			<- paste(cmd,"\nfi",sep='')
	cmd			<- paste(cmd,"\n#######################################################
# end: check if all ExaML boostrap trees have been computed and if yes create ExaML bootstrap tree
#######################################################",sep='')
cmd			<- paste(cmd,"\n#######################################################
# start: zip and copy ExaML output
#######################################################",sep='')	
	#outside if:	zip and copy ExaML output to outdir just in case something went wrong
	cmd			<- paste(cmd,paste("\nzip ",infile,'_examlout_',signat.out,".zip  ExaML_result.",infile,'_',signat.out,".* ExaML_info.",infile,'_',signat.out,".*",sep=''),sep='')
	cmd			<- paste(cmd,paste("\ncp -f ",infile,'_examlout_',signat.out,".zip",' ',outdir,sep=''),sep='')
	cmd			<- paste(cmd,paste("\nrm ExaML_result.",infile,'_',signat.out,".* ExaML_info.",infile,'_',signat.out,".*",sep=''),sep='')
	cmd			<- paste(cmd,"\ncd $CWD",sep='')
	cmd			<- paste(cmd,"\n#######################################################
# end: zip and copy ExaML output
#######################################################\n",sep='')				
	cmd
}

#' @export
hivc.cmd.examl.bootstrap<- function(indir, infile, signat.in, signat.out, bs.from=0, bs.to=99, bs.n=bs.to-bs.from+ifelse(bs.from==0,1,0), outdir=indir, prog.parser= PR.EXAML.PARSER, prog.starttree= PR.EXAML.STARTTREE, prog.examl=PR.EXAML.EXAML, opt.bootstrap.by="codon", args.examl="-m GAMMA -D", prog.supportadder=PR.EXAML.BS, tmpdir.prefix="examl", resume=1, verbose=1)
{
	hpcsys			<- hivc.get.hpcsys()
	#hpcsys			<- "cx1.hpc.ic.ac.uk"
	#create number of seeds for the number of runs being processed, which could be less than bs.n
	bs.id			<- seq.int(bs.from,bs.to)
	bs.seeds		<- floor( runif(length(bs.id), 1e4, 1e5-1) )
	tmpdir.prefix	<- paste(tmpdir.prefix,'_',format(Sys.time(),"%y-%m-%d-%H-%M-%S"),sep='')
	lapply(seq_along(bs.seeds), function(i)
			{
				cmd			<- ''				
				if(hpcsys=="debug")						#my MAC - don t use scratch
				{
					cmd		<- paste(cmd,hivc.cmd.examl.bsalignment(indir, infile, signat.in, signat.out, bs.id[i], opt.bootstrap.by=opt.bootstrap.by, outdir=indir, verbose=verbose),sep='\n')
					cmd		<- paste(cmd,hivc.cmd.examl(indir, infile, signat.in, signat.out, outdir=outdir, prog.parser= prog.parser, prog.starttree= prog.starttree, args.starttree.seed=bs.seeds[i], args.starttree.bsid= bs.id[i], prog.examl=prog.examl, args.examl=args.examl, resume=resume, verbose=verbose),sep='\n')
				}
				else if(hpcsys=="cx1.hpc.ic.ac.uk")		#imperial - use scratch directory
				{										
					if(resume)
					{
						cmd	<- paste(cmd,"\nif ", sep='')
						cmd	<- paste(cmd,"[ ! -s ",outdir,'/ExaML_result.',infile,'_',signat.in,".finaltree.",sprintf("%03d",bs.id[i])," ]", sep='')
						cmd	<- paste(cmd," || ", sep='')
						cmd	<- paste(cmd,"[ ! -s ",outdir,'/ExaML_info.',infile,'_',signat.in,".finaltree.",sprintf("%03d",bs.id[i])," ];", sep='')
						cmd	<- paste(cmd," then\n",sep='')
						cmd	<- paste(cmd,"#######################################################
# start: not indented if statement -- don t do anything if ExaML output exists already
#######################################################",sep='')
					}					
					cmd		<- paste(cmd,"CWD=$(pwd)\n",sep='\n')
					cmd		<- paste(cmd,"echo $CWD\n",sep='')
					tmpdir	<- paste("$CWD/",tmpdir.prefix,sep='')
					cmd		<- paste(cmd,"mkdir -p ",tmpdir,'\n',sep='')
					tmp		<- paste(indir,'/',infile,'_',gsub('/',':',signat.in),".R",sep='')
					cmd		<- paste(cmd,"cp ",tmp," ",tmpdir,sep='')
					cmd		<- paste(cmd,hivc.cmd.examl.bsalignment(tmpdir, infile, signat.in, signat.out, bs.id[i], opt.bootstrap.by=opt.bootstrap.by, outdir=tmpdir, verbose=verbose),sep='\n')
					cmd		<- paste(cmd,hivc.cmd.examl(tmpdir, infile, signat.in, signat.out, outdir=tmpdir, prog.parser= prog.parser, prog.starttree= prog.starttree, args.starttree.seed=bs.seeds[i], args.starttree.bsid= bs.id[i], prog.examl=prog.examl, args.examl=args.examl, resume=resume, verbose=verbose),sep='\n')
					cmd		<- paste(cmd,"cp -f ",tmpdir,"/ExaML_result.",infile,'_',gsub('/',':',signat.in),".finaltree.",sprintf("%03d",bs.id[i]),' ', outdir,'\n',sep='')
					cmd		<- paste(cmd,"cp -f ",tmpdir,"/ExaML_info.",infile,'_',gsub('/',':',signat.in),".finaltree.",sprintf("%03d",bs.id[i]),' ', outdir,'\n',sep='')
					if(resume)
					{
					cmd		<- paste(cmd,"#######################################################
# end: not indented if statement -- don t do anything if ExaML output exists already
#######################################################\nelse\n\techo 'resumed bootstrap number ",bs.id[i],"'\nfi\n",sep='')					
					}
				}				
				cmd			<- paste(cmd,"#######################################################
# start: check if all ExaML boostrap trees have been computed and if yes create ExaML bootstrap tree
#######################################################",sep='')	
				cmd			<- paste(cmd,"\nCWD=$(pwd)\n",sep='')
				cmd			<- paste(cmd,"cd ",outdir,sep='') 
				cmd			<- paste(cmd,paste("\necho \'check if all bootstrap samples have been computed\'",sep=''))
				tmp			<- paste("\nif [ $(find . -name 'ExaML_result.",infile,'_',gsub('/',':',signat.in),".finaltree*' | wc -l) -eq ",bs.n," ]; then",sep='')
				cmd			<- paste(cmd,tmp,sep='')
				cmd			<- paste(cmd,paste("\n\techo \'all bootstrap samples computed -- find best tree and add bootstrap support values\'",sep=''))				
				tmp			<- c(	paste("ExaML_result.",infile,'_',signat.out,".finaltree",sep=''), 	paste("ExaML_result.",infile,'_',signat.out,".bstree",sep=''))
				#add all bootstrap final trees into bs tree file
				cmd			<- paste(cmd,"\n\tfor i in $(seq 0 ",bs.n-1,"); do cat ",tmp[1],".$(printf %03d $i) >> ",tmp[2],"; done",sep='')
				#identify suffix finaltree.XXX of final tree with highest likelihood
				cmd			<- paste(cmd,"\n\tBSBEST=$(grep 'Likelihood of best tree' ExaML_info.* | awk '{print $5,$1;}' | sort -n | tail -1 | grep -o 'finaltree.*' | cut -d':' -f 1)",sep='')
				cmd			<- paste(cmd,paste("\n\techo \"best tree is $BSBEST\"",sep=''))		
				#create final tree with bootstrap support values
				tmp			<- c(	paste(infile,'_',signat.out,".phylip.examl.binary",sep=''),	paste("ExaML_result.",infile,'_',signat.out,".$BSBEST",sep=''),	paste("ExaML_result.",infile,'_',signat.out,".bstree",sep=''), paste(infile,'_',signat.out,".supporttree",sep=''))
				cmd			<- paste(cmd,"\n\t",prog.supportadder," -f b -m GTRCAT -s ",tmp[1]," -t ",tmp[2]," -z ",tmp[3]," -n ",tmp[4],sep='' )
				cmd			<- paste(cmd,paste("\n\techo \'all bootstrap samples computed -- found best tree and added bootstrap support values\'",sep=''))										
				#delete ExaML output that is not further needed
				cmd			<- paste(cmd,paste("\n\techo \'start cleanup\'",sep=''))
				cmd			<- paste(cmd,"\n\trm RAxML_info*",sep=' ')								
				cmd			<- paste(cmd,paste("\n\tmv RAxML_bipartitions.",infile,'_',signat.out,".supporttree ",infile,"_examlbs",bs.n,'_',signat.out,".newick",sep=''),sep='')				
				cmd			<- paste(cmd,paste("\n\trm RAxML_bipartitionsBranchLabels.",infile,'_',signat.out,".supporttree",sep=''),sep='')
				cmd			<- paste(cmd,paste("\n\trm ExaML_result.",infile,'_',signat.out,".bstree",sep=''),sep='')								
				cmd			<- paste(cmd,paste("\n\tzip ",infile,'_examlout_',signat.out,".zip  ExaML_result.",infile,'_',signat.out,".* ExaML_info.",infile,'_',signat.out,".*",sep=''),sep='')
				cmd			<- paste(cmd,paste("\n\trm ExaML_result.",infile,'_',signat.out,".* ExaML_info.",infile,'_',signat.out,".*",sep=''),sep='')
				cmd			<- paste(cmd,paste("\n\techo \'end cleanup\'",sep=''))
				cmd			<- paste(cmd,"\nfi",sep='')
				cmd			<- paste(cmd,"\ncd $CWD",sep='')
				cmd			<- paste(cmd,"\n#######################################################
# end: check if all ExaML boostrap trees have been computed and if yes create ExaML bootstrap tree
#######################################################\n",sep='')				
			})
}

#' @export
hivc.cmd.examl.bsstarttree<- function(indir, infile, signat.in, signat.out, bs.from=0, bs.to=99, bs.n=bs.to-bs.from+ifelse(bs.from==0,1,0),outdir=indir, prog.parser= PR.EXAML.PARSER, prog.starttree= PR.EXAML.STARTTREE, prog.examl=PR.EXAML.EXAML, args.examl="-m GAMMA -D", prog.supportadder=PR.EXAML.BS, resume=1, verbose=1)
{
	#create number of seeds for the number of runs being processed, which could be less than bs.n
	bs.id	<- seq.int(bs.from,bs.to)
	bs.seeds<- floor( runif(length(bs.id), 1e4, 1e5-1) )
	lapply(seq_along(bs.seeds), function(i)
			{
				cmd			<- paste("cp ",paste(indir,paste(infile,'_',signat.in,".phylip",sep=''),sep='/'),sep='')
				cmd			<- paste(cmd,' ',paste(indir,paste(infile,'_',signat.in,".phylip.",sprintf("%03d",bs.id[i]),sep=''),sep='/'),"\n",sep='')
				cmd			<- paste(cmd, hivc.cmd.examl(indir, infile, signat.in, signat.out, outdir=outdir, prog.parser= prog.parser, prog.starttree= prog.starttree, args.starttree.seed=bs.seeds[i], args.starttree.bsid= bs.id[i], prog.examl=prog.examl, args.examl=args.examl, resume=resume, verbose=verbose),sep='')
				curr.dir	<- getwd()
				cmd			<- paste(cmd,"#######################################################
# start: check if all ExaML boostrap trees have been computed and if yes create ExaML bootstrap tree
#######################################################",sep='')				
				cmd			<- paste(cmd,"\ncd ",outdir,sep='')
				cmd			<- paste(cmd,paste("\necho \'check if all bootstrap samples have been computed\'",sep=''))
				tmp			<- paste("\nif [ $(find . -name 'ExaML_result*' | wc -l) -eq ",bs.n," ]; then",sep='')
				cmd			<- paste(cmd,tmp,sep='')
				cmd			<- paste(cmd,paste("\n\techo \'all bootstrap samples computed -- find best tree and add bootstrap support values\'",sep=''))				
				tmp			<- c(	paste("ExaML_result.",infile,'_',signat.out,".finaltree",sep=''), 	paste("ExaML_result.",infile,'_',signat.out,".bstree",sep=''))
				#add all bootstrap final trees into bs tree file
				cmd			<- paste(cmd,"\n\tfor i in $(seq 0 ",bs.n-1,"); do cat ",tmp[1],".$(printf %03d $i) >> ",tmp[2],"; done",sep='')
				#identify suffix finaltree.XXX of final tree with highest likelihood
				cmd			<- paste(cmd,"\n\tBSBEST=$(grep 'Likelihood of best tree' ExaML_info.* | awk '{print $5,$1;}' | sort -n | tail -1 | grep -o 'finaltree.*' | cut -d':' -f 1)",sep='')
				cmd			<- paste(cmd,paste("\n\techo \"best tree is $BSBEST\"",sep=''))		
				#create final tree with bootstrap support values
				tmp			<- c(	paste(infile,'_',signat.out,".phylip.examl.binary",sep=''),	paste("ExaML_result.",infile,'_',signat.out,".$BSBEST",sep=''),	paste("ExaML_result.",infile,'_',signat.out,".bstree",sep=''), paste(infile,'_',signat.out,".supporttree",sep=''))
				cmd			<- paste(cmd,"\n\t",prog.supportadder," -f b -m GTRCAT -s ",tmp[1]," -t ",tmp[2]," -z ",tmp[3]," -n ",tmp[4],sep='' )
				cmd			<- paste(cmd,paste("\n\techo \'all bootstrap samples computed -- found best tree and added bootstrap support values\'",sep=''))										
				#delete ExaML output that is not further needed
				cmd			<- paste(cmd,paste("\n\techo \'start cleanup\'",sep=''))
				cmd			<- paste(cmd,"\n\trm RAxML_info*",sep=' ')
				cmd			<- paste(cmd,paste("\n\trm ",infile,'_',signat.out,".phylip.examl.binary",sep=''),sep='')				
				cmd			<- paste(cmd,paste("\n\tmv RAxML_bipartitions.",infile,'_',signat.out,".supporttree ",infile,"_examlbs",bs.n,'_',signat.out,".newick",sep=''),sep='')				
				cmd			<- paste(cmd,paste("\n\trm RAxML_bipartitionsBranchLabels.",infile,'_',signat.out,".supporttree",sep=''),sep='')
				cmd			<- paste(cmd,paste("\n\trm ExaML_result.",infile,'_',signat.out,".bstree",sep=''),sep='')								
				cmd			<- paste(cmd,paste("\n\ttar -zcf ",infile,'_examlout_',signat.out,".tar.gz  ExaML_result.",infile,'_',signat.out,".* ExaML_info.",infile,'_',signat.out,".*",sep=''),sep='')
				cmd			<- paste(cmd,paste("\n\trm ExaML_result.",infile,'_',signat.out,".* ExaML_info.",infile,'_',signat.out,".*",sep=''),sep='')
				cmd			<- paste(cmd,paste("\n\techo \'end cleanup\'",sep=''))
				cmd			<- paste(cmd,"\nfi",sep='')
				cmd			<- paste(cmd,"\ncd ",curr.dir,sep='')
				cmd			<- paste(cmd,"\n#######################################################
# end: check if all ExaML boostrap trees have been computed and if yes create ExaML bootstrap tree
#######################################################\n",sep='')
				#cat(cmd)
				#stop()
				#if [ $(find -E . -name 'ExaML_result*' | wc -l)==2 ]; then echo 'hello'; fi
			})
}

#' @export
hivc.cmd.examl.cleanup<- function(outdir, prog= PR.EXAML.EXAML)
{
	cmd<- "#######################################################
# clean up after ExaML tree
#######################################################"
	cmd<- paste(cmd,paste("\necho \'clean after ",prog,"\'\n",sep=''))	
	
	tmp<- list.files(outdir, full.names=1)
	#tmp<- tmp[c(grep("examlstarttree",tmp), grep("examlbin",tmp),grep("examl.binary",tmp),grep("phylip",tmp))]
	#cmd<- paste(cmd, "\nrm ", paste(tmp,collapse=' ',sep=''), '\n', sep='')	
	cmd<- paste(cmd,"rm ",paste(paste(outdir,c("RAxML_info*","RAxML_parsimonyTree*","ExaML_binaryCheckpoint*"),sep='/'),collapse=' ',sep=''),sep=' ')
	
	cmd<- paste(cmd,paste("\necho \'cleaned up after ",prog,"\'\n",sep=''))
	cmd
}

#' @export
hivc.cmd.beast.runxml<- function(indir, infile, insignat, prog.beast=PR.BEAST, prog.beastmcc=PR.BEASTMCC, beastmcc.burnin=500, beastmcc.heights="median", hpc.tmpdir.prefix="beast", hpc.ncpu=1)
{
	cmd		<- "#######################################################
# start: run BEAST
#######################################################"	
	hpcsys	<- hivc.get.hpcsys()
	#hpcsys<- "cx1.hpc.ic.ac.uk"
	cmd		<- paste(cmd,paste("\necho \'run ",prog.beast[hpcsys],"\'\n",sep=''))
	if(hpcsys=="debug")						#my MAC - don t use scratch
	{		
		tmp		<- paste(indir,'/',infile,'_',gsub('/',':',insignat),".xml",sep='')
		cmd		<- paste(cmd,prog.beast[hpcsys]," -strict -overwrite -working ",tmp,'\n',sep='')
	}
	else if(hpcsys=="cx1.hpc.ic.ac.uk")		#imperial - use scratch directory
	{
		cmd		<- paste(cmd,"CWD=$(pwd)\n",sep='')
		cmd		<- paste(cmd,"echo $CWD\n",sep='')		
		tmpdir	<- paste("$CWD/",hpc.tmpdir.prefix,'_',format(Sys.time(),"%y-%m-%d-%H-%M-%S"),sep='')
		cmd		<- paste(cmd,"mkdir -p ",tmpdir,'\n',sep='')
		tmp		<- paste(indir,'/',infile,'_',gsub('/',':',insignat),".xml",sep='')
		cmd		<- paste(cmd,"cp ",tmp," ",tmpdir,'\n',sep='')
		tmp		<- paste(tmpdir,'/',infile,'_',gsub('/',':',insignat),".xml",sep='')
		cmd		<- paste(cmd,prog.beast[hpcsys]," -strict -working -threads ",hpc.ncpu," ",tmp,'\n',sep='')	
		cmd		<- paste(cmd,"cp -f ",tmpdir,"/* ", indir,'\n',sep='')		
	}
	cmd		<- paste(cmd,"echo \'end ",prog.beast[hpcsys],"\'\n",sep='')
	cmd		<- paste(cmd,"#######################################################
# end: run BEAST
#######################################################\n",sep='')
	cmd<- paste(cmd,"#######################################################
# start: run TREEANNOTATOR
#######################################################\n",sep='')
	cmd		<- paste(cmd,"echo \'run ",prog.beastmcc[hpcsys],"\'\n",sep='')
	if(hpcsys=="debug")						#my MAC - don t use scratch
	{	
		tmp		<- paste(indir,'/',infile,'_',gsub('/',':',insignat),".timetrees",sep='')
		cmd		<- paste(cmd,prog.beastmcc[hpcsys]," -burnin ",beastmcc.burnin," -heights ",beastmcc.heights," ",tmp,sep='')
		tmp		<- paste(indir,'/',infile,"_mcc_",gsub('/',':',insignat),".nex",sep='')
		cmd		<- paste(cmd," ",tmp,"\n",sep='')
	}
	else if(hpcsys=="cx1.hpc.ic.ac.uk")		#imperial - use scratch directory
	{
		tmp		<- paste(tmpdir,'/',infile,'_',gsub('/',':',insignat),".timetrees",sep='')
		cmd		<- paste(cmd,prog.beastmcc[hpcsys]," -burnin ",beastmcc.burnin," -heights ",beastmcc.heights," ",tmp,sep='')
		tmp		<- paste(tmpdir,'/',infile,"_mcc_",gsub('/',':',insignat),".nex",sep='')
		cmd		<- paste(cmd," ",tmp,"\n",sep='')
		cmd		<- paste(cmd,"cp -f ",tmp," ", indir,'\n',sep='')
	}
	cmd		<- paste(cmd,"echo \'end ",prog.beastmcc[hpcsys],"\'\n",sep='')
	cmd<- paste(cmd,"#######################################################
# end: run TREEANNOTATOR
#######################################################\n",sep='')		
	cmd
}

#' @export
hivc.cmd.beast2.runxml<- function(indir, infile, insignat, prog.beast=PR.BEAST2, prog.opt.Xms="64m", prog.opt.Xmx="400m", hpc.tmpdir.prefix="beast2", hpc.ncpu=1)
{
	cmd		<- "#######################################################
# start: run BEAST2
#######################################################"	
	hpcsys	<- hivc.get.hpcsys()
	#hpcsys<- "cx1.hpc.ic.ac.uk"
	cmd		<- paste(cmd,paste("\necho \'run ",prog.beast,"\'\n",sep=''))
	if(hpcsys=="debug")						#my MAC - don t use scratch
	{		
		tmp		<- paste(indir,'/',infile,'_',gsub('/',':',insignat),".xml",sep='')
		cmd		<- paste(cmd,"java -Xms",prog.opt.Xms," -Xmx",prog.opt.Xmx," -jar ",prog.beast," -overwrite -working -threads ",hpc.ncpu," ",tmp,'\n',sep='')
	}
	else if(hpcsys=="cx1.hpc.ic.ac.uk")		#imperial - use scratch directory
	{
		cmd		<- paste(cmd,"CWD=$(pwd)\n",sep='')
		cmd		<- paste(cmd,"echo $CWD\n",sep='')		
		tmpdir	<- paste("$CWD/",hpc.tmpdir.prefix,'_',format(Sys.time(),"%y-%m-%d-%H-%M-%S"),sep='')
		cmd		<- paste(cmd,"mkdir -p ",tmpdir,'\n',sep='')
		tmp		<- paste(indir,'/',infile,'_',gsub('/',':',insignat),".xml",sep='')
		cmd		<- paste(cmd,"cp ",tmp," ",tmpdir,'\n',sep='')
		tmp		<- paste(tmpdir,'/',infile,'_',gsub('/',':',insignat),".xml",sep='')
		cmd		<- paste(cmd,"java -Xms",prog.opt.Xms," -Xmx",prog.opt.Xmx," -jar ",prog.beast," -working -threads ",hpc.ncpu," ",tmp,'\n',sep='')	
		cmd		<- paste(cmd,"cp -f ",tmpdir,"/* ", indir,'\n',sep='')		
	}
	cmd		<- paste(cmd,"echo \'end ",prog.beast,"\'\n",sep='')
	cmd		<- paste(cmd,"#######################################################
# end: run BEAST2
#######################################################\n",sep='')
	cmd
}

#add additional high performance computing information 
#' @export
hivc.cmd.hpcwrapper<- function(cmd, hpcsys= hivc.get.hpcsys(), hpc.walltime=24, hpc.mem=HPC.MEM, hpc.nproc=HPC.NPROC[hpcsys], hpc.q=NA)
{
	wrap<- "#!/bin/sh"
	#hpcsys<- HPC.CX1.IMPERIAL
	if(hpcsys==HPC.CX1.IMPERIAL)
	{				
		tmp	<- paste("#PBS -l walltime=",hpc.walltime,":59:59,pcput=",hpc.walltime,":45:00",sep='')
		wrap<- paste(wrap, tmp, sep='\n')		
		tmp	<- paste("#PBS -l select=1:ncpus=",hpc.nproc,":mem=",hpc.mem,sep='')
		wrap<- paste(wrap, tmp, sep='\n')
		wrap<- paste(wrap, "#PBS -j oe", sep='\n')
		if(!is.na(hpc.q))
			wrap<- paste(wrap, paste("#PBS -q",hpc.q), sep='\n\n')
		wrap<- paste(wrap, HPC.LOAD, sep='\n')
	}
	else if(hpcsys=='debug')
		cat(paste("\ndetected no HPC system and no hpcwrapper generated, domain name is",hpcsys))
	else
		stop(paste("unknown hpc system with domain name",hpcsys))
	
	cmd<- lapply(seq_along(cmd),function(i){	paste(wrap,cmd[[i]],sep='\n')	})
	if(length(cmd)==1)
		cmd<- unlist(cmd)
	cmd	
}

#create high performance computing qsub file and submit
#' @export
hivc.cmd.hpccaller<- function(outdir, outfile, cmd)
{
	if( nchar( Sys.which("qsub") ) )
	{
		file	<- paste(outdir,'/',outfile,'.qsub',sep='')
		cat(paste("\nwrite HPC script to",file,"\n"))
		cat(cmd,file=file)
		cmd		<- paste("qsub",file)
		cat( cmd )
		cat( system(cmd, intern=TRUE) )
		Sys.sleep(1)
	}
	else
	{
		file	<- paste(outdir,'/',outfile,'.sh',sep='')
		cat(paste("\nwrite Shell script to\n",file,"\nStart this shell file manually\n"))
		cat(cmd,file=file)
		Sys.chmod(file, mode = "777")		
	}
	
}