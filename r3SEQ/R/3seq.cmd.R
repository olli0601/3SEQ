#' This file contains R functions that provide an interface to HPC scripting, and the 3SEQ and ExaML programs
#' Shell scripts are generated that can be either run directly, or submitted to an HPC system.
#' At present, cmd.hpccaller provides an interface to the CX1B high performance system that runs PBS (portable batch system)

PR.PACKAGE					<- "r3seq"
PR.STARTME					<- system.file(package=PR.PACKAGE, "misc", "3seq.startme.R") 
PR.RECOMB.3SEQ				<- system.file(package=PR.PACKAGE, "ext", "3seq") 
PR.RECOMB.PROCESS3SEQOUTPUT	<- paste(PR.STARTME,"-exe=RECOMB.PROCESS3SEQOUT",sep=' ')
PR.RECOMB.CHECKCANDIDATES	<- paste(PR.STARTME,"-exe=RECOMB.CHECKCANDIDATES",sep=' ')
PR.RECOMB.PLOTINCONGRUENCE	<- paste(PR.STARTME,"-exe=RECOMB.PLOTINCONGRUENCE",sep=' ')
HPC.NPROC					<- {tmp<- c(1,4); names(tmp)<- c("debug","cx1.hpc.ic.ac.uk"); tmp}
HPC.MPIRUN					<- {tmp<- c("mpirun","mpiexec"); names(tmp)<- c("debug","cx1.hpc.ic.ac.uk"); tmp}
HPC.CX1.IMPERIAL			<- "cx1.hpc.ic.ac.uk"		#this is set to system('domainname',intern=T) for the hpc cluster of choice
HPC.MEM						<- "1750mb"
HPC.CX1.IMPERIAL.LOAD		<- "module load intel-suite mpi R/2.15"

#generate 3seq command
#' @export
r3seq.cmd.run.3seq<- function(infile, outfile=paste(infile,".3s.rec",sep=''), recomb.3seq.siglevel=0.1, recomb.3seq.testvsall.beginatseq=NA, recomb.3seq.testvsall.endatseq=NA, prog= PR.RECOMB.3SEQ, nproc=1, verbose=1)
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
r3seq.cmd.process.3SEQ.output<- function(indir, infile, insignat, prog= PR.RECOMB.PROCESS3SEQOUTPUT, resume=1, verbose=1)
{
	cmd<- "#######################################################
# start: run r3seq.prog.process.3SEQ.output
#######################################################"
	cmd<- paste(cmd,paste("\necho \'run ",prog,"\'\n",sep=''))
	#default commands
	cmd<- paste(cmd,prog," -v=",verbose," -resume=",resume,sep='')
	cmd<- paste(cmd," -indir=",indir," -infile=",infile," -insignat=",insignat,sep='')
	#verbose stuff
	cmd<- paste(cmd,paste("\necho \'end ",prog,"\'\n",sep=''))
	cmd<- paste(cmd,"#######################################################
# end: run r3seq.prog.process.3SEQ.output
#######################################################\n",sep='')
	cmd	
}

#check 3seq candidate recombinants
#' @export
r3seq.cmd.check.candidates<- function(indir, infile, insignat, triplet.id, prog= PR.RECOMB.CHECKCANDIDATES, resume=1, verbose=1,hpc.walltime=NA, hpc.q=NA, hpc.mem=NA, hpc.nproc=NA)
{
	cmd<- "#######################################################
# start: r3seq.prog.get.incongruence
#######################################################"
	cmd<- paste(cmd,paste("\necho \'run ",prog,"\'\n",sep=''))
	#default commands
	cmd<- paste(cmd,prog," -v=",verbose," -resume=",resume,sep='')
	cmd<- paste(cmd," -indir=",indir," -infile=",infile," -insignat=",insignat," -tripletid=",triplet.id,sep='')	
	if(!is.na(hpc.walltime))	cmd<- paste(cmd," -hpc.walltime=",hpc.walltime,sep='')
	if(!is.na(hpc.q))			cmd<- paste(cmd," -hpc.q=",hpc.q,sep='')
	if(!is.na(hpc.mem))			cmd<- paste(cmd," -hpc.mem=",hpc.mem,sep='') 
	if(!is.na(hpc.nproc))		cmd<- paste(cmd," -hpc.nproc=",hpc.nproc,sep='')	
	#verbose stuff
	cmd<- paste(cmd,paste("\necho \'end ",prog,"\'\n",sep=''))
	cmd<- paste(cmd,"#######################################################
# end: r3seq.prog.get.incongruence
#######################################################\n",sep='')
	cmd	
}

r3seq.cmd.plot.incongruence<- function(indir, infile, insignat, triplet.id=NA, prog= PR.RECOMB.PLOTINCONGRUENCE, opt.select=NA,verbose=1)
{
	cmd<- "#######################################################
# start: r3seq.prog.plot.incongruence
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
# end: r3seq.prog.plot.incongruence
#######################################################\n",sep='')
	cmd	
}

######################################################################################
cmd.hpcsys<- function()
{
	tmp<- system('domainname',intern=T)
	if(!nchar(tmp))	tmp<- "debug"
	tmp
}

#' @export
cmd.hpcwrapper.cx1.ic.ac.uk<- function(hpc.walltime=24, hpc.mem=HPC.MEM, hpc.nproc=1, hpc.q=NA)
{
	wrap<- "#!/bin/sh"
	tmp	<- paste("#PBS -l walltime=",hpc.walltime,":59:59,pcput=",hpc.walltime,":45:00",sep='')
	wrap<- paste(wrap, tmp, sep='\n')		
	tmp	<- paste("#PBS -l select=1:ncpus=",hpc.nproc,":mem=",hpc.mem,sep='')
	wrap<- paste(wrap, tmp, sep='\n')
	wrap<- paste(wrap, "#PBS -j oe", sep='\n')
	if(!is.na(hpc.q))
		wrap<- paste(wrap, paste("#PBS -q",hpc.q), sep='\n\n')
	wrap<- paste(wrap, HPC.CX1.IMPERIAL.LOAD, sep='\n')
	wrap
}

#add additional high performance computing information 
#' @export
cmd.hpcwrapper<- function(cmd, hpc.sys= cmd.hpcsys(), hpc.walltime=24, hpc.mem=HPC.MEM, hpc.nproc=1, hpc.q=NA)
{	
	#hpc.sys<- HPC.CX1.IMPERIAL
	if(hpc.sys==HPC.CX1.IMPERIAL)
		wrap<- cmd.hpcwrapper.cx1.ic.ac.uk(hpc.walltime=hpc.walltime, hpc.mem=hpc.mem, hpc.nproc=hpc.nproc, hpc.q=hpc.q)
	else
	{
		wrap<- "#!/bin/sh"
		cat(paste("\ndetected no HPC system and no hpcwrapper generated, domain name is",hpc.sys))
	}
	cmd<- lapply(seq_along(cmd),function(i){	paste(wrap,cmd[[i]],sep='\n')	})
	if(length(cmd)==1)
		cmd<- unlist(cmd)
	cmd	
}

#create high performance computing qsub file and submit
#' @export
cmd.hpccaller<- function(outdir, outfile, cmd)
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
		cat(paste("\nwrite Shell script to\n",file,"\nNo 'qsub' function detected to submit the shell script automatically.\nStart this shell file manually\n"))
		cat(cmd,file=file)
		Sys.chmod(file, mode = "777")		
	}
	
}