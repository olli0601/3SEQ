#' this file contains all R functions of the R package
#' @import ape
#' @import data.table



#' @export
hivc.seq.create.referencepairs<- function(dir.name= DATA)
{
	if(0)	#generate ATHENA_2013_hptn052.rda
	{
		gban				<- c( paste("JN",seq.int(247047,247075),sep=''),paste("JN",seq.int(634296,634492),sep='') )	
		tmp					<- hivc.read.GenBank(gban, as.character=0, attributes= c("isolate","country","collection_date"))		
		file				<- "ATHENA_2013_hptn052.fa"
		write.dna(tmp, paste(DATA,file,sep='/'), format = "fasta" )
		cmd					<- hiv.cmd.clustalo(paste(dir.name,"tmp",sep='/'), file, signat='', outdir=paste(dir.name,"tmp",sep='/'))
		cat(cmd)
		if(0) system(cmd)		
		file				<- paste(DATA,"tmp/ATHENA_2013_hptn052.fa.clustalo",sep='/')
		hptn052				<- read.dna(file, format = "fasta" )		
		rownames(hptn052)	<- attr(tmp,"isolate")		
		file				<- paste(DATA,"tmp/ATHENA_2013_hptn052.rda",sep='/')
		save(hptn052, file= file)
	}
	file<- paste(DATA,"tmp/ATHENA_2013_hptn052.rda",sep='/')		
	load(file)
	#select control sequences and compute pairwise distances between them
	hptn052.control		<- grep("control", rownames(hptn052))	
	hptn052.control.d	<- hivc.pwdist( hptn052[hptn052.control,] )
	hptn052.control.d	<- hptn052.control.d[ upper.tri(hptn052.control.d) ]
	#select positive sequences and compute pairwise distance between them	
	hptn052.sdc			<- c(grep("A", rownames(hptn052)), grep("B", rownames(hptn052)))
	hptn052.sdc			<- hptn052[hptn052.sdc,]
	hptn052.sdc.ind		<- sapply(strsplit(rownames(hptn052.sdc),'-'), function(x)
			{ 
				as.numeric(x[2]) + ifelse(x[3]=='I',0,0.5)
			})
	hptn052.sdc.ind		<- sapply( unique( hptn052.sdc.ind ), function(x){		hptn052.sdc[hptn052.sdc.ind==x,]	})
	hptn052.sdc.ind		<- lapply( which( sapply(hptn052.sdc.ind,nrow)==2 ), function(i)  hptn052.sdc.ind[[i]] )
	hptn052.sdc.ind.d	<- sapply( hptn052.sdc.ind, function(x) hivc.pwdist(x)[1,2] )
	
	list( gd.islink= hptn052.sdc.ind.d, gd.unlinked= hptn052.control.d )
}

#' @export
hivc.seq.write.dna.nexus<- function(seq.DNAbin.mat, ph=NULL, file=NULL, nexus.format="DNA",nexus.gap='-', nexus.missing='?', nexus.interleave="NO")
{		
	tmp		<- cbind( rownames(seq.DNAbin.mat), apply( as.character( seq.DNAbin.mat ), 1, function(x) paste(x,sep='',collapse='')  ) )
	tmp		<- apply(tmp, 1, function(x) paste(x, collapse='\t', sep=''))
	tmp		<- paste(tmp, collapse='\n',sep='')
	header	<- paste( "#NEXUS\nBEGIN DATA;\nDIMENSIONS NTAX=",nrow(seq.DNAbin.mat)," NCHAR=",ncol(seq.DNAbin.mat),";\nFORMAT DATATYPE=",nexus.format," MISSING=",nexus.missing," GAP=",nexus.gap," INTERLEAVE=",nexus.interleave,";\nMATRIX\n", collapse='',sep='')
	tmp		<- paste(header, tmp, "\n;\nEND;\n", sep='')
	if(!is.null(ph))
	{
		tmp	<- paste(tmp,"#BEGIN TAXA;\nTAXLABELS ", paste(ph$tip.label, sep='',collapse=' '), ';\nEND;\n\n',sep='')
		tmp	<- paste(tmp, "BEGIN TREES;\nTREE tree1 = ", write.tree(ph), "\nEND;\n", sep='')
	}
	if(!is.null(file))
		cat(tmp, file=file)
	tmp
}		

#' @export
hivc.seq.write.dna.phylip<- function(seq.DNAbin.mat, file)
{		
	tmp<- cbind( rownames(seq.DNAbin.mat), apply( as.character( seq.DNAbin.mat ), 1, function(x) paste(x,sep='',collapse='')  ) )
	tmp<- paste(t(tmp),collapse='\n',sep='')	
	tmp<- paste( paste(c(nrow(seq.DNAbin.mat),ncol(seq.DNAbin.mat)),sep='',collapse=' '),'\n',tmp,'\n',collapse='',sep='' )
	cat(tmp, file=file)
}

hivc.seq.find<- function(char.matrix, pos0= NA, from= c(), verbose=1)
{
	if(is.na(pos0)) 	stop("start position of token to be replaced is missing")
	if(!length(from))	stop("token to be replaced is missing")
	query.colidx	<- seq.int(pos0,pos0+length(from)-1)
	query.yes		<- which( apply(char.matrix, 1, function(x)	all(x[query.colidx]==from) ) )
	query.yes	
}

#' @export
hivc.seq.unique<- function(seq.DNAbin.matrix)
{
	x<- as.character(seq.DNAbin.matrix)
	x<- apply(x, 1, function(z) paste(z,collapse=''))
	seq.DNAbin.matrix[!duplicated(x),]			
}

#' @export
hivc.seq.dist<- function(seq.DNAbin.matrix, verbose=1)
{
	if(0)
	{
		require(ape)
		ans<- dist.dna(seq.DNAbin.matrix, model="raw", as.matrix=1)
	}
	if(1)
	{
		library(bigmemory)
		options(bigmemory.typecast.warning=FALSE)
		big.matrix.charmax<- 127
		dummy	<- 0
		ans<- big.matrix(nrow(seq.DNAbin.matrix),nrow(seq.DNAbin.matrix), dimnames=list(rownames(seq.DNAbin.matrix),c()), type="char", init=NA)
		if(nrow(seq.DNAbin.matrix)>5)
		{
			ans[1,1:6]<- 2^seq.int(6,11)-1
			if(is.na(ans[1,2]))		stop("unexpected behaviour of bigmemory")
			ans[1,1:6]<- NA
		}						
		for(i1 in seq.int(1,nrow(seq.DNAbin.matrix)-1))
		{			
			seq1<- seq.DNAbin.matrix[i1,]
			time<- system.time	(
							tmp	<- 1 - sapply(seq.int(i1+1,nrow(seq.DNAbin.matrix)),function(i2){		.C("hivc_dist_ambiguous_dna", seq1, seq.DNAbin.matrix[i2,], ncol(seq1), dummy )[[4]]			})
						)[3]
			if(verbose)	cat(paste("\ncompute distance of row",i1,"entries",nrow(seq.DNAbin.matrix)-i1,"took",time))			
			tmp												<- round(tmp*1e3,d=0)			
			tmp[tmp>big.matrix.charmax]						<- big.matrix.charmax
			ans[i1, seq.int(i1+1,nrow(seq.DNAbin.matrix))]	<- tmp
		}		
	}
	ans
}

hivc.seq.replace<- function(seq.DNAbin.matrix, code.from='?', code.to='n', verbose=0)
{
	seq.DNAbin.matrix	<- as.character(seq.DNAbin.matrix)		
	seq.DNAbin.matrix	<- apply(seq.DNAbin.matrix, 2, function(col) 		gsub(code.from,code.to,col,fixed=1)			)	
	as.DNAbin( seq.DNAbin.matrix )
}

#' @export
hivc.clu.geneticdist.cutoff<- function(dir.name= DATA, plot=1, verbose=1, level.retain.unlinked=0.05)
{	
	refpairs			<- hivc.seq.create.referencepairs(dir.name)
	refpairs			<- lapply(refpairs,function(x) x*100)
	
	xlim				<- range(c(refpairs[[1]], refpairs[[2]]))
	width				<- c(1.5,0.75)
	breaks				<- seq(xlim[1],xlim[2],len=40)	
	refpairs.h			<- lapply(seq_along(refpairs),function(i)
			{ 				
				tmp<- hist(refpairs[[i]],breaks=breaks,plot=0)
				tmp$kde<- density(refpairs[[i]], kernel="biweight",from=breaks[1],to=breaks[length(breaks)],width = max(EPS,width[i]*diff(summary(refpairs[[i]])[c(2,5)])))
				tmp
			})
	names(refpairs.h)	<- names(refpairs)
	ylim				<- range(c(refpairs.h[[1]]$intensities,refpairs.h[[2]]$intensities))
	cols				<- c(my.fade.col("black",0.2),my.fade.col("black",0.6))
	plot(1,1,type='n',xlim=xlim,ylim=ylim,xlab="% genetic distance",ylab="density")		
	lapply(seq_along(refpairs.h),function(i)
			{
				plot(refpairs.h[[i]],add=1,freq=0,col=cols[i],border=NA)				
			})
	lapply(seq_along(refpairs.h),function(i)
			{
				lines(refpairs.h[[i]][["kde"]]$x,refpairs.h[[i]][["kde"]]$y,col=cols[i])								
			})
	
	ref.gd.unlinked.cdf	<- cumsum( refpairs.h[["gd.unlinked"]][["kde"]]$y / sum(refpairs.h[["gd.unlinked"]][["kde"]]$y) )		
	gd.cutoff			<- refpairs.h[["gd.unlinked"]][["kde"]]$x[ which.max( ref.gd.unlinked.cdf[ref.gd.unlinked.cdf<=level.retain.unlinked] ) ]
	abline(v=gd.cutoff,lty=2)
	gd.specificity<- refpairs.h[["gd.islink"]][["kde"]]$y / sum(refpairs.h[["gd.islink"]][["kde"]]$y)
	gd.specificity<- sum( gd.specificity[ refpairs.h[["gd.islink"]][["kde"]]$x<=gd.cutoff ] )
	
	if(verbose)	cat(paste("\ncutoff for high sensitivity of",1-level.retain.unlinked,"of detecting unlinked is",gd.cutoff))
	if(verbose)	cat(paste("\nspecificity for detecting linked at cutoff",gd.cutoff," is",gd.specificity))
	
	list(gd.cutoff= gd.cutoff, gd.specificity=gd.specificity)	
}
######################################################################################
hivc.clu.collapse.monophyletic.withinpatientseq<- function(cluphy.subtrees, df.cluinfo, verbose=1)
{
	cluphy.subtrees.names	<- names(cluphy.subtrees)		
	cluphy.subtrees			<- lapply(seq_along(cluphy.subtrees), function(i)
								{
									#i			<- 409
									x			<- cluphy.subtrees[[i]]
									setkey(df.cluinfo, "FASTASampleCode")		
									tmp		<- subset(df.cluinfo[x$tip.label, ][, list(nseq=length(FASTASampleCode), FASTASampleCode=FASTASampleCode, PosSeqT=PosSeqT), by="Patient"], nseq>1, c(Patient, FASTASampleCode, PosSeqT))	#determine Patients with >1 seq in cluster
									if(nrow(tmp))	#check if any within patient seqs in this cluster are monophyletic
									{
										tmp			<- tmp[,	{
																	z	<- extract.clade(x, hivc.clu.mrca(x, FASTASampleCode), root.edge=1)
																	list( FASTASampleCode=FASTASampleCode, PosSeqT=PosSeqT, isMonophyletic= Ntip(z)==length(FASTASampleCode), root.edge= z$root.edge)
																},by="Patient"]
										tmp			<- subset(tmp, isMonophyletic)
									}
									if(nrow(tmp))	#drop all within patient seqs from cluster except the one closest to root
									{			
										x.d2r		<- node.depth.edgelength(x)
										x.rootedge	<- x$root.edge
										names(x.d2r)<- x$tip.label			 
										tmp			<- tmp[,  list(FASTASampleCode=FASTASampleCode, root.edge=root.edge, order=order(x.d2r[FASTASampleCode], PosSeqT))	,by="Patient"]
										x			<- drop.tip(x,subset(tmp,order>1)[,FASTASampleCode])						#drop within patient seqs furthest away from root
										x$root.edge	<- x.rootedge
										tmp			<- subset(tmp,order==1)
										x.tip		<- match( tmp[,FASTASampleCode], x$tip.label)			
										x$edge.length[ sapply( x.tip, function(z) which(x$edge[,2]==z) ) ]	<- tmp[,root.edge]	#set edge.length of new tip to root length of within patient clade			
									}
									x								
								})		
	names(cluphy.subtrees)	<- cluphy.subtrees.names
	if(verbose) cat(paste("\nnumber of clusters after dropping monophyletic within patient seqs (Should not change), n=",length(cluphy.subtrees)))
	#collect tip labels and update df.cluinfo
	cluphy.df				<- merge(data.table( FASTASampleCode= unlist( lapply(cluphy.subtrees, function(x)	x$tip.label ) ) ), df.cluinfo, by="FASTASampleCode")
	if(verbose) cat(paste("\nnumber of seq in clusters after dropping monophyletic within patient seqs, n=",nrow(cluphy.df)))
	cluphy					<- hivc.clu.polyphyletic.clusters(cluphy.df, cluphy.subtrees)$cluphy
	cluphy.clustering		<- hivc.clu.clusterbycluphy(cluphy.subtrees, cluphy)		
	list(cluphy=cluphy, cluphy.df=cluphy.df, cluphy.clustering=cluphy.clustering, cluphy.subtrees=cluphy.subtrees)
}	
######################################################################################
hivc.clu.clusterbycluphy<- function(cluphy.subtrees, cluphy)	
{
	require(phangorn)
	cluphy.subtrees.names					<- as.numeric(names(cluphy.subtrees))	
	clustering								<- list()	
	clustering[["ntips"]]					<- Ntip(cluphy)
	tmp										<- rbindlist( lapply( seq_along(cluphy.subtrees), function(i)
												{
													x			<- cluphy.subtrees[[i]]
													x.tip		<- sapply(cluphy.subtrees[[i]]$tip.label, function(z) 	match(z,cluphy$tip.label) )
													rev(Ancestors(cluphy, x.tip[1]))[2]
													x.tip.anc	<- unique( unlist( lapply(x.tip, function(z) rev(Ancestors(cluphy, z))[-1] ) ) )
													data.table(Node=c(x.tip,x.tip.anc), MRCA=x.tip.anc[1],cluster=cluphy.subtrees.names[i] )
												}))
	setkey(tmp, Node)	
	clustering[["clu.mem"]]					<- rep(NA, Nnode(cluphy,internal=0))
	clustering[["clu.mem"]][tmp[,Node]]		<- tmp[,cluster]					
	clustering[["clu.idx"]]					<- rep(NA, max(tmp[,cluster],na.rm=1))
	clustering[["clu.idx"]][tmp[,cluster]]	<- tmp[,MRCA]
	clustering[["size.tips"]]				<- table(table(subset(tmp, Node<=Ntip(cluphy))[,cluster]))
	clustering
}
######################################################################################
#' Compute a set of phylogenetic clusters according to several criteria; taken from http://dx.doi.org/10.6084/m9.figshare.97225
#' Perform DFS. Label a subtree a cluster if both
#'   -its median pairwise patristic ditance (MPPD) is below a percentile of the whole-tree p-distance and
#'   -it is not a leaf. 
#' @export
#' @param ph					phylogenetic tree with branch lengths
#' @param thresh.brl 			MPPD threshold
#' @param dist.brl 				precomputed MPPD values for each node
#' @param thresh.nodesupport	boostrap reliability threshold
#' @param nodesupport			precomputed bootstrap values
#' @param retval				evaluate MPPD by tip nodes only or for all nodes
#' @return vector indicating, for each internal node, which cluster it belongs to
hivc.clu.clusterbythresh<- function(ph,thresh.brl=NULL,dist.brl=NULL,thresh.nodesupport=NULL,nodesupport=NULL,retval="all")
{
	require(ape)
	require(igraph)
	require(geiger)	
	if(all( is.null(c(thresh.brl, thresh.nodesupport))) ) 		stop("all threshold criteria NULL - provide at least one")
	if(!is.null(thresh.nodesupport) && is.null(nodesupport))	stop("node support threshold set but no node support values provided")
	if(!is.null(thresh.brl) && is.null(dist.brl))				stop("branch length threshold set but no branch length distances provided")
	if(!is.null(nodesupport) && any(nodesupport>1+EPS))			warning("Found nodesupport values above 1")
	if(is.null(thresh.brl))
	{		
		dist.brl			<- rep(1,Nnode(ph))
		thresh.brl			<- 1
	}	
	if(is.null(thresh.nodesupport))
	{
		nodesupport			<- rep(1,Nnode(ph))
		thresh.nodesupport	<- 1		
	}
#print(nodesupport[1:10]); print(thresh.nodesupport)
	## set up clustering
	ntips		<- Ntip(ph)	
	clu.i 		<- 0 ## cluster number
	clu.mem 	<- rep(NA,ntips+ph$Nnode) ## cluster member assignment
	clu.idx		<- rep(NA,ntips+ph$Nnode) ## cluster index assignment
	igraph.ph	<- graph.edgelist(ph$edge) ## ph in igraph form
	dfs 		<- graph.dfs(igraph.ph,root=ntips+1,neimode='out',order=TRUE,dist=TRUE)
#print(c(length(dfs$order),Ntip(ph),Nnode(ph)))
	## travese the ph in depth first order
	for(i in 1:length(dfs$order))
	{
		node <- dfs$order[i]
#print( c(node,node-ntips, is.na(clu.mem[node]), thresh.brl, dist.brl[node-ntips], dist.brl[node-ntips]<=thresh.brl, thresh.nodesupport, nodesupport[node-ntips], nodesupport[node-ntips]<thresh.nodesupport) )
		if(	node > ntips	&&											## skip leaves
			is.na(clu.mem[node]) &&										## only consider unassigned nodes
			nodesupport[node-ntips]>=thresh.nodesupport-EPS &&
			dist.brl[node-ntips]<=thresh.brl+EPS
			)	 
		{				
#print(nodesupport[node-ntips]);print(c(node,i))			
			clu.i 			<- clu.i+1
			subtree 		<- graph.dfs(igraph.ph,node,neimode='out',unreachable=FALSE)$order
			subtree 		<- subtree[! is.nan(subtree)]
			clu.mem[subtree]<- clu.i
			clu.idx[node]	<- clu.i
		}
	}
	ans <- list( clu.mem=clu.mem, clu.idx=which(!is.na(clu.idx)), size.all=table(clu.mem), size.tips=table(clu.mem[1:ntips]), ntips=ntips, thresh.brl=thresh.brl, thresh.nodesupport=thresh.nodesupport)
	if(retval=="tips")
		ans$clu.mem 		<- ans$clu.mem[1:ntips]
	return(ans)
}
######################################################################################
hivc.clu.exp.typeIerror.randomclu<- function(ph, dist.brl, nodesupport, ph.unlinked.seroneg, ph.unlinked.dead, thresh.brl, thresh.bs)
{
	ph.tips.n				<- Ntip(ph)
	#seroneg.n				<- nrow(ph.unlinked.seroneg)
	#ph.unlinked.seroneg.v	<- ph.unlinked.seroneg[,PhNode]
	clustering				<- hivc.clu.clusterbythresh(ph, thresh.nodesupport=thresh.bs, thresh.brl=thresh.brl, dist.brl=dist.brl, nodesupport=nodesupport, retval="all")
	clusters				<- unique(clustering[["clu.mem"]][!is.na(clustering[["clu.mem"]])])
	clu.fp.exp	<- sapply(clusters,function(clu)		
			{		
				tmp							<- which(clustering[["clu.mem"]]==clu)
				clu.leaves					<- tmp[ tmp<=ph.tips.n ]
				clu.leaves.seroneg			<- sapply(clu.leaves, function(x)
						{
							as.numeric( subset(ph.unlinked.seroneg, x==PhNode, c(PhNodeUnlinked,NegT) ) )	 	#O(logn) search; returns NA if subset dt is empty																 
						})				
#print(clu); print(tmp); print(clu.leaves); print(clu.leaves.seroneg)
				#if all leaves not seroneg, return 0								
				if(all(is.na(clu.leaves.seroneg[1,])))	return(0)
				#if at least one seroneg, compute prob that 1-P(X=0) where X is number of seroconverters dead at SC of latest seroconverter				
				clu.leaves.seroneg			<- rbind(clu.leaves,clu.leaves.seroneg)[,!is.na(clu.leaves.seroneg[2,]),drop=0]								
				clu.leaves.latestseroneg.ix	<- which.max( clu.leaves.seroneg[3,] )
				tmp							<- sapply(clu.leaves.seroneg[2,],function(x) length(ph.unlinked.dead[[x]]) )
				latest.seroneg.unlinked		<- ph.unlinked.dead[[ clu.leaves.seroneg[2,clu.leaves.latestseroneg.ix] ]]				
				if( max(tmp)!=length(latest.seroneg.unlinked) )	
					stop("expect latest seroneg to have most unlinked HIV+ dead associated")	#which.max does not work because .ix is not necessarily the first																
				#m is "number false pos" is length(latest.seroneg.unlinked)
				#n is "number false pos not known" is ph.tips.n-m
				#k is "number of draws" is length(clu.leaves)-1 
				ans<- 1-dhyper(0, length(latest.seroneg.unlinked), ph.tips.n-length(latest.seroneg.unlinked), length(clu.leaves)-1)
				#c(ans, length(clu.leaves)-1, length(latest.seroneg.unlinked), ph.tips.n)
				ans
			})
	list(fp.data=clu.fp.exp, fp.n=sum(clu.fp.exp), fp.rate= sum(clu.fp.exp)/length(clusters), clustering=clustering)
}
######################################################################################
hivc.clu.truepos<- function(clustering, ph.linked, ph.tips.n, verbose=0)
{
	clusters	<- unique(clustering[["clu.mem"]][!is.na(clustering[["clu.mem"]])])
	tmp			<- which(!is.na(clustering[["clu.mem"]][seq_len(ph.tips.n)]))
	clu.linked	<- merge( data.table(Node=tmp), ph.linked, by="Node" )	
	clu.linked	<- clu.linked[,list(nTP.clu= length(FASTASampleCode)), by="Patient"]
	clu.linked	<- merge(ph.linked, clu.linked, all.x=1, by="Patient")
	set(clu.linked, which(is.na(clu.linked[,nTP.clu])) ,"nTP.clu", 1L )
	setkey(clu.linked, Node)
	
	if(!length(clusters))
	{
		clu.tp		<- cbind( subset(clu.linked,select=c(Patient,nTP,nTP.clu)), TP.clu=0)
		setnames(clu.tp, c("nTP","nTP.clu"),c("TP.n","TP.n.clu"))				
	}
	else
	{				
		clu.tp	<- lapply(c(-1,clusters), function(clu)		
			{		
				#clu<- 48
				if(clu==-1)					#bug in data.table; if first data.table empty, rbind fails		http://lists.r-forge.r-project.org/pipermail/datatable-help/2012-November/001372.html
					return( as.data.table(data.frame(Patient="XXX",TP.clu=0,TP.n=0,TP.n.clu=0,clu=0,clu.tips=0,stringsAsFactors=F)) )
				tmp							<- which(clustering[["clu.mem"]]==clu)
				clu.tips					<- tmp[ tmp<=ph.tips.n ]				
				clu.tips.linked				<- na.omit(merge( data.table(Node=clu.tips, key="Node"), clu.linked, all.x=1))	#extract clu.tips that should be linked
				if(!nrow(clu.tips.linked))	
					ans						<- as.data.table(na.omit(data.frame(Patient=NA,TP.clu=NA,TP.n=NA, TP.n.clu=NA, clu=NA, clu.tips=NA)))
				else	
				{
					ans						<- clu.tips.linked[, list(TP.clu=length(FASTASampleCode), TP.n=nTP[1], TP.n.clu=nTP.clu[1], clu=clu, clu.tips=length(clu.tips) ), by=Patient]
					set(ans, NULL, "Patient", as.character(ans[,Patient]))
				}
				ans
			})
		#rbind and remove bugfix first row
		clu.tp		<- rbindlist(clu.tp)[-1,]	
		
	}	
	# account for different seq numbers per patient with a combinatorial argument 
	# -- only important when the number of correctly clustered within-patient sequences is evaluated
	clu.tp[, pairs.TP.n:= choose(TP.n,2)]
	clu.tp[, pairs.TP.n.clu:= choose(TP.n.clu,2)]
	clu.tp[, pairs.TP.clu:= choose(TP.clu,2)]
	# patients in wrong cluster
	clu.diffclu								<- subset( clu.tp[,list(nclu=length(na.omit(unique(clu)))),by="Patient"], nclu>1, Patient )
	clu.diffclu								<- merge( clu.diffclu, clu.tp, all.x=1, by="Patient" )

	# evaluate among all sequences in tree
	patients.allseqinsameclu				<- nrow( subset(clu.tp, TP.clu==TP.n) )							#if TP.clu==TP.n then Patient can occur only once in subset
	patients.total							<- length(unique(clu.linked[,Patient]))
	tmp										<- clu.tp[, list(TP.clu= max(TP.clu), pairs.TP.clu= max(pairs.TP.clu), TP.n=TP.n[1], pairs.TP.n=pairs.TP.n[1]),by=Patient]	#if TP.clu!=TP.n then Patient occurs multiple times
	patients.freqinsameclu					<- sum(tmp[, TP.clu]) / sum(tmp[, TP.n])
	# preferred way to quantify  the number of correctly clustered within-patient sequences
	pairs.freqinsameclu						<- sum(tmp[, pairs.TP.clu]) / sum(tmp[, pairs.TP.n])
	
	# evaluate among all sequences that cluster
	patients.clustering.allseqindiffclu		<- length(unique(clu.diffclu[,Patient]))
	patients.clustering.allseqinsameclu		<- nrow( subset(clu.tp, TP.clu==TP.n.clu & TP.n.clu>1) )
	patients.clustering.total				<- length(unique(subset(clu.linked,nTP.clu>1)[,Patient]))
	tmp										<- subset(clu.tp, TP.n.clu>1)[, list(TP.clu= max(TP.clu), pairs.TP.clu=max(pairs.TP.clu), TP.n.clu=TP.n.clu[1], pairs.TP.n.clu=pairs.TP.n.clu[1]),by=Patient]
	patients.clustering.freqinsameclu		<- sum(tmp[, TP.clu]) / sum(tmp[,TP.n.clu])
	# preferred way to quantify  the number of correctly clustered within-patient sequences	
	pairs.clustering.freqinsameclu			<- sum(tmp[, pairs.TP.clu]) / sum(tmp[,pairs.TP.n.clu])
	
	list(	clu.n				= length(clusters),	
			clu.onlytp			= subset(clu.tp, TP.clu==TP.n),
			clu.missedtp		= subset(clu.tp, TP.clu!=TP.n),
			clu.diffclu 		= clu.diffclu,
			tp.by.all			= patients.allseqinsameclu / patients.total,
			tp.by.sum			= pairs.freqinsameclu,
			tpclu.by.all		= patients.clustering.allseqinsameclu / patients.clustering.total,			
			tpclu.by.diffclu	= (patients.clustering.total - patients.clustering.allseqindiffclu) / patients.clustering.total,
			tpclu.by.sum		= pairs.clustering.freqinsameclu			
			)	
}
######################################################################################
hivc.clu.trueneg<- function(clustering, ph.unlinked.info, ph.unlinked, ph.tips.n, verbose=0)
{
	clusters			<- unique(clustering[["clu.mem"]][!is.na(clustering[["clu.mem"]])])
	clu.fp				<- sapply(clusters, function(clu)		
							{		
								tmp							<- which(clustering[["clu.mem"]]==clu)
								clu.tips					<- tmp[ tmp<=ph.tips.n ]
								clu.tips					<- merge( data.table(Node=clu.tips, key="Node"), subset(ph.unlinked.info,select=c(Node,NodeIdx,NegT)) )						
								if(all( is.na(clu.tips[,NodeIdx]) ))	return(rep(NA,2))
								if(!nrow(clu.tips))	return(rep(NA,2))				
								tmp							<- clu.tips[,NegT]
								tmp							<- ifelse(all(is.na(tmp)),1,which.max(tmp))	
								clu.tips.unlinked			<- ph.unlinked[[ clu.tips[tmp,NodeIdx] ]]			
								#print(clu.tips)
								c(nrow( merge(clu.tips, clu.tips.unlinked, by="Node") ), nrow(clu.tips.unlinked))																					
							})
	colnames(clu.fp)	<- clusters		
	clu.fp				<-	clu.fp[, apply(!is.na(clu.fp),2,all),drop=0]
	list( 	clu.n		= length(clusters), 
			clu.fp		= clu.fp, 
			fpn.by.all	= length(which(clu.fp[1,]>0)), 
			fpn.by.sum	= sum(clu.fp[1,]), 
			fp.by.all	= length(which(clu.fp[1,]>0))/length(clusters)
			)
}
######################################################################################
hivc.clu.clusterbytruepos<- function(ph, dist.brl, nodesupport, ph.linked, thresh.brl=NULL, thresh.bs=NULL, level= 0.95, tol= 0.005, mxit= 10, thresh.bs.lower= min(nodesupport), thresh.brl.upper=max(dist.brl), method.tp="tp.rate.tpclu",verbose=0)
{		
	if(is.null(thresh.brl) && is.null(thresh.bs))	stop("either thresh.bs or thresh.brl required")
	if(!is.null(thresh.brl) && !is.null(thresh.bs))	stop("nothing to do")
	
	error								<- 1
	ph.tips.n							<- Ntip(ph)
	nit									<- 0	
	if(is.null(thresh.bs))
	{
		srch.find.thresh.bs				<- 1
		srch.upper						<- 1
		srch.lower						<- thresh.bs	<- thresh.bs.lower		
	}
	else
	{
		srch.find.thresh.bs				<- 0
		srch.lower						<- thresh.brl 	<- 0
		srch.upper						<- thresh.brl.upper		
	}
	if(verbose)
		cat(paste("\ncalibrate bs?",srch.find.thresh.bs,"\nsrch.lower",srch.lower,"srch.upper",srch.upper))
	#binary search algorithm
	while(	abs(error)>tol && srch.lower<srch.upper && nit<mxit)
	{	
		if(nit)
		{
			if(srch.find.thresh.bs)	
					thresh.bs			<- ( srch.upper + srch.lower ) / 2
			else
					thresh.brl			<- ( srch.upper + srch.lower ) / 2
		}
		clustering					<- hivc.clu.clusterbythresh(ph, thresh.nodesupport=thresh.bs, thresh.brl=thresh.brl, dist.brl=dist.brl, nodesupport=nodesupport, retval="all")
		tmp							<- hivc.clu.truepos(clustering, ph.linked, ph.tips.n)
print(tmp)		
		if(verbose)	
			cat(paste("\nit=",nit,", #clu=",tmp[["clu.n"]],", BS=",thresh.bs,", BRL=",thresh.brl,", avg of %TP in cluster=",round(tmp[[method.tp]],digits=4),sep=''))			
		error						<- tmp[[method.tp]]	- level	
#print(c(srch.upper,srch.lower, thresh.bs, thresh.brl))
		if(srch.find.thresh.bs)	
		{
			if(error<0)
				srch.upper			<- thresh.bs
			else
				srch.lower			<- thresh.bs							
		}
		else
		{
			if(error<0)
				srch.lower			<- thresh.brl				
			else
				srch.upper			<- thresh.brl				
		}		
		nit							<- nit+1
	}
	ans<- list(thresh.bs=thresh.bs, thresh.brl=thresh.brl, clu=clustering, clu.fp.n=error+level, srch.nit= nit, srch.error=error, srch.tol=tol)
	ans		
}
######################################################################################
hivc.clu.clusterbytrueneg<- function(ph, dist.brl, nodesupport, ph.unlinked.info, ph.unlinked, thresh.brl=NULL, thresh.bs=NULL, level= 0.01, tol= 0.005, mxit= 20, thresh.bs.lower= min(nodesupport), thresh.brl.upper=max(dist.brl), verbose=0)
{		
	if(is.null(thresh.brl) && is.null(thresh.bs))	stop("either thresh.bs or thresh.brl required")
	if(!is.null(thresh.brl) && !is.null(thresh.bs))	stop("nothing to do")
	
	error								<- 1
	ph.tips.n							<- Ntip(ph)
	nit									<- 0	
	if(is.null(thresh.bs))
	{
		srch.find.thresh.bs				<- 1
		srch.upper						<- 1
		srch.lower						<- thresh.bs	<- thresh.bs.lower		
	}
	else
	{
		srch.find.thresh.bs				<- 0
		srch.lower						<- 0
		srch.upper						<- thresh.brl	<- thresh.brl.upper		
	}
	if(verbose)
		cat(paste("\ncalibrate bs?",srch.find.thresh.bs,"\nsrch.lower",srch.lower,"srch.upper",srch.upper,"thresh.bs",thresh.bs,"thresh.brl",thresh.brl))
	#binary search algorithm
	while(abs(error)>tol && srch.lower<srch.upper && nit<mxit)
	{	
		if(nit)
		{
			if(srch.find.thresh.bs)	
				thresh.bs			<- ( srch.upper + srch.lower ) / 2
			else
				thresh.brl			<- ( srch.upper + srch.lower ) / 2				
		}
			
		clustering					<- hivc.clu.clusterbythresh(ph, thresh.nodesupport=thresh.bs, thresh.brl=thresh.brl, dist.brl=dist.brl, nodesupport=nodesupport, retval="all")
		tmp							<- hivc.clu.trueneg(clustering, ph.unlinked.info, ph.unlinked, ph.tips.n, verbose=0)
		clustering.fp				<- tmp$clu.fp
		if(verbose)	
			cat(paste("\nit=",nit,", #clu=",tmp[["clu.n"]],", BS=",thresh.bs,", BRL=",thresh.brl,", #unlinked in any cluster=",tmp[["fp.n.sum"]],", #FP=",tmp[["fp.n"]],", %FP=",round(tmp[["fp.rate"]],d=5),sep=''))			
		error						<- tmp[["fp.rate"]]	- level	
#print(c(srch.upper,srch.lower, thresh.bs, thresh.brl))
		if(srch.find.thresh.bs)	
		{
			if(error<0)
				srch.upper			<- thresh.bs
			else
				srch.lower			<- thresh.bs							
		}
		else
		{
			if(error<0)
				srch.lower			<- thresh.brl				
			else
				srch.upper			<- thresh.brl				
		}
		nit							<- nit+1
	}
	ans<- list(thresh.bs=thresh.bs, thresh.brl=thresh.brl, clu=clustering, clu.fp=clustering.fp, srch.nit= nit, srch.error=error, srch.tol=tol)
	ans		
}
######################################################################################
hivc.clu.polyphyletic.clusters<- function(cluphy.df, cluphy.subtrees=NULL, ph=NULL, clustering=NULL, verbose=1, plot.file=NA, pdf.scaley=25, pdf.xlim=NULL, cex.nodelabel=0.2, cex.tiplabel=0.2, adj.tiplabel= c(-0.15,0.5))
{
	if(!is.null(ph) && !is.null(clustering))
	{
		#get node corresponding to index of selected clusters		
		cluphy.cluidx				<- clustering[["clu.idx"]][ unique( cluphy.df[,cluster] ) ]		
		if(any(is.na(cluphy.cluidx)))	stop("unexcpected NA in cluphy.cluidx")
		#get selected clusters into single phylogeny
		cluphy.subtrees				<- lapply(cluphy.cluidx, function(x)		extract.clade(ph, x, root.edge= 1, interactive = FALSE) 		)
		names(cluphy.subtrees)		<- unique( cluphy.df[,cluster] )
	}
	else if(is.null(cluphy.subtrees))
		stop("expect either ('ph' and 'clustering') or 'cluphy.subtrees' as input")
	cluphy							<- eval(parse(text=paste('cluphy.subtrees[[',seq_along(cluphy.subtrees),']]', sep='',collapse='+')))
	print(cluphy)
	plot.coordinates				<- NULL
	#plot selected clusters
	if(!is.na(plot.file))
	{		
		cluphy.tiplabels						<- hivc.clu.get.tiplabels( cluphy, copy(cluphy.df) )
		if(verbose) cat(paste("\nwrite tree to file",plot.file))
		color.except.rootedge					<- rep(1, Nnode(cluphy, internal.only=F))
		color.except.rootedge[Ntip(cluphy)+1]	<- NA
		plot.coordinates						<- hivc.clu.plot(cluphy, color.except.rootedge, file=plot.file, pdf.scaley=pdf.scaley, pdf.off=0, pdf.xlim=pdf.xlim, cex.nodelabel=cex.nodelabel )		
		hivc.clu.plot.tiplabels( seq_len(Ntip(cluphy)), cluphy.tiplabels$text, cluphy.tiplabels$col, cex=cex.tiplabel, adj=adj.tiplabel, add.xinch=0, add.yinch=0 )
		dev.off()
	}	
	list(cluphy=cluphy, cluphy.subtrees=cluphy.subtrees, plot.coordinates=plot.coordinates)
}	
######################################################################################
hivc.clu.brl.bwpat<- function(cluphy.subtrees, cluphy.df, rtn.val.for.na.patient=FALSE)
{
	library(adephylo)
	#get branch length matrices for each subtree 
	cluphy.disttips				<- lapply(cluphy.subtrees, distTips)
	#select branch lengths between seqs of different patients 
	setkey(cluphy.df, FASTASampleCode)
	cluphy.brl.bwpat			<- lapply(cluphy.disttips, function(x)
			{
				x.patient		<- subset( cluphy.df[attr(x, "Labels"),], select=Patient )
				tmp				<- CJ(p1= seq_len(nrow(x.patient)), p2= seq_len(nrow(x.patient)))
				tmp[,lower.tri:=tmp[,p1<p2]]
				x.bwpat.matidx	<- tmp[, lower.tri & x.patient[p1,Patient]!=x.patient[p2,Patient]]		#if Patient ID not know, matidx is NA
				if(rtn.val.for.na.patient)					
					x.bwpat.matidx[is.na(x.bwpat.matidx)]<- TRUE										
				as.matrix(x)[x.bwpat.matidx]															#if Patient ID not know, returns NA
			}) 
	names(cluphy.brl.bwpat)		<- names(cluphy.subtrees)		
	cluphy.brl.bwpat
}
######################################################################################
hivc.clu.getplot.incountry<- function(ph, clustering, df.cluinfo, verbose=1, plot.file= NA, char.frgn='CountryInfection=="FRGN"', char.frgntn='CountryInfection=="FRGNTN"')
{
	clut						<- parse(text=paste("!is.na(CountryInfection) & (",char.frgn," | ",char.frgntn,") ", sep=''))		
	clut						<- table( df.cluinfo[,cluster, eval(clut) ] )
	#
	# select clusters containing only in-country seq
	#	
	clut.nofrgn					<- apply(clut, 2, function(x)	x[1]>0 & x[2]==0		)
	clut.nofrgn					<- as.numeric( colnames(clut[,clut.nofrgn]) )
	if(verbose) cat(paste("\nnumber of clusters with only in-country seq is n=", length(clut.nofrgn)))
	cluphy.nofrgn				<- merge(data.table(cluster=clut.nofrgn), df.cluinfo, by="cluster")
	if(verbose) cat(paste("\nnumber of seq in clusters with only in-country seq is n=", nrow(cluphy.nofrgn)))
	#
	# get subtrees corresponding to clusters containing only in-country seq
	#	
	cluphy.nofrgn.subtrees		<- hivc.clu.polyphyletic.clusters(cluphy.nofrgn , ph=ph, clustering=clustering)$cluphy.subtrees		
	#
	# select clusters containing at least one in-country seq
	#
	clut.mixed					<- apply(clut, 2, function(x)	all(x>0)		)		
	clut.mixed					<- as.numeric( colnames(clut[,clut.mixed]) )
	if(verbose) cat(paste("\nnumber of clusters with at least one in-country seq is n=", length(clut.mixed)))
	cluphy.mixed				<- merge(data.table(cluster=clut.mixed), df.cluinfo, by="cluster")
	if(verbose) cat(paste("\nnumber of seq in clusters with at least one in-country seq is n=", nrow(cluphy.mixed)))
	#
	# split mixed clusters
	#
	#exclude pairs -- split would result in singleton
	tmp							<- cluphy.mixed[, list(size=length(unique(FASTASampleCode))), by="cluster"]
	cluphy.mixed.rm				<- subset(tmp,size<=2, cluster)[,cluster]		
	cluphy.mixed				<- merge( subset(tmp,size>2, cluster), df.cluinfo, all.x=1, by="cluster" )
	if(verbose) cat(paste("\nnumber of splitting clusters is n=", length(unique(cluphy.mixed[,cluster])) ))
	if(length(intersect( unique(cluphy.mixed[,cluster]),unique(cluphy.nofrgn[,cluster]) )))
		stop("unexpected overlap between 'cluphy.mixed' and 'cluphy.nofrgn'")		
	#extract selected clusters		
	cluphy.cluidx				<- clustering[["clu.idx"]][ unique( cluphy.mixed[,cluster] ) ]					
	cluphy.split.subtrees		<- lapply(cluphy.cluidx, function(x)		extract.clade(ph, x, root.edge= 1, interactive = FALSE) 		)
	names(cluphy.split.subtrees)<- unique( cluphy.mixed[,cluster] )				
	#split selected clusters				
	#subset(cluphy.mixed, select=c(cluster,FASTASampleCode,Patient,CountryInfection))
	#splitexpr<- parse(text=paste( char.frgn, char.frgntn, sep=" | " ))
	#cluphy.df<- df.cluinfo
	cluphy.split.subtrees		<- hivc.clu.splitcluster(cluphy.split.subtrees, df.cluinfo, parse(text=paste( char.frgn, char.frgntn, sep=" | " )) )
	#z<- subset(tmp$cluphy.df[cluphy.mixed[,FASTASampleCode],], !is.na(cluster),select=c(cluster,FASTASampleCode,Patient,CountryInfection))
	#length(unique(z[,cluster]))	
	if(verbose) cat(paste("\nnumber of splitted clusters is n=", length(cluphy.split.subtrees)))
	#
	# collect all subtrees and verify that at least two patients in subtree 
	#	
	cluphy.subtrees				<- c(cluphy.split.subtrees, cluphy.nofrgn.subtrees)
	setkey(df.cluinfo, FASTASampleCode)
	tmp							<- which(sapply(cluphy.subtrees, function(x)		nrow(subset( df.cluinfo[x$tip.label,], length(unique(na.omit(Patient)))>1 ))>0			))	
	if(verbose)	cat(paste("\nnumber of final in-country clusters is n=",length(tmp)))
	cluphy.subtrees.names		<- as.numeric(names(cluphy.subtrees)[tmp])
	cluphy.subtrees				<- lapply(tmp, function(i) cluphy.subtrees[[i]]	)
	names(cluphy.subtrees)		<- cluphy.subtrees.names	
	# build new clustering data.table
	cluphy.df					<- rbindlist( lapply( seq_along(cluphy.subtrees), function(i)	data.table(FASTASampleCode=cluphy.subtrees[[i]]$tip.label, cluster=cluphy.subtrees.names[i] )		))	
	setkey(cluphy.df, FASTASampleCode)
	df.cluinfo[,cluster:=NULL]
	df.cluinfo					<- merge(cluphy.df,df.cluinfo,by="FASTASampleCode")
	if(verbose) cat(paste("\nnumber of seq in in-country clusters is n=", nrow(df.cluinfo)))
	#
	# build polyphyletic tree from clusters and construct 'clustering' for this tree
	#
	tmp										<- hivc.clu.polyphyletic.clusters(df.cluinfo, cluphy.subtrees=cluphy.subtrees, plot.file=plot.file )
	cluphy									<- tmp$cluphy
	clustering								<- hivc.clu.clusterbycluphy(cluphy.subtrees, cluphy)

	list(clustering=clustering, df.cluinfo=df.cluinfo, cluphy=cluphy, cluphy.subtrees=cluphy.subtrees)
}
######################################################################################
hivc.clu.getplot.msmexposuregroup<- function(ph, clustering, df.cluinfo, verbose=1, plot.file= NA, levels.msm=c("BI","MSM","IDU","NA"), levels.het=c("BI","HET","IDU","NA"), levels.mixed=c("BI","MSM","HET","IDU","NA"),levels.oth="OTH", split.clusters=0)
{	
	clut						<- table( df.cluinfo[,cluster,Trm], useNA= "always" )
	rownames(clut)[nrow(clut)]	<- "NA"
	clut						<- clut[,-ncol(clut)]
	clut						<- rbind(clut, apply(clut, 2, sum) )
	rownames(clut)[nrow(clut)]	<- "sum"
	#
	# clusters with MSM only
	#
	clut.onlymsm				<- apply(clut, 2, function(x)		sum( x[levels.msm] )==x["sum"]		)
	clut.onlymsm				<- as.numeric(colnames(clut[,clut.onlymsm]))
	cluphy.onlymsm				<- merge(data.table(cluster=clut.onlymsm), df.cluinfo, by="cluster")	
	if(verbose)	cat(paste("\nnumber of clusters with MSM only seq is n=",length(unique(cluphy.onlymsm[,cluster]))))		
	#
	# clusters with MSM/HET
	#
	clut.het					<- apply(clut, 2, function(x)		sum( x[levels.msm] )!=x["sum"] & sum( x[levels.het] )!=x["sum"] & sum( x[levels.mixed] )==x["sum"]		)
	clut.het					<- as.numeric(colnames(clut[,clut.het]))
	cluphy.het					<- merge(data.table(cluster=clut.het), df.cluinfo, by="cluster")
	if(verbose)	cat(paste("\nnumber of clusters with MSM/HET seq is n=",length(unique(cluphy.het[,cluster]))))
	#	
	# clusters with MSM/HET-M
	#	
	cluphy.hetM 				<- cluphy.het[,list( HetM.only= all(Sex=="M") ),by=cluster]
	if(verbose)	cat(paste("\nnumber of clusters with MSM/HET-M seq is n=",nrow(subset(cluphy.hetM,HetM.only))))
	cluphy.hetM					<- merge(subset(cluphy.hetM,HetM.only,select=cluster), df.cluinfo, all.x=1, by="cluster")					
	#
	# get subtrees corresponding to MSM only and MSM/HET-M clusters
	#
	cluphy.msm					<- rbind(cluphy.onlymsm, cluphy.hetM)
	cluphy.msm.subtrees			<- hivc.clu.polyphyletic.clusters(cluphy.msm , ph=ph, clustering=clustering)$cluphy.subtrees
	#
	# clusters with MSM/HET-F
	#
	cluphy.hetF 				<- cluphy.het[,list( with.HETF= any(Sex=="F") ),by=cluster]
	if(verbose)	cat(paste("\nnumber of clusters with MSM/HET-F seq is n=",nrow(subset(cluphy.hetF,with.HETF))))
	cluphy.hetF					<- merge(subset(cluphy.hetF,with.HETF,select=cluster), df.cluinfo, all.x=1, by="cluster")
	if(length(intersect( unique(cluphy.hetF[,cluster]),unique(cluphy.hetM[,cluster]) )))
		stop("unexpected overlap between cluphy.hetF and cluphy.hetM")
	if(length(intersect( unique(cluphy.onlymsm[,cluster]),unique(cluphy.hetM[,cluster]) )))
		stop("unexpected overlap between cluphy.onlymsm and cluphy.hetM")		
	#
	# clusters with OTH
	#	
	clut.oth							<- apply(clut, 2, function(x)		x["MSM"]>0 & sum( x[levels.msm] )!=x["sum"] & x[levels.oth]>0 & sum( x[levels.oth] )!=x["sum"] & sum( x[levels.het] )!=x["sum"] & sum( x[unique(c(levels.msm,levels.het,levels.oth))] )==x["sum"]		)
	clut.oth							<- as.numeric(colnames(clut[,clut.oth]))
	cluphy.oth							<- merge(data.table(cluster=clut.oth), df.cluinfo, by="cluster")
	if(verbose)	cat(paste("\nnumber of clusters with OTH seq is n=",length(unique(cluphy.oth[,cluster]))))
	# combine cluphy.oth and cluphy.hetF for splitting
	if(length(intersect( unique(cluphy.hetF[,cluster]),unique(cluphy.oth[,cluster]) )))
		stop("unexpected overlap between cluphy.hetF and cluphy.oth")
	cluphy.split						<- rbind(cluphy.hetF, cluphy.oth)
	# exclude pairs -- split would result in singleton
	tmp									<- cluphy.split[, list(size=length(unique(FASTASampleCode))), by="cluster"]
	cluphy.split.rm						<- subset(tmp,size<=2, cluster)[,cluster]		
	cluphy.split						<- merge( subset(tmp,size>2, cluster), df.cluinfo, all.x=1, by="cluster" )	
	# extract subtrees for MSM/HET-F clusters		
	cluphy.cluidx						<- clustering[["clu.idx"]][ unique( cluphy.split[,cluster] ) ]					
	cluphy.split.subtrees				<- lapply(cluphy.cluidx, function(x)		extract.clade(ph, x, root.edge= 1, interactive = FALSE) 		)
	names(cluphy.split.subtrees)		<- unique( cluphy.split[,cluster] )	
	if(split.clusters)	#split cluster at HETF or OTH sequence
	{
		# split HET-F, OTH clusters					
		#subset(cluphy.hetF, select=c(cluster,FASTASampleCode,Patient,Sex))
		#splitexpr<- parse(text='Sex=="F"')
		#cluphy.df<- df.cluinfo
		if(verbose)	cat(paste("\nnumber of splitting clusters is n=",length(unique(cluphy.split[,cluster]))))
		cluphy.split.subtrees			<- hivc.clu.splitcluster(cluphy.split.subtrees, df.cluinfo, parse(text='Sex=="F" | Trm=="OTH"'))
		#z<- subset(tmp$cluphy.df[cluphy.hetF[,FASTASampleCode],], !is.na(cluster),select=c(cluster,FASTASampleCode,Patient,Sex,Trm))			
		if(verbose)	cat(paste("\nnumber of retained clusters is n=",length(cluphy.split.subtrees)))		
	}
	else				#drop HETF or OTH sequence from cluster
	{		
		if(verbose)	cat(paste("\nnumber of droptip clusters is n=",length(unique(cluphy.split[,cluster]))))
		cluphy.split.subtrees			<- hivc.clu.droptipincluster(cluphy.split.subtrees, df.cluinfo, parse(text='Sex=="F" | Trm=="OTH"') )		
		if(verbose)	cat(paste("\nnumber of retained clusters is n=",length(cluphy.split.subtrees)))						
	}
	#
	# collect all subtrees and verify that at least one Trm=="MSM" and at least two patients in subtree 
	#	
	cluphy.subtrees				<- c(cluphy.split.subtrees, cluphy.msm.subtrees)
	setkey(df.cluinfo, FASTASampleCode)
	tmp							<- which(sapply(cluphy.subtrees, function(x)		nrow(subset( df.cluinfo[x$tip.label,], any(Trm=="MSM") & length(unique(na.omit(Patient)))>1 ))>0			))	
	if(verbose)	cat(paste("\nnumber of final MSM clusters is n=",length(tmp)))
	cluphy.subtrees.names		<- as.numeric(names(cluphy.subtrees)[tmp])
	cluphy.subtrees				<- lapply(tmp, function(i) cluphy.subtrees[[i]]	)
	names(cluphy.subtrees)		<- cluphy.subtrees.names	
	# build new clustering data.table
	cluphy.df					<- rbindlist( lapply( seq_along(cluphy.subtrees), function(i)	data.table(FASTASampleCode=cluphy.subtrees[[i]]$tip.label, cluster=cluphy.subtrees.names[i] )		))	
	setkey(cluphy.df, FASTASampleCode)
	df.cluinfo[,cluster:=NULL]
	df.cluinfo					<- merge(cluphy.df,df.cluinfo,by="FASTASampleCode")
	if(verbose) cat(paste("\nnumber of seq in in-country clusters is n=", nrow(df.cluinfo)))
	#
	# build polyphyletic tree from clusters and construct 'clustering' for this tree
	#
	tmp										<- hivc.clu.polyphyletic.clusters(df.cluinfo, cluphy.subtrees=cluphy.subtrees, plot.file=plot.file )
	cluphy									<- tmp$cluphy
	clustering								<- hivc.clu.clusterbycluphy(cluphy.subtrees, cluphy)
	
	list(cluphy=cluphy, cluphy.clustering=clustering, cluphy.df=df.cluinfo, cluphy.subtrees=cluphy.subtrees)
}	
######################################################################################
hivc.clu.getplot.mixedexposuregroup<- function(ph, clustering, df.cluinfo, verbose=1, plot.file= NA, levels.msm=c("BI","MSM","IDU","NA"), levels.het=c("BI","HET","IDU","NA"), levels.oth=c("OTH","NA"))
{	
	clut						<- table( df.cluinfo[,cluster,Trm], useNA= "always" )
	rownames(clut)[nrow(clut)]	<- "NA"
	clut						<- clut[,-ncol(clut)]
	clut						<- rbind(clut, apply(clut, 2, sum) )
	rownames(clut)[nrow(clut)]	<- "sum"
	clut.sametype				<- apply(clut, 2, function(x)		any( c( sum( x[levels.msm] ),sum( x[levels.het] ),sum( x[levels.oth] )	)==x["sum"])		)
	clut.onlymsm				<- apply(clut, 2, function(x)		sum( x[levels.msm] )==x["sum"]		)	
	if(verbose) print( table(clut.sametype) )
	if(verbose) print( table(clut.onlymsm) )
	#get selected clusters into single phylogeny
	cluphy.df					<- merge(data.table(cluster=which(!clut.sametype)), df.cluinfo, by="cluster")
	tmp							<- hivc.clu.polyphyletic.clusters(cluphy.df, ph=ph, clustering=clustering, plot.file=plot.file )
	
	ans<- list(cluphy=tmp$cluphy, cluphy.df=cluphy.df )
	ans
}	
######################################################################################
hivc.clu.getplot.excludeallmultifrgninfection<- function(ph, clustering, df.cluinfo, verbose=1, plot.file= NA, char.select= "!all(!is.na(CountryInfection) & CountryInfection!='NL')", pdf.scaley=25, pdf.xlim=0, cex.nodelabel=0.2, cex.tiplabel=0.2, adj.tiplabel= c(-0.15,0.5))
{		
	require(phangorn)
	expr.select	<- parse(text=char.select)
	tmp			<- subset( df.cluinfo[,eval(expr.select),by="cluster"], V1, cluster )
	clu.exlc	<- setdiff( df.cluinfo[,cluster], tmp[,cluster] )
	if(verbose) cat(paste("\nnumber of clusters with at least one patient with NL or NA CountryInfection, n=", nrow(tmp)))
	cluphy.df	<- merge(tmp, df.cluinfo, by="cluster")
	if(verbose) cat(paste("\nseq in clusters with at least one patient with NL or NA CountryInfection, n=", nrow(cluphy.df)))
	
	tmp						<- hivc.clu.polyphyletic.clusters(cluphy.df, ph=ph, clustering=clustering, plot.file=plot.file, pdf.scaley=pdf.scaley, pdf.xlim=pdf.xlim, cex.nodelabel=cex.nodelabel, cex.tiplabel=cex.tiplabel )
	cluphy					<- tmp$cluphy
	clustering				<- hivc.clu.clusterbycluphy(tmp$cluphy.subtrees, cluphy)	#produce clustering for cluphy	
	
	list(cluphy=cluphy, cluphy.df=cluphy.df, cluphy.clustering=clustering )
}
######################################################################################
hivc.clu.getplot.multifrgninfection<- function(ph, clustering, df.cluinfo, verbose=1, plot.file=NA, pdf.xlim=0.3, pdf.scaley=4, cex.nodelabel=0.4, cex.tiplabel=0.4, adj.tiplabel= c(-0.15,0.5))		
{
	tmp			<- subset(df.cluinfo, !is.na(Patient) & !is.na(CountryInfection) & CountryInfection!="NL", c(cluster, Patient))
	tmp[,key:=tmp[, paste(Patient,cluster,sep='.')]]		
	setkey(tmp,key)
	tmp			<- unique(tmp)
	cluphy.df	<- subset( tmp[, list(clu.nFrgnInfection=length(Patient)), by="cluster"], clu.nFrgnInfection>1, cluster )
	cluphy.df	<- merge( cluphy.df, df.cluinfo, by="cluster" )
	cluphy						<- NULL
	if(!is.na(plot.file))
	{
		tmp						<- hivc.clu.polyphyletic.clusters(cluphy.df, ph=ph, clustering=clustering, plot.file=plot.file, pdf.scaley=pdf.scaley, pdf.xlim=pdf.xlim, cex.nodelabel=cex.nodelabel, cex.tiplabel=cex.tiplabel )
		cluphy					<- tmp$cluphy
	}
	list(cluphy=cluphy, cluphy.df=cluphy.df )
}	
######################################################################################
hivc.clu.getplot.potentialsuperinfections<- function(ph, clustering, cluphy.df, verbose=1, plot.file=NA, pdf.xlim=0.3, pdf.scaley=4, cex.nodelabel=0.4, cex.tiplabel=0.4, adj.tiplabel= c(-0.15,0.5))
{
	require(phangorn)
	cluphy.df			<- subset(df.cluinfo, !is.na(Patient))[,list(n.clu=length(unique(na.omit(cluster))), cluster=unique(na.omit(cluster)), cluster.merged=unique(na.omit(cluster))[1]),by="Patient"]
	cluphy.df			<- subset(cluphy.df, n.clu>1)
	cluphy.merged		<- lapply( unique(cluphy.df[,Patient]), function(x)
			{
				x				<- subset(cluphy.df,Patient==x)[,cluster]
				x.mrcas			<- clu$clustering$clu.idx[x]
				x.anc			<- lapply(x.mrcas, function(z) Ancestors(ph, z))
				x.anc.jnt		<- do.call("intersect", x.anc)
				x.mrca.depth	<- sapply(x.anc, function(z) which(z %in% x.anc.jnt)[1] )
				ans				<- lapply(seq_along(x.mrcas),function(i)	extract.clade(ph, x.mrcas[i], root.edge= x.mrca.depth[i], interactive = FALSE)		)	
				ans				<- ans[[1]]+ans[[2]] 
				ans$root.edge	<- ph$edge.length[ which( ph$edge[,2]==x.anc.jnt[1] ) ]
				ans					
			})
	cluphy				<- eval(parse(text=paste('cluphy.merged[[',seq_along(cluphy.merged),']]', sep='',collapse='+')))
	#tmp					<- cluphy.df
	cluphy.df			<- merge(data.table(cluster=unique(cluphy.df[,cluster])), df.cluinfo, by="cluster")
	if(!is.na(plot.file))
	{
		cluphy.tiplabels	<- hivc.clu.get.tiplabels( cluphy, cluphy.df )
		if(verbose) cat(paste("\nwrite tree to file",plot.file))
		hivc.clu.plot(cluphy, cluphy.df[,cluster], file=plot.file, pdf.scaley=pdf.scaley, pdf.off=0, pdf.xlim= pdf.xlim, cex.nodelabel=cex.nodelabel )											
		hivc.clu.plot.tiplabels( seq_len(Ntip(cluphy)), cluphy.tiplabels$text, cluphy.tiplabels$col, cex=cex.tiplabel, adj=adj.tiplabel, add.xinch=0, add.yinch=0 )
		dev.off()
	}
	list(cluphy=cluphy, cluphy.subtrees=cluphy.merged, cluphy.df=cluphy.df)
}	
######################################################################################
hivc.clu.getplot.female2female<- function( ph, clustering, df.cluinfo, plot.file=NA )
{	
	#select clusters for HET transmission
	clut						<- table( df.cluinfo[,cluster,Trm], useNA= "always")
	rownames(clut)[nrow(clut)]	<- "NA"
	clut						<- clut[,-ncol(clut)]
	clut						<- rbind(clut, apply(clut, 2, sum) )
	rownames(clut)[nrow(clut)]	<- "sum"		
	onlyhet.clut				<- apply(clut, 2, function(x)		sum( x["HET"] )==x["sum"]		)		
	onlyhet.df					<- subset(df.cluinfo, cluster%in%which(onlyhet.clut) )
	#select clusters for gender: all F
	tmp							<- onlyhet.df[, { 
				tmp<- length(unique(Patient)); 
				list(clu.bwpat=tmp>1, clu.genderF=all(Sex=="F") ) 
			}, by="cluster"]
	onlyhet.select				<- subset( tmp, clu.bwpat & clu.genderF, select=cluster )		
	cluphy.df					<- merge( onlyhet.select, onlyhet.df, by="cluster" )		
	#select clusters for country infection: at most one out of NL origin		
	tmp							<- cluphy.df[, list(cluster=cluster[1], Patient=Patient[1], CountryInfection=CountryInfection[1]), by="Patient"]
	tmp							<- table(tmp[,cluster,CountryInfection], useNA= "always")
	tmp							<- tmp[, -ncol(tmp)]
	rownames(tmp)[nrow(tmp)]	<- "NA"
	tmp							<- rbind(tmp, apply(tmp, 2, sum))
	rownames(tmp)[nrow(tmp)]	<- "sum"
	tmp							<- as.numeric( colnames(tmp)[ tmp["NL",] + tmp["NA",] + 1 >= tmp["sum",] ] )
	cluphy.df					<- subset( cluphy.df, cluster%in%tmp )
	tmp							<- hivc.clu.polyphyletic.clusters(cluphy.df, ph=ph, clustering=clustering, plot.file=plot.file, pdf.scaley=4, cex.nodelabel=0.4, cex.tiplabel=0.4)
		
	ans<- list(cluphy=tmp$cluphy, cluphy.df=cluphy.df, cluphy.subtrees=tmp$cluphy.subtrees )
	ans
}	
######################################################################################
hivc.clu.getplot.mixedfrgngroup<- function(ph, clustering, df.cluinfo, verbose=1, plot.file= NA, char.frgn='CountryInfection=="FRGN"', char.frgntn='CountryInfection=="FRGNTN"' )							
{	
	clut						<- parse(text=paste("!is.na(CountryInfection) & (",char.frgn," | ",char.frgntn,") ", sep=''))		
	clut						<- table( df.cluinfo[,cluster, eval(clut) ] )
	clut.mixed					<- apply(clut, 2, function(x)	all(x>0)		)
	clut.mixed					<- as.numeric( colnames(clut[,clut.mixed]) )
	if(verbose) cat(paste("\nnumber of mixed clusters (in-country / frgn) is n=", length(clut.mixed) ))
	cluphy.df					<- merge(data.table(cluster=clut.mixed), df.cluinfo, by="cluster")		
	cluphy			<- NULL
	
	if(!is.na(plot.file))
	{
		tmp			<- hivc.clu.polyphyletic.clusters(cluphy.df, ph=ph, clustering=clustering, plot.file=plot.file, pdf.scaley=4, cex.nodelabel=0.4, cex.tiplabel=0.4)
		cluphy		<- tmp$cluphy
	}
	list(cluphy=cluphy, cluphy.df=cluphy.df )			
}
######################################################################################
hivc.clu.get.tiplabels<- function(ph, 	df.info, col.notmsm="#4EB3D3", col.Early="#EF9708", col.highVL="#FEE391", col.AfterTreat="#D4B9DA", col.green="#D9F0A3", col.latePres="#FA9FB5",
										select=c("CountryInfection","Trm","Sex","isAcute","lRNA.early","NegT","AnyPos_T1","PosSeqT","lRNA.hb4tr_LT","lRNA_bTS","lRNA_TS","lRNA_aTS","lRNAi_bTS","lRNAi_aTS","AnyT_T1","TrImo_bTS","TrImo_aTS","PosCD4_T1","CD4_T1","CD4_bTS","CD4_TS","CD4_aTS","Patient","RegionHospital") )
{
	require(colorspace)
	require(RColorBrewer)
	#
	#	PATIENT
	#
	#set colors CountryInfection
	tmp					<- rep("transparent",nrow(df.info))
	df.info[,CountryInfection.col:=tmp]
	set(df.info, which(df.info[,!is.na(CountryInfection) & CountryInfection!="NL"]), "CountryInfection.col", col.notmsm)
	#set colors Patient
	tmp					<- unique( df.info[,Patient] )
	#tmp2				<- diverge_hcl(length(tmp), h = c(246, 40), c = 96, l = c(85, 90))		
	tmp					<- data.table(Patient=tmp, Patient.col="transparent", key="Patient")
	df.info				<- merge(  df.info, tmp, all.x=1, by="Patient" )
	#set colors Sex
	df.info[,Sex.col:="transparent"]						
	set(df.info, which(df.info[,!is.na(Sex) & Sex!="M"]), "Sex.col", col.notmsm)
	#set colors MSM or BI	
	df.info[,Trm.col:="transparent"]						
	set(df.info, which(df.info[,!is.na(Trm) & Trm!="MSM"]), "Trm.col", col.notmsm)
	#set colors RegionHospital
	tmp					<- levels( df.info[,RegionHospital] )		
	tmp2				<- brewer.pal(length(tmp), "Dark2")
	tmp					<- data.table(RegionHospital=tmp, RegionHospital.col=tmp2, key="RegionHospital")
	df.info				<- merge(  df.info, tmp, all.x=1, by="RegionHospital" )
	#set colors isAcute
	df.info[,isAcute.col:="transparent"]						
	set(df.info, which(df.info[,isAcute%in%c("Yes","Maybe")]), "isAcute.col", col.Early)
	#
	#	TIMELINE
	#
	#	set transparent colors:		PosSeqT NegT AnyPos_T1	
	df.info[, NegT.col:="transparent"]
	df.info[, AnyPos_T1.col:="transparent"]
	#
	#	TREATMENT
	#
	#	set AnyT_T1.col				if 	PosSeqT<AnyT_T1 	orange 		else pink
	df.info[, AnyT_T1.col:="transparent"]
	select.seqb4tr		<- which(df.info[, !is.na(PosSeqT) & !is.na(AnyT_T1) & PosSeqT<=AnyT_T1])
	select.seqatr		<- which(df.info[, !is.na(PosSeqT) & !is.na(AnyT_T1) & PosSeqT>AnyT_T1])
	set(df.info, select.seqb4tr, "AnyT_T1.col", col.Early)
	set(df.info, select.seqatr, "AnyT_T1.col", col.AfterTreat)
	df.info[, TrImo_bTS.col:="transparent"]
	df.info[, TrImo_aTS.col:="transparent"]
	set(df.info, which(df.info[, TrImo_aTS>4]), "TrImo_aTS.col", col.AfterTreat)
	set(df.info, which(df.info[, TrImo_bTS>4]), "TrImo_bTS.col", col.AfterTreat)
	df.info[, PosSeqT.col:="transparent"]
	set(df.info, select.seqb4tr, "PosSeqT.col", col.Early)
	set(df.info, select.seqatr, "PosSeqT.col", col.AfterTreat)	
	#
	#	VIRAL LOAD
	#	
	df.info[, lRNA_bTS.col:="transparent"]
	set(df.info,	which(df.info[,lRNA_bTS>5]),	"lRNA_bTS.col",		col.highVL) 
	df.info[, lRNA_TS.col:="transparent"]
	set(df.info,	which(df.info[,lRNA_TS>5]),	"lRNA_TS.col",		col.highVL) 
	df.info[, lRNA_aTS.col:="transparent"]
	set(df.info,	which(df.info[,lRNA_aTS>3.5]),	"lRNA_aTS.col",		col.highVL) 
	df.info[, lRNAi_bTS.col:="transparent"]
	set(df.info,	which(df.info[,lRNAi_bTS>0.75]),	"lRNAi_bTS.col",		col.highVL) 
	df.info[, lRNAi_aTS.col:="transparent"]
	set(df.info,	which(df.info[,lRNAi_aTS>0.25]),	"lRNAi_aTS.col",		col.highVL) 
	df.info[, lRNA.hb4tr_LT.col:="transparent"]
	set(df.info,	which(df.info[,!is.na(lRNA.hb4tr_LT) & PosSeqT<=lRNA.hb4tr_LT]),	"lRNA.hb4tr_LT.col",		col.highVL) 
	df.info[, lRNA.early.col:="transparent"]
	set(df.info,	which(df.info[, lRNA.early]),	"lRNA.early.col",		col.Early) 
	#
	#	CD4
	#		
	df.info[, CD4_T1.col:="transparent"]
	set(df.info,	which(df.info[,CD4_T1>350]),	"CD4_T1.col",		col.green)	
	df.info[, PosCD4_T1.col:="transparent"]
	set(df.info,	which(df.info[,CD4_T1>350]),	"PosCD4_T1.col",		col.green)	
	df.info[, CD4_bTS.col:="transparent"]
	set(df.info,	which(df.info[,CD4_bTS>350]),	"CD4_bTS.col",		col.green)	
	df.info[, CD4_TS.col:="transparent"]
	set(df.info,	which(df.info[,CD4_TS>350]),	"CD4_TS.col",		col.green)	
	df.info[, CD4_aTS.col:="transparent"]
	set(df.info,	which(df.info[,CD4_aTS>350]),	"CD4_aTS.col",		col.green)
	#
	#	color late presenter
	#
	tmp<- which(df.info[,CD4_T1<350])
	set(df.info,	tmp,	"CD4_T1.col",			col.latePres)
	set(df.info,	tmp,	"PosCD4_T1.col",		col.latePres)
	set(df.info,	tmp,	"lRNA_bTS.col",			col.latePres)
	set(df.info,	tmp,	"lRNAi_bTS.col",		col.latePres)
	#	
	#	convert time to string	& handle inaccurate NegT or AnyPos_T1
	#
	tmp		<- which(df.info[, as.POSIXlt(NegT)$mday==1 & as.POSIXlt(NegT)$mon==0])				#since NegT were reset, these are the ones with inaccurate month
	set(df.info,NULL,	"NegT", 		as.character( df.info[,NegT], "%y.%m" ))
	set(df.info, tmp,	"NegT", 		paste(df.info[tmp,substr(NegT,1,3)],'??',sep='') )
	tmp		<- which(df.info[, as.POSIXlt(AnyPos_T1)$mday==31 & as.POSIXlt(AnyPos_T1)$mon==11])	#since AnyPos_T1 were reset, these are the ones with inaccurate month
	set(df.info,NULL,	"AnyPos_T1",	as.character( df.info[,AnyPos_T1], "%y.%m" ))
	set(df.info, tmp, 	"AnyPos_T1", 	paste(df.info[tmp,substr(AnyPos_T1,1,3)],'??',sep='') )
	set(df.info,NULL,	"PosSeqT", 		as.character( df.info[,PosSeqT], "%y.%m" ))
	set(df.info,NULL,	"lRNA.hb4tr_LT",as.character( df.info[,lRNA.hb4tr_LT], "%y.%m" ))
	set(df.info,NULL,	"PosCD4_T1", 	as.character( df.info[,PosCD4_T1], "%y.%m" ))
	set(df.info,NULL,	"AnyT_T1", 		as.character( df.info[,AnyT_T1], "%y.%m" ))
	#	
	#	set isAcute to either Y or M or N
	set(df.info,NULL,"isAcute", 		as.character( df.info[,isAcute]))
	set(df.info,which(df.info[,isAcute=="No"]),"isAcute",'N')
	set(df.info,which(df.info[,isAcute=="Maybe"]),"isAcute",'M')
	set(df.info,which(df.info[,isAcute=="Yes"]),"isAcute",'Y')
	#	set lRNA.early to either HVLE or ----
	set(df.info,NULL,"lRNA.early", 		as.character( df.info[, lRNA.early]))
	set(df.info,which(df.info[,lRNA.early=="TRUE"]),"lRNA.early","Y")
	set(df.info,which(df.info[,lRNA.early=="FALSE"]),"lRNA.early","N")
	#	set lRNAi_bTS lRNAi_aTS
	set(df.info,NULL,"lRNAi_bTS", 		as.character( round(df.info[, lRNAi_bTS],d=2)))
	set(df.info,NULL,"lRNAi_aTS", 		as.character( round(df.info[, lRNAi_aTS],d=2)))
	#	set TrImo_bTS TrImo_aTS
	set(df.info,NULL,"TrImo_bTS", 		as.character( round(df.info[, TrImo_bTS],d=1)))
	set(df.info,NULL,"TrImo_aTS", 		as.character( round(df.info[, TrImo_aTS],d=1)))
	#	set CD4_T1
	set(df.info,NULL,"CD4_T1", 		as.character( round(df.info[, CD4_T1],d=0)))
	#	set CD4_TS CD4_bTS CD4_aTS
	set(df.info,NULL,"CD4_TS", 		as.character( df.info[, CD4_TS]))
	set(df.info,NULL,"CD4_bTS", 		as.character( df.info[, CD4_bTS]))
	set(df.info,NULL,"CD4_aTS", 		as.character( df.info[, CD4_aTS]))
	#	set 'Amst' to 'A'
	set(df.info,NULL,"RegionHospital", 		as.character( df.info[,RegionHospital]))
	set(df.info,which(df.info[,RegionHospital=="Amst"]),"RegionHospital",'A')	
	#
	#	handle missing entries -- ensure that alignment is OK
	#
	setkey(df.info, Patient)
	tmp					<- which(is.na(df.info[,Patient]))
	set(df.info, tmp, "Patient", '')
	set(df.info, tmp, "Patient.col", "transparent")
	tmp					<- which(is.na(df.info[,RegionHospital]))
	set(df.info, tmp, "RegionHospital", '-')
	set(df.info, tmp, "RegionHospital.col", "transparent")		
	tmp					<- which(is.na(df.info[,CountryInfection]))
	set(df.info, tmp, "CountryInfection", "--")
	set(df.info, tmp, "CountryInfection.col", "transparent")				
	tmp					<- which(is.na(df.info[,Trm]))
	set(df.info, tmp, "Trm", '--')
	set(df.info, tmp, "Trm.col", "transparent")
	tmp					<- which(is.na(df.info[,Sex]))
	set(df.info, tmp, "Sex", '-')
	set(df.info, tmp, "Sex.col", "transparent")	
	set(df.info, which(df.info[,is.na(isAcute)]), "isAcute", 	'-')	
	tmp					<- which(is.na(df.info[,AnyPos_T1]))
	set(df.info, tmp, "AnyPos_T1", "--.--")	
	set(df.info, tmp, "AnyPos_T1.col", "transparent")
	tmp					<- which(is.na(df.info[,NegT]))
	set(df.info, tmp, "NegT", "--.--")	
	set(df.info, tmp, "NegT.col", "transparent")
	tmp					<- which(is.na(df.info[,PosSeqT]))
	set(df.info, tmp, "PosSeqT", "--.--")	
	set(df.info, tmp, "PosSeqT.col", "transparent")		
	tmp					<- which(is.na(df.info[,lRNA.hb4tr_LT]))
	set(df.info, tmp, "lRNA.hb4tr_LT", "--.--")		
	tmp					<- which(is.na(df.info[,PosCD4_T1]))
	set(df.info, tmp, "PosCD4_T1", "--.--")		
	tmp					<- which(is.na(df.info[,AnyT_T1]))
	set(df.info, tmp, "AnyT_T1", "--.--")		
	set(df.info, which(df.info[,is.na(TrImo_aTS)]), "TrImo_aTS", 	'-')
	set(df.info, which(df.info[,is.na(TrImo_bTS)]), "TrImo_bTS", 	'-')	
	set(df.info, which(df.info[,is.na(CD4_T1)]), "CD4_T1", 			'---')
	set(df.info, which(df.info[,is.na(PosCD4_T1)]), "PosCD4_T1", 	'---')
	set(df.info, which(df.info[,is.na(CD4_TS)]), "CD4_TS", 			'---')
	set(df.info, which(df.info[,is.na(CD4_bTS)]), "CD4_bTS", 		'---')
	set(df.info, which(df.info[,is.na(CD4_aTS)]), "CD4_aTS", 		'---')
	#
	#	add suffixes 
	#	
	set(df.info,NULL,"AnyPos_T1",		paste("d:",df.info[,AnyPos_T1],sep=''))
	set(df.info,NULL,"NegT", 			paste("n:",df.info[,NegT],sep=''))
	set(df.info,NULL,"PosSeqT", 		paste("s:",df.info[,PosSeqT],"   ",sep=''))
	set(df.info,NULL,"lRNA.hb4tr_LT",	paste("VLeh:",df.info[,lRNA.hb4tr_LT],sep=''))
	set(df.info,NULL,"AnyT_T1", 		paste("TR+:",df.info[,AnyT_T1],sep=''))
	set(df.info,NULL,"TrImo_bTS", 		paste("TR-:",df.info[,TrImo_bTS],sep=''))
	set(df.info,NULL,"TrImo_aTS", 		paste(":",df.info[,TrImo_aTS],"   ",sep=''))			
	set(df.info,NULL,"isAcute", 		paste("AC:",df.info[,isAcute],sep=''))
	set(df.info,NULL,"lRNA.early", 		paste(":",df.info[,lRNA.early],"   ",sep=''))
	set(df.info,NULL,"CountryInfection",paste("inf:",df.info[,CountryInfection],sep=''))
	set(df.info,NULL,"RegionHospital",	paste("",df.info[,RegionHospital],sep=''))
	set(df.info,NULL,"lRNA_bTS", 		paste("VL:",df.info[,lRNA_bTS],sep=''))
	set(df.info,NULL,"lRNA_TS", 		paste(":",df.info[,lRNA_TS],sep=''))
	set(df.info,NULL,"lRNA_aTS", 		paste(":",df.info[,lRNA_aTS],sep=''))
	set(df.info,NULL,"lRNAi_bTS", 		paste("VI:",df.info[,lRNAi_bTS],sep=''))
	set(df.info,NULL,"lRNAi_aTS", 		paste(":",df.info[,lRNAi_aTS],"   ",sep=''))	
	set(df.info,NULL,"PosCD4_T1", 		paste("CD4:",df.info[,PosCD4_T1],sep=''))
	set(df.info,NULL,"CD4_T1", 			paste(":",df.info[,CD4_T1],sep=''))	
	set(df.info,NULL,"CD4_bTS", 		paste(":",df.info[,CD4_bTS],sep=''))
	set(df.info,NULL,"CD4_TS", 			paste(":",df.info[,CD4_TS],sep=''))
	set(df.info,NULL,"CD4_aTS", 		paste(":",df.info[,CD4_aTS],"   ",sep=''))
	#
	#get df.info into order of tips
	#
	setkey(df.info, FASTASampleCode)
	df.info				<- df.info[ph$tip.label,]
	#select text and col matrix 
	text				<- t( as.matrix( df.info[,select, with=0] ) )
	colnames(text)		<- ph$tip.label
	col					<- t( as.matrix( df.info[,paste(select,".col",sep=''),with=0] ) )
	colnames(col)		<- ph$tip.label
	ans<- list(text=text, col=col)
	ans
}
######################################################################################
#prepare standard format of tip labels -- requires df.info to be sorted along the tips as they appear in a phylogeny
hivc.clu.get.tiplabels.v1<- function(ph, df.info)
{
	require(colorspace)
	require(RColorBrewer)
	#set colors CountryInfection
	tmp					<- rep("transparent",nrow(df.info))
	df.info[,CountryInfection.col:=tmp]
	set(df.info, which(df.info[,CountryInfection=="NL"]), "CountryInfection.col", "#EF9708")
	#set colors Patient
	tmp					<- unique( df.info[,Patient] )
	tmp2				<- diverge_hcl(length(tmp), h = c(246, 40), c = 96, l = c(85, 90))		
	tmp					<- data.table(Patient=tmp, Patient.col=tmp2, key="Patient")
	df.info				<- merge(  df.info, tmp, all.x=1, by="Patient" )
	#set colors Sex
	tmp					<- rep("transparent",nrow(df.info))
	df.info[,Sex.col:=tmp]						
	set(df.info, which(df.info[,Sex=="M"]), "Sex.col", "#EF9708")
	#set colors MSM or BI
	tmp					<- rep("transparent",nrow(df.info))
	df.info[,Trm.col:=tmp]						
	set(df.info, which(df.info[,Trm%in%c("MSM","BI")]), "Trm.col", "#EF9708")
	#set colors RegionHospital
	tmp					<- levels( df.info[,RegionHospital] )		
	tmp2				<- brewer.pal(length(tmp), "Dark2")
	tmp					<- data.table(RegionHospital=tmp, RegionHospital.col=tmp2, key="RegionHospital")
	df.info				<- merge(  df.info, tmp, all.x=1, by="RegionHospital" )		
	#set colors time		
	tmp					<- range( range(df.info[, NegT],na.rm=1), range(df.info[, AnyPos_T1],na.rm=1) )
	tmp					<- as.POSIXlt( seq.Date(tmp[1],tmp[2]+365,by="years") )$year
	tmp2				<- heat_hcl(length(tmp), h = c(0, -100), l = c(75, 40), c = c(40, 80), power = 1)
	yearcols			<- data.table(Year=tmp, Year.col=tmp2)
	#set colors PosSeqT
	tmp					<- data.table(PosSeqT= unique( df.info[,PosSeqT] ), key="PosSeqT" )
	tmp2				<- tmp[, as.POSIXlt(PosSeqT)$year]
	tmp[,"Year":=tmp2]
	tmp					<- subset( merge(tmp, yearcols, all.x=1, by="Year"), select=c(PosSeqT,Year.col))
	setnames(tmp,"Year.col","PosSeqT.col")
	df.info			<- merge(  df.info, tmp, all.x=1, by="PosSeqT" )
	#set colors NegT
	tmp					<- data.table(NegT= unique( df.info[,NegT] ), key="NegT" )
	tmp2				<- tmp[, as.POSIXlt(NegT)$year]
	tmp[,"Year":=tmp2]
	tmp					<- subset( merge(tmp, yearcols, all.x=1, by="Year"), select=c(NegT,Year.col))
	setnames(tmp,"Year.col","NegT.col")
	df.info			<- merge(  df.info, tmp, all.x=1, by="NegT" )
	#set colors AnyPos_T1
	tmp					<- data.table(AnyPos_T1= unique( df.info[,AnyPos_T1] ), key="AnyPos_T1" )
	tmp2				<- tmp[, as.POSIXlt(AnyPos_T1)$year]
	tmp[,"Year":=tmp2]
	tmp					<- subset( merge(tmp, yearcols, all.x=1, by="Year"), select=c(AnyPos_T1,Year.col))
	setnames(tmp,"Year.col","AnyPos_T1.col")
	df.info			<- merge(  df.info, tmp, all.x=1, by="AnyPos_T1" )									
	
	#convert time to string
	set(df.info,NULL,"AnyPos_T1",substr(as.character( df.info[,AnyPos_T1] ),1,7))
	set(df.info,NULL,"NegT",substr(as.character( df.info[,NegT] ),1,7))
	set(df.info,NULL,"PosSeqT",substr(as.character( df.info[,PosSeqT] ),1,7))
	
	#handle missing entries
	setkey(df.info, Patient)
	tmp					<- which(is.na(df.info[,Patient]))
	set(df.info, tmp, "Patient", '')
	set(df.info, tmp, "Patient.col", "transparent")
	tmp					<- which(is.na(df.info[,RegionHospital]))
	set(df.info, tmp, "RegionHospital", '-')
	set(df.info, tmp, "RegionHospital.col", "transparent")		
	tmp					<- which(is.na(df.info[,CountryInfection]))
	set(df.info, tmp, "CountryInfection", "--")
	set(df.info, tmp, "CountryInfection.col", "transparent")				
	tmp					<- which(is.na(df.info[,Trm]))
	set(df.info, tmp, "Trm", '')
	set(df.info, tmp, "Trm.col", "transparent")
	tmp					<- which(is.na(df.info[,Sex]))
	set(df.info, tmp, "Sex", '')
	set(df.info, tmp, "Sex.col", "transparent")
	tmp					<- which(is.na(df.info[,AnyPos_T1]))
	set(df.info, tmp, "AnyPos_T1", "-------")
	set(df.info, NULL, "AnyPos_T1", paste("HIV+:",df.info[,AnyPos_T1],sep=''))
	set(df.info, tmp, "AnyPos_T1.col", "transparent")
	tmp					<- which(is.na(df.info[,NegT]))
	set(df.info, tmp, "NegT", "-------")
	set(df.info, NULL, "NegT", paste("HIV-:",df.info[,NegT],sep=''))
	set(df.info, tmp, "NegT.col", "transparent")
	tmp					<- which(is.na(df.info[,PosSeqT]))
	set(df.info, tmp, "PosSeqT", "-------")
	set(df.info, NULL, "PosSeqT", paste("HIVS:",df.info[,PosSeqT],sep=''))
	set(df.info, tmp, "PosSeqT.col", "transparent")		
	#get df.info into order of tips
	setkey(df.info, FASTASampleCode)
	df.info				<- df.info[ph$tip.label,]
	#select text and col matrix 
	text				<- t( as.matrix( subset(df.info,select=c(CountryInfection, Trm, Sex, NegT, AnyPos_T1, PosSeqT, Patient, RegionHospital)) ) )
	colnames(text)		<- ph$tip.label
	col					<- t( as.matrix( subset(df.info,select=c(CountryInfection.col, Trm.col, Sex.col, NegT.col, AnyPos_T1.col, PosSeqT.col, Patient.col, RegionHospital.col)) ) )
	colnames(col)		<- ph$tip.label
	ans<- list(text=text, col=col)
	ans
}	
######################################################################################
hivc.clu.plot.tiplabels<- function (tip, text, col, xx=NULL, adj = c(-0.05, 0.5), cex=1, add.xinch= 0.03, add.yinch= 0.02) 
{		
	lastPP 			<- get("last_plot.phylo", envir = .PlotPhyloEnv)
	if(is.null(xx))
		xx 			<- lastPP$xx[tip]
	yy 				<- lastPP$yy[tip]
	if(length(tip)==1)
	{
		#assume vector of text and col of same length
		if(!is.vector(text) || !is.vector(col))	stop("expect vector text and vector col")
		if(length(text)!=length(col))			stop("expect same length of text and col")
		wh				<- sapply(text, function(x){	 c(xinch(strwidth(x, units = "inches", cex = cex)), yinch(strheight(x, units = "inches", cex = cex))) 		})
		wh.total		<- c(sum(wh[1,]),max(wh[2,]))
		coord			<- matrix(NA,6,length(text),dimnames=list(c("xl","xr","yb","yt","xx","yy"),c()))
		coord["xl",]	<- xx - wh.total[1] * adj[1] - xinch(add.xinch)
		tmp				<- cumsum(wh[1,]) + xinch(add.xinch) / ncol(wh)
		coord["xl",]	<- coord["xl",] + c(0, tmp[-ncol(coord)])
		coord["xr",]	<- coord["xl",] + wh[1,] + xinch(add.xinch)/ ncol(wh)
		coord["yb",]	<- yy - wh.total[2] * adj[2] - yinch(add.yinch)
		coord["yt",]	<- coord["yb",] + wh.total[2] + yinch(add.yinch)
		coord["xx",]	<- coord["xl",] + c( diff( coord[1,] ), diff(coord[1:2,ncol(wh)])  ) / 2
		coord["yy",]	<- rep(yy, ncol(coord))
#print(wh); print(wh.total); print(coord)		
	}
	else
	{
		#assume matrix of text and col of same ncol
		if(!is.matrix(text) || !is.matrix(col))				stop("expect matrix text and vector col")
		if(ncol(text)!=ncol(col) || nrow(text)!=nrow(col))	stop("expect same dimensions of text and col")
		if(ncol(text)!=length(tip))							stop("expect cols of text and col to correspond to tips")
		coord	<- lapply( seq_along(xx),function(i)
				{
					wh				<- sapply(text[,i], function(x){	 c(xinch(strwidth(x, units = "inches", cex = cex)), yinch(strheight(x, units = "inches", cex = cex))) 		})
					wh.total		<- c(sum(wh[1,]),max(wh[2,]))
					coord			<- matrix(NA,6,nrow(text),dimnames=list(c("xl","xr","yb","yt","xx","yy"),c()))
					coord["xl",]	<- xx[i] - wh.total[1] * adj[1] - xinch(add.xinch)
					tmp				<- cumsum(wh[1,]) + xinch(add.xinch) / ncol(wh)
					coord["xl",]	<- coord["xl",] + c(0, tmp[-ncol(coord)])
					coord["xr",]	<- coord["xl",] + wh[1,] + xinch(add.xinch)/ ncol(wh)
					coord["yb",]	<- yy[i] - wh.total[2] * adj[2] - yinch(add.yinch)
					coord["yt",]	<- coord["yb",] + wh.total[2] + yinch(add.yinch)
					coord["xx",]	<- coord["xl",] + c( diff( coord[1,] ), diff(coord[1:2,ncol(wh)])  ) / 2
					coord["yy",]	<- rep(yy[i], ncol(coord))
					coord
				})
		coord	<- do.call("cbind",coord)
		text	<- as.vector(text)
		col		<- as.vector(col)
#print(coords); print(text); print(col)
	}
	rect(coord["xl",], coord["yb",], coord["xr",], coord["yt",], col = col, border=NA)
	text(coord["xx",], coord["yy",], text, cex=cex)
}
######################################################################################
hivc.clu.plot<- function(	ph, clu, edge.col.basic="black", show.node.label= T, show.tip.label=F, file=NULL,  
							highlight.edge.of.tiplabel=NULL, highlight.edge.of.tiplabel.col="red", 
							highlight.cluster=NULL, highlight.cluster.col="red",							
							pdf.scaley=10, pdf.width= 7, pdf.height=pdf.scaley*7, pdf.off=1, pdf.xlim=NULL,
							cex.nodelabel=0.5, cex.edge.incluster=1, cex.edge.outcluster= cex.edge.incluster/3, no.margin=T, ...)
{
	require(colorspace)	
	clu.edge							<- clu[ ph$edge[,1] ]
	clu.edge							<- clu.edge+1				#set col index for non-clustering edges to 1
	clu.edge[is.na(clu.edge)]			<- 1
	#cols.n								<- length(unique(clu.edge))-1	
	#cols								<- c("black",rainbow_hcl(cols.n, start = 30, end = 300))	
	clu.edge.col						<- rep(edge.col.basic, nrow(ph$edge))	#cols[clu.edge]
	clu.edge.col[clu.edge==1]			<- "grey50" 
	if(!is.null(highlight.cluster))
	{
		if(length(highlight.cluster.col)==1)
			highlight.cluster.col<- rep(highlight.cluster.col,length(highlight.cluster))
		for(i in seq_along(highlight.cluster))
			clu.edge.col[clu.edge%in%(highlight.cluster[[i]]+1)]<- 	highlight.cluster.col[i]	
	}	
	if(!is.null(highlight.edge.of.tiplabel))
	{
		highlight.edge<- lapply(highlight.edge.of.tiplabel, function(x)		which( ph$edge[,2]%in%which( substr(ph$tip.label, 1, nchar(x))==x ) )		)					
		if(length(highlight.edge.of.tiplabel.col)==1)
			highlight.edge.of.tiplabel.col	<- rep(highlight.edge, length(highlight.edge.of.tiplabel.col) )		
		for(i in seq_along(highlight.edge))
			clu.edge.col[highlight.edge[[i]]]		<- highlight.edge.of.tiplabel.col[i]		
	}	
	clu.edge.width						<- rep(cex.edge.outcluster, length(clu.edge))
	clu.edge.width[clu.edge!=1]			<- cex.edge.incluster
	clu.edge.lty						<- rep(1, length(clu.edge))
	clu.edge.lty[clu.edge!=1]			<- 1
	if(class(file)=="character")
		pdf(file,width=pdf.width,height=pdf.height)
	if(no.margin) 	par(mar=c(0,0,0,0))
	plot.coordinates					<- plot(ph, show.tip.label=show.tip.label, show.node.label=show.node.label, cex=cex.nodelabel, edge.color=clu.edge.col, edge.width=clu.edge.width, edge.lty=clu.edge.lty, x.lim=pdf.xlim, ...)
	if(class(file)=="character" && pdf.off)
		dev.off()
	plot.coordinates
}
######################################################################################
hivc.clu.plot.tptn<- function(stat.cmp.x, stat.cmp.y, plotfile, cols, xlab= "#FP (among all)", ylab= "%TP (among all)", xlim= range(stat.cmp.x), labels= as.numeric(colnames(stat.cmp.x)), verbose=1)
{		
	if(verbose) cat(paste("\nplot file",plotfile))
	pdf(file=plotfile,width=5,height=5)
	plot(1,1,type='n',bty='n',xlim=xlim,ylim=range(stat.cmp.y),xlab=xlab,ylab=ylab)
	dummy	<- sapply(seq_len(nrow(stat.cmp.x)),function(i)
			{
				points(stat.cmp.x[i,],stat.cmp.y[i,],col=cols[i],type='b')
				text(stat.cmp.x[i,],stat.cmp.y[i,],col=cols[i],labels=labels,adj=c(-0.8,0.5),cex=0.5)
			})
	legend("bottomright",border=NA,bty='n',fill=cols,legend=names(cols))
	dev.off()
}	
######################################################################################
hivc.clu.plot.withinpatientseq.not.samecluster<- function(missed.ph, clustering, clusters.tp, missed.df, plotfile, verbose=1)
{
	clusters.tp.missed									<- merge( data.table(Patient=unique( clusters.tp$clu.missedtp[,Patient] ) ), subset(missed.df,select=c(Patient, FASTASampleCode, cluster, Node)), by="Patient" )				
	# highlight sequences of patients that are not clustering in blue
	node.missed											<- subset(clusters.tp.missed,is.na(cluster),Node)[,Node]
	missed.ph$tip.label[ node.missed ]					<- paste("TPm",missed.ph$tip.label[ node.missed ],sep='-') 
	tmp													<- which( missed.df[,Node]%in%node.missed )
	set(missed.df, tmp, "FASTASampleCode", paste("TPm",missed.df[tmp,FASTASampleCode],sep='-') )		
	# highlight sequences of patients that are not in the same cluster in red
	node.wrong											<- data.table(Patient=unique(clusters.tp$clu.diffclu[,Patient]))
	node.wrong											<- merge(node.wrong , subset(missed.df,select=c(Patient, FASTASampleCode, cluster, Node)), by="Patient" )
	node.wrong											<- subset(node.wrong, !is.na(cluster), Node)[,Node]
	missed.ph$tip.label[ node.wrong ]					<- paste("TPw",missed.ph$tip.label[ node.wrong ],sep='-') 
	tmp													<- which( missed.df[,Node]%in%node.wrong )
	set(missed.df, tmp, "FASTASampleCode", paste("TPw",missed.df[tmp,FASTASampleCode],sep='-') )
	# highlight clusters with missed within patient seq 
	missed.clumem										<- clustering[["clu.mem"]]
	tmp													<- na.omit( unique(clusters.tp.missed[, cluster]) )
	missed.clumem[ !missed.clumem%in%tmp ]	 			<- NA
	#plot tree to file	
	tiplabels.missed									<- hivc.clu.get.tiplabels( missed.ph, missed.df )					
	if(verbose) cat(paste("write tree to file",plotfile))
	hivc.clu.plot(missed.ph, missed.clumem, file=plotfile, pdf.scaley=25, pdf.off=0, highlight.edge.of.tiplabel=c("TPw","TPm"), highlight.edge.of.tiplabel.col= c("red","blue"), cex.nodelabel=0.1 )
	hivc.clu.plot.tiplabels( seq_len(Ntip(missed.ph)), tiplabels.missed$text, tiplabels.missed$col, cex=0.12, adj=c(-0.15,0.5), add.xinch=0, add.yinch=0 )
	dev.off()
}	
######################################################################################
hivc.clu.mrca<- function(ph, tiplabel, x.tip=NULL)
{
	require(phangorn)
	if(is.null(x.tip))
		x.tip	<- match(tiplabel, ph$tip.label)				
	x.tip.anc	<- lapply(x.tip, function(z) Ancestors(ph, z) )
	x.tip.jnt	<- my.intersect.n( x.tip.anc )			
	x.tip.anc[[1]][min(which(x.tip.anc[[1]] %in% x.tip.jnt))]				
}	
######################################################################################
hivc.clu.min.transmission.cascade<- function(brlmat)
{
	if(is.matrix(brlmat))
		ans	<- .Call("hivc_clu_mintransmissioncascade", brlmat[upper.tri(brlmat)])  						
	else
		ans	<- .Call("hivc_clu_mintransmissioncascade", brlmat)
	ans		
}
######################################################################################
#' Compute a statistic of patristic distances for the subtree starting at \code{node}
#' @param node 				internal node at which the subtree starts
#' @param tree 				phylogenetic tree
#' @param distmat 			matrix of patristic distances, either precomputed or NULL
#' @param eval.dist.btw 	string, either "leaf" or "all". Specifies the nodes of the subtree on which patristic distances are considered.
#' @param stat.fun 			function, statistic to be computed from a set of patristic distances.
#' @return statistic of patristic distances		
hivc.clu.brdist.stats.subtree<- function(node, tree, distmat, eval.dist.btw="leaf", stat.fun=max)
{
	require(ape)
	require(igraph)
	require(geiger)	
	if(eval.dist.btw=="leaf")
	{
		nlist	<- tips(tree,node)
		foo 	<- distmat[nlist,nlist] 		
	}
	else if(eval.dist.btw=="all")
	{
		nlist	<- tips(tree,node)
		elist 	<- tree$edge[which.edge(tree,nlist),2]
		foo 	<- distmat[elist,elist] 	
	}
	else	
		stop("unrecognized eval.dist.btw")
	return( stat.fun(foo[upper.tri(foo,diag=FALSE)]) )
}
######################################################################################
#' Compute a statistic of patristic distances for each subtree in the tree
hivc.clu.brdist.stats<- function(tree, distmat=NULL, eval.dist.btw="leaf", stat.fun=max)
{
	require(phytools)
	if(is.null(distmat))
	{
		if(eval.dist.btw=='leaf')
			distmat	<- cophenetic.phylo(tree)
		else if(eval.dist.btw=='all')
			distmat <- dist.nodes(tree)
		else	stop("unrecognized eval.dist.btw")
	}
	ntips			<- Ntip(tree)
	nint 			<- tree$Nnode 		## number of internal nodes
	return(sapply(seq.int(ntips+1,ntips+nint), function(x) 		hivc.clu.brdist.stats.subtree(x,tree,distmat,eval.dist.btw=eval.dist.btw, stat.fun=stat.fun)		))	
}
######################################################################################
hivc.clu.droptipincluster.old<- function(ph, clustering, cluphy.subtrees, df.cluinfo, splitexpr)
{
	require(phangorn)
	setkey(df.cluinfo, FASTASampleCode)
	for(i in seq_along(cluphy.subtrees))
	{
		x																	<- cluphy.subtrees[[i]]
		x.cluster															<- as.numeric(names(cluphy.subtrees)[i])
		#drop selected tips from cluster					
		x.droplabel															<- subset( df.cluinfo[x$tip.label,], eval(splitexpr), FASTASampleCode )[,FASTASampleCode]
		if(length(x.droplabel)+1==length(x$tip.label))
		{			
			x.droplabel														<- x$tip.label
			x.tip.mrca														<- NA
			x																<- numeric()					#cannot do NULL which would automatically shorten 'cluphy.subtrees' and give an outofbounds error			
		}
		else	
		{
			x																<- drop.tip(x, x.droplabel, root.edge=1)			
			x.tip															<- sapply(x$tip.label, function(z) match(z, ph$tip.label) )				
			x.tip.anc														<- lapply(x.tip, function(z) Ancestors(ph, z) )
			x.tip.jnt														<- my.intersect.n( x.tip.anc )			
			x.tip.mrca														<- x.tip.anc[[1]][min(which(x.tip.anc[[1]] %in% x.tip.jnt))]							
		}
		#drop selected tips from data.table
		set(df.cluinfo, which(df.cluinfo[,FASTASampleCode%in%x.droplabel ]), "cluster", NA)												
		#drop selected tips from clustering
		clustering[["clu.mem"]][ clustering[["clu.mem"]]==x.cluster ]		<- NA				
		if(!is.na(x.tip.mrca))		
			clustering[["clu.mem"]][ c(x.tip.mrca, x.tip) ]					<- x.cluster					#WARNING: might miss inner nodes, so size.all would not be accurate
		clustering[["clu.idx"]][ x.cluster ]								<- x.tip.mrca		
		cluphy.subtrees[[i]]												<- x							
	}	
	cluphy.subtrees	<- lapply( which(sapply(cluphy.subtrees, length)>0), function(i) 	cluphy.subtrees[[i]] 	)
	
	list(cluphy.subtrees=cluphy.subtrees, clustering=clustering, cluphy.df=df.cluinfo)
}
######################################################################################
hivc.clu.droptipincluster<- function(cluphy.subtrees, df.cluinfo, splitexpr)
{
	require(phangorn)
	setkey(df.cluinfo, FASTASampleCode)
	for(i in seq_along(cluphy.subtrees))
	{
		x																	<- cluphy.subtrees[[i]]
		x.cluster															<- as.numeric(names(cluphy.subtrees)[i])
		#drop selected tips from cluster					
		x.droplabel															<- subset( df.cluinfo[x$tip.label,], eval(splitexpr), FASTASampleCode )[,FASTASampleCode]
		if(length(x.droplabel)+1==length(x$tip.label))
			x																<- numeric()					#cannot do NULL which would automatically shorten 'cluphy.subtrees' and give an outofbounds error					
		else	
			x																<- drop.tip(x, x.droplabel, root.edge=1)					
		cluphy.subtrees[[i]]												<- x							
	}	
	cluphy.subtrees	<- lapply( which(sapply(cluphy.subtrees, length)>0), function(i) 	cluphy.subtrees[[i]] 	)	
	cluphy.subtrees
}
######################################################################################
#modify clustering by splitting existing clusters based on 'splitexpr'. This is useful when additional meta-information is available.
hivc.clu.splitcluster.old<- function(ph, clustering, cluphy.subtrees, cluphy.df, splitexpr)
{
	require(phangorn)
	setkey(cluphy.df, FASTASampleCode)
	tmp	<- list()
	for(j in seq_along(cluphy.subtrees))		#need for loop because we change 'cluphy.df' and 'clustering'
	{
		#print(j)
		#j<- 17
		x					<- cluphy.subtrees[[j]]			
		x.tiplabel			<- x$tip.label
		x.cluster			<- cluphy.df[x.tiplabel[1],cluster][,cluster]
		x.splitlabel		<- subset( cluphy.df[x.tiplabel,], eval(splitexpr), FASTASampleCode )[,FASTASampleCode]
		#print(x); print(x.splitlabel); print(cluphy.df[x$tip.label,], )
		splitlabel	<- x.splitlabel
		x.splittrees		<- hivc.phy.splittip(x, x.splitlabel)
		#from x.splittrees, select between patient clusters
		#print(x.splittrees)
		x.clusters			<- c()
		if(length(x.splittrees))
			x.clusters		<- which( sapply( seq_along(x.splittrees), function(i) length(unique(cluphy.df[x.splittrees[[i]][["tip.label"]],Patient][,Patient]))>1 	) ) 	
		x.splittrees		<- lapply( x.clusters, function(i)	x.splittrees[[i]]	)		
		#update 'clustering' and 'cluphy.df'
		if(!length(x.splittrees))			#Splitting breaks x into singletons. In 'df', set 'cluster' to NA. In 'clustering', update 'clu.idx' and 'clu.mem'.
		{
			x.clusters		<- c()
			clustering[["clu.idx"]][ x.cluster ]							<- NA
			clustering[["clu.mem"]][ clustering[["clu.mem"]]==x.cluster ]	<- NA			
			set(cluphy.df, which(cluphy.df[,FASTASampleCode%in%x.tiplabel]), "cluster", NA)
		}
		else	
		{				
			x.clusters		<- c( x.cluster, seq.int(length(clustering[["clu.idx"]])+1, len=length(x.splittrees)-1) )
			#print("HERE"); print(x.clusters); print(x.splittrees)
			for(i in seq_along(x.splittrees))
			{
				#print(i)
				#i<- 1
				#single smaller cluster. In 'df', set 'cluster' for lost x.tiplabels to NA. In 'clustering', update 'clu.idx' and 'clu.mem'.
				x.tip																<- sapply(x.splittrees[[i]][["tip.label"]], function(z) match(z, ph$tip.label) )				
				x.tip.anc															<- lapply(x.tip, function(z) Ancestors(ph, z) )
				x.tip.jnt															<- my.intersect.n( x.tip.anc )			
				x.tip.mrca															<- x.tip.anc[[1]][min(which(x.tip.anc[[1]] %in% x.tip.jnt))]				
				clustering[["clu.mem"]][ clustering[["clu.mem"]]==x.clusters[i] ]	<- NA				
				clustering[["clu.mem"]][ c(x.tip.mrca, x.tip) ]						<- x.clusters[i]												#WARNING: might miss inner nodes, so size.all would not be accurate
				clustering[["clu.idx"]][ x.clusters[i] ]							<- x.tip.mrca													#ok to assign beyond current length
				x.tiplabel															<- setdiff(x.tiplabel, x.splittrees[[i]][["tip.label"]])
				set(cluphy.df, which(cluphy.df[,FASTASampleCode%in%x.tiplabel ]), "cluster", NA)													#set lost tip labels to NA
				set(cluphy.df, which(cluphy.df[,FASTASampleCode%in%x.splittrees[[i]][["tip.label"]]]), "cluster", x.clusters[i])					#reset tip labels to x.clusters[i], as this might be a new cluster					
			}
		}
		tmp[[j]]<- list(clu= x.splittrees, clu.names=x.clusters)					
	}
	cluphy.cluidx				<- which(sapply(seq_along(tmp), function(i)		length(tmp[[i]][["clu"]])>0	))	
	cluphy.subtrees				<- eval(parse(text=paste("c(",paste('tmp[[',seq_along(tmp),']][["clu"]]', sep='',collapse=','),")",sep='')))	
	names(cluphy.subtrees)		<- eval(parse(text=paste("c(",paste('tmp[[',seq_along(tmp),']][["clu.names"]]', sep='',collapse=','),")",sep='')))	
	clustering[["size.all"]]	<- NULL
	clustering[["size.tips"]]	<- table( clustering[["clu.mem"]][1:Ntip(ph)] )
	
	list(cluphy.subtrees=cluphy.subtrees, clustering=clustering, cluphy.df=cluphy.df)
}	
######################################################################################
#modify clustering by splitting existing clusters based on 'splitexpr'. This is useful when additional meta-information is available.
hivc.clu.splitcluster<- function(cluphy.subtrees, cluphy.df, splitexpr)
{
	require(phangorn)
	setkey(cluphy.df, FASTASampleCode)
	tmp	<- list()
	for(j in seq_along(cluphy.subtrees))		#need for loop because we change 'cluphy.df' and 'clustering'
	{
		#print(j)
		#j<- 17
		x					<- cluphy.subtrees[[j]]			
		x.tiplabel			<- x$tip.label
		x.cluster			<- cluphy.df[x.tiplabel[1],cluster][,cluster]
		x.splitlabel		<- subset( cluphy.df[x.tiplabel,], eval(splitexpr), FASTASampleCode )[,FASTASampleCode]
		#print(x); print(x.splitlabel); print(cluphy.df[x$tip.label,], )
		splitlabel	<- x.splitlabel
		x.splittrees		<- hivc.phy.splittip(x, x.splitlabel)
		#from x.splittrees, select between patient clusters
		#print(x.splittrees)
		x.clusters			<- c()
		if(length(x.splittrees))
			x.clusters		<- which( sapply( seq_along(x.splittrees), function(i) length(unique(cluphy.df[x.splittrees[[i]][["tip.label"]],Patient][,Patient]))>1 	) ) 	
		x.splittrees		<- lapply( x.clusters, function(i)	x.splittrees[[i]]	)		
		
		tmp[[j]]<- list(clu= x.splittrees, clu.names=x.clusters)					
	}
	cluphy.cluidx				<- which(sapply(seq_along(tmp), function(i)		length(tmp[[i]][["clu"]])>0	))	
	cluphy.subtrees				<- eval(parse(text=paste("c(",paste('tmp[[',seq_along(tmp),']][["clu"]]', sep='',collapse=','),")",sep='')))	
	names(cluphy.subtrees)		<- eval(parse(text=paste("c(",paste('tmp[[',seq_along(tmp),']][["clu.names"]]', sep='',collapse=','),")",sep='')))	
	
	cluphy.subtrees
}	
######################################################################################
hivc.phy.splittip<- function(x, splitlabel)
{
	#str(x)
	x.splittip			<- which(x$tip.label%in%splitlabel)								#take first matching label as this 'x.splittip' and put rest on stack
	if(!length(x.splittip))		
		return(list(x))
	x.splittip			<- x.splittip[1]
	splitlabel.stack	<- setdiff(splitlabel, x$tip.label[ x.splittip ])	
	#print(x$tip.label); print(x.splittip); print(splitlabel.stack)	
	if(Ntip(x)<3)																		#tip is among tiplabels, so splitting will break pair
		return(list())
	x.splitnode			<- x$edge[x$edge[,2]==x.splittip, 1]							#ancestor of HETF tip		
	#x$tip.label[x.splittip]
	x.splittree2		<- extract.clade(x, x.splitnode, root.edge = 1)
	
	#replace 'x.splitnode' with 'x.splittip' and drop dummy tip	
	x.rmnode							<- x.splitnode
	tmp									<- x$edge[x$edge[,1]==x.splitnode,2] 			
	while(any(tmp>Ntip(x)))
	{
		x.rmnode	<- c(x.rmnode, tmp[ tmp>Ntip(x) ])
		tmp			<- x$edge[x$edge[,1]%in%tmp,2]
	}
	x$node.label						<- x$node.label[ -(x.rmnode-Ntip(x)) ]
	tmp									<- x$edge[,1]%in%x.rmnode
	x$edge.length						<- x$edge.length[ !tmp ]	
	x.rmtip								<- x$edge[tmp,2]
	x.rmtip								<- x.rmtip[x.rmtip!=x.splittip & x.rmtip<=Ntip(x)]	 
	x$tip.label							<- x$tip.label[-x.rmtip]
	x$Nnode								<- x$Nnode-length(x.rmnode)	
	x$edge								<- x$edge[!x$edge[,1]%in%x.rmnode,]				#drop clade
	x.splitedge							<- x$edge[,2]==x.splitnode						#find edge leading to clade	
	x$edge[x.splitedge,2]				<- x.splittip									#set edge to dummy tip
	if(nrow(x$edge))
	{
		tmp								<- as.vector(x$edge)							#need to get node numbers right again
		tmp.rm							<- which(!seq_len(max(tmp))%in%tmp)
		for(z in rev(tmp.rm))
			tmp[tmp>=z]					<- tmp[tmp>=z]-1			
		x$edge							<- matrix(tmp, nrow=nrow(x$edge), ncol=ncol(x$edge) )
		x.splittip						<- x$edge[x.splitedge,2]						#renumbering might have changed dummy tip
	}
	else
	{
		x.splittip						<- 1
	}
	ans							<- list()
	#print("HERE"); str(x); print(x.splittip); str(x.splittree2); print(x$tip.label[x.splittip]); print("HERE2")
	if(Ntip(x)>2)
		ans[[1]]				<- drop.tip(x, x.splittip, root.edge=1)
	if(Ntip(x.splittree2)>2)
		ans[[length(ans)+1]]	<- drop.tip(x.splittree2, x$tip.label[x.splittip], root.edge=1)
	if(length(splitlabel.stack) && length(ans))
	{
		tmp	<- paste("c(",paste('hivc.phy.splittip(ans[[',seq_along(ans),']], splitlabel.stack)', sep='',collapse=','),")",sep='')		
		tmp	<- eval(parse(text=tmp))																	#annoyingly, sapply etc is not the same as c( etc )		
		#tmp	<- sapply(seq_along(ans), function(i)	hivc.phy.splittip(ans[[i]], splitlabel.stack) )		#c(list()) is [[1]] list(), need a hack
		#ans	<- list()
		#ans	<- lapply( which(sapply(tmp, length)>0), function(i) tmp[[i]] )
		ans <- tmp
	}
	ans
}


