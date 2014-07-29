######################################################################################
r3seq.guess.related.sequences<- function(triplet, seq, seq.select.n, seq.select.f=5, verbose=10)
{
	#	create sequence matrices corresponding to the two breakpoint regions 
	seq.in		<- seq[,seq.int(triplet[,bp1.1],triplet[,bp1.2])]
	seq.out		<- if(triplet[,child.start]<triplet[,bp1.1]-1) seq.int(triplet[,child.start],triplet[,bp1.1]-1) else numeric(0) 
	seq.out		<- if(triplet[,bp1.2]+1<triplet[,child.len]) c(seq.out,seq.int(triplet[,bp1.2]+1, triplet[,child.len]))	else 	seq.out
	seq.out		<- seq[,seq.out]
	seq.select.f<- ifelse(min(ncol(seq.out),ncol(seq.in))<150, seq.select.f, seq.select.f)
	if(verbose)	cat(paste("\nsetting inflation factor to",seq.select.f))
	if(is.na(seq.select.n))
		seq.select.n<- 10 * seq.select.f		
	#	select background sequences for child based on sequence similarity
	if(verbose)	cat(paste("\ncompute genetic distances for parent1 parent2 child"))
	tmp				<- which( rownames(seq)==triplet[,child] )			
	seq.dist					<- sapply(seq_len(nrow(seq.in))[-tmp],function(i){		seq.dist.pairwise( seq.in[tmp,], seq.in[i,])	})	
	seq.dist[is.nan(seq.dist)]	<- Inf
	seq.df			<- data.table( FASTASampleCode=rownames(seq.in)[-tmp] , dist=seq.dist, group="child", region="in" ) 		
	seq.dist					<- 1 - sapply(seq_len(nrow(seq.out))[-tmp],function(i){		seq.dist.pairwise(seq.out[tmp,], seq.out[i,])	})	
	seq.dist[is.nan(seq.dist)]	<- Inf
	seq.df			<- rbind(seq.df, data.table( FASTASampleCode=rownames(seq.out)[-tmp] , dist=seq.dist, group="child", region="out" )) 
	#	select background sequences for parent1 based on sequence similarity
	tmp				<- which( rownames(seq)==triplet[,parent1] )				
	seq.dist					<- 1 - sapply(seq_len(nrow(seq.in))[-tmp],function(i){		seq.dist.pairwise(seq.in[tmp,], seq.in[i,])		})	
	seq.dist[is.nan(seq.dist)]	<- Inf
	seq.df			<- rbind(seq.df,data.table( FASTASampleCode=rownames(seq.in)[-tmp] , dist=seq.dist, group="parent1", region="in" )) 		
	seq.dist					<- 1 - sapply(seq_len(nrow(seq.out))[-tmp],function(i){		seq.dist.pairwise( seq.out[tmp,], seq.out[i,])		})	
	seq.dist[is.nan(seq.dist)]	<- Inf
	seq.df			<- rbind(seq.df, data.table( FASTASampleCode=rownames(seq.out)[-tmp] , dist=seq.dist, group="parent1", region="out" )) 
	#	select background sequences for parent2 based on sequence similarity
	tmp				<- which( rownames(seq)==triplet[,parent2] )				
	seq.dist					<- 1 - sapply(seq_len(nrow(seq.in))[-tmp],function(i){		seq.dist.pairwise(seq.in[tmp,], seq.in[i,])			})	
	seq.dist[is.nan(seq.dist)]	<- Inf
	seq.df			<- rbind(seq.df,data.table( FASTASampleCode=rownames(seq.in)[-tmp] , dist=seq.dist, group="parent2", region="in" )) 		
	seq.dist					<- 1 - sapply(seq_len(nrow(seq.out))[-tmp],function(i){		seq.dist.pairwise(seq.out[tmp,], seq.out[i,])	})	
	seq.dist[is.nan(seq.dist)]	<- Inf
	seq.df			<- rbind(seq.df, data.table( FASTASampleCode=rownames(seq.out)[-tmp] , dist=seq.dist, group="parent2", region="out" ))
	#	first pass:
	#	select closest FASTASampleCode by group and region to verify recombination breakpoint by phylogenetic incongruence
	setkeyv(seq.df, c("group","region","dist"))
	seq.df			<- subset(seq.df, dist>0)
	if(verbose)	cat(paste("\nFound related sequences with dist>0, n=",length(seq.df[,unique(FASTASampleCode)])))
	#	for each group and region, select the n closest sequence names (non-unique)
	seq.df			<- seq.df[	,	list(FASTASampleCode=FASTASampleCode[seq_len(seq.select.n)], dist=dist[seq_len(seq.select.n)]), by=c("group","region")]
	#	keep each sequence name once, for the group it is closest to		
	seq.df			<- seq.df[, {
				tmp<- which.min(dist)
				list(dist= dist[tmp], group=group[tmp], region=region[tmp])
			}, by=c("FASTASampleCode")]
	setkeyv(seq.df, c("group","region","dist"))					
	if(verbose)	cat(paste("\ndetermined candidates for balancing filler sequences, n=",nrow(seq.df)))
	#	select unique sequences
	tmp				<- c( triplet[,parent1],triplet[,parent2],triplet[,child], seq.df[, FASTASampleCode] )
	tmp				<- seq.unique(seq.in[tmp,])
	seq.df			<- merge( data.table(FASTASampleCode=rownames(tmp)), seq.df, by="FASTASampleCode" )
	tmp				<- c( triplet[,parent1],triplet[,parent2],triplet[,child], seq.df[, FASTASampleCode] )
	tmp				<- seq.unique(seq.out[tmp,])
	seq.df			<- merge( data.table(FASTASampleCode=rownames(tmp)), seq.df, by="FASTASampleCode" )
	seq.df			<- subset(seq.df, !FASTASampleCode%in%c( triplet[,parent1],triplet[,parent2],triplet[,child] ) )		
	if(verbose)	cat(paste("\nfound candidates for filler sequences that are unique on both recombinant regions, n=",nrow(seq.df)))
	if(verbose)	print( seq.df[	,	list(n=length(FASTASampleCode)) ,by=c("group","region")] )
	#		
	seq.select.n	<- seq.select.n/seq.select.f 
	#	select 'seq.select.n' unique filler sequences for 'in' region, balancing by group as much as possible
	if(verbose)	cat(paste("\nSelect filler sequences for recombinant region 'in'"))
	tmp						<- subset( seq.df, region=='in' )[,FASTASampleCode]
	tmp						<- seq.unique(seq.in[tmp ,])
	seq.in.df				<- merge( data.table(FASTASampleCode=rownames(tmp)), subset(seq.df,region=="in"), by="FASTASampleCode" )
	setkey(seq.in.df, dist)
	if(nrow(seq.in.df)<seq.select.n)	cat(paste("\ncan only select less than the requested number of sequences, n=",nrow(seq.in.df)))
	tmp						<- rbind( seq.in.df, data.table(FASTASampleCode=NA, dist=NA, group=c("child","parent1","parent2"), region=NA) )
	seq.in.order			<- tmp[	,	list(n=length(na.omit(FASTASampleCode))) ,by=c("group")]		
	seq.in.order			<- seq.in.order[order(n),]
	overflow				<- 0
	ans						<- data.table(FASTASampleCode=NA, dist=NA, group=NA, region=NA)
	for(x in seq.in.order[,group])
	{			
		#print(x)
		tmp					<- subset(seq.in.df, group==x)
		#print(tmp)
		ans					<- rbind(tmp[seq_len( min(seq.select.n+overflow, nrow(tmp)) ),], ans )
		overflow			<- ifelse(seq.select.n+overflow<nrow(tmp), 0, seq.select.n+overflow-nrow(tmp))
	}
	if(overflow>0)	cat(paste("\nNot as many filler sequences as requested for recombinant region 'in', n=",nrow(seq.in.df)))
	seq.in.df				<- ans[-nrow(ans),]
	if(verbose)	cat(paste("\nSelected balancing sequences for recombinant region 'in', n=",nrow(seq.in.df)))
	if(verbose)	print( seq.in.df[	,	list(n=length(FASTASampleCode)) ,by=c("group","region")] )
	#
	#	select 'seq.select.n' unique filler sequences for 'out' region, balancing by group as much as possible
	#
	if(verbose)	cat(paste("\nSelect sequences for recombinant region 'out'"))
	tmp						<- subset( seq.df, region=='out' )[,FASTASampleCode] 
	tmp						<- seq.unique(seq.out[tmp,])		
	seq.out.df				<- merge( data.table(FASTASampleCode=rownames(tmp)), subset(seq.df,region=="out"), by="FASTASampleCode" )
	setkey(seq.out.df, dist)
	if(nrow(seq.out.df)<seq.select.n)	cat(paste("\ncan only select less than the requested number of sequences, n=",nrow(seq.out.df)))
	tmp						<- rbind( seq.out.df, data.table(FASTASampleCode=NA, dist=NA, group=c("child","parent1","parent2"), region=NA) )
	seq.out.order			<- tmp[	,	list(n=length(FASTASampleCode)) ,by=c("group")]		
	seq.out.order			<- seq.out.order[order(n),]
	overflow				<- 0
	ans						<- data.table(FASTASampleCode=NA, dist=NA, group=NA, region=NA)		
	for(x in seq.out.order[,group])
	{			
		#print(x)
		tmp					<- subset(seq.out.df, group==x)
		#print(tmp)
		ans					<- rbind(tmp[seq_len( min(seq.select.n+overflow, nrow(tmp)) ),], ans )
		overflow			<- ifelse(seq.select.n+overflow<nrow(tmp), 0, seq.select.n+overflow-nrow(tmp))
	}
	if(overflow>0)	cat(paste("\nNot as many filler sequences as requested for recombinant region 'out', n=",nrow(seq.out.df)))
	seq.out.df				<- ans[-nrow(ans),]
	if(verbose)	cat(paste("\nSelected balancing sequences for recombinant region 'out', n=",nrow(seq.out.df)))
	if(verbose)	print( seq.out.df[	,	list(n=length(FASTASampleCode)) ,by=c("group","region")] )
	#	combine unique filler sequences
	seq.df					<- rbind( seq.in.df, seq.out.df )
	if(verbose)	cat(paste("\nSelected balancing set of closest filler sequences"))
	if(verbose)	print( seq.df[	,	list(n=length(FASTASampleCode)) ,by=c("group","region")] )		
	if( length(unique( seq.df[, FASTASampleCode] ))!=nrow(seq.df) )		stop("Unexpected non-unique sequence names")
	if( nrow(seq.unique( seq.in[ seq.df[, FASTASampleCode], ] ))!=nrow(seq.df) ) stop("Unexpected non-unique 'in' sequences")
	if( nrow(seq.unique( seq.out[ seq.df[, FASTASampleCode], ] ))!=nrow(seq.df) ) stop("Unexpected non-unique 'out' sequences")
	#	could be that the triplet sequences are not unique among each other 
	#	if so, make change to one of the triplet sequences
	seq.select				<- c( triplet[,parent1],triplet[,parent2],triplet[,child], seq.df[, FASTASampleCode] )
	tmp						<- seq.unique( seq.in[ seq.select, ] )
	if( !length(setdiff(c( triplet[,parent1],triplet[,parent2],triplet[,child] ),rownames(tmp))) )
		seq.in				<- tmp
	if( length(setdiff(c( triplet[,parent1],triplet[,parent2],triplet[,child] ),rownames(tmp))) )
	{
		tmp					<- setdiff(c( triplet[,parent1],triplet[,parent2],triplet[,child] ),rownames(tmp) )		#name of sequence in triplet that is identical with one other sequence in triplet
		if(verbose)	cat(paste("\nFound identical triplet sequence for region 'in'", tmp))
		seq.in				<- as.character(seq.in)
		seq.in[tmp,1]		<- ifelse(seq.in[tmp,1]=='t','c',ifelse(seq.in[tmp,1]=='c','t',ifelse(seq.in[tmp,1]=='a','g','a')))
		seq.in				<- as.DNAbin(seq.in)
		tmp					<- seq.unique( seq.in[ seq.select, ] )
		if( length(setdiff(c( triplet[,parent1],triplet[,parent2],triplet[,child] ),rownames(tmp))) )	stop("Unexpected duplicate for 'in'")
		seq.in				<- tmp
	}			
	seq.out					<- seq.out[seq.select,]
	tmp						<- seq.unique(seq.out)
	if( !length(setdiff(c( triplet[,parent1],triplet[,parent2],triplet[,child] ),rownames(tmp))) )
		seq.out				<- tmp
	if( length(setdiff(c( triplet[,parent1],triplet[,parent2],triplet[,child] ),rownames(tmp))) )
	{
		tmp					<- setdiff(c( triplet[,parent1],triplet[,parent2],triplet[,child] ),rownames(tmp) )		#name of sequence in triplet that is identical with one other sequence in triplet
		if(verbose)	cat(paste("\nFound identical triplet sequence for region 'out'", tmp))
		seq.out				<- as.character(seq.out)
		seq.out[tmp,1]		<- ifelse(seq.out[tmp,1]=='t','c',ifelse(seq.out[tmp,1]=='c','t',ifelse(seq.out[tmp,1]=='a','g','a')))
		seq.out				<- as.DNAbin(seq.out)
		tmp					<- seq.unique(seq.out[seq.select,])
		if( length(setdiff(c( triplet[,parent1],triplet[,parent2],triplet[,child] ),rownames(tmp))) )	stop("Unexpected duplicate for 'out'")
		seq.out				<- tmp
	}
	if(any(rownames(seq.in)!=rownames(seq.out)))	stop("Unexpected unequal sequences selected")
	#	reset rownames
	tmp						<- c( paste(c("tparent1","tparent2","tchild"),seq.select[1:3],sep='_'), seq.df[,list(label= paste(region,group,FASTASampleCode,sep='_')), by="FASTASampleCode"][,label] )
	rownames(seq.in)		<- tmp
	rownames(seq.out)		<- tmp
	
	list(seq.in=seq.in, seq.out=seq.out)
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



