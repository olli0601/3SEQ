#' this file contains all R functions of the hivclust package
#' @useDynLib hivc
#' @import ape
#' @import phytools
#' @import igraph
#' @import geiger
#' @import data.table

#' @export
HIVC.db.locktime	<-  as.Date("30/03/2013", format="%d/%m/%Y")


######################################################################################
#' @export
hivc.seq.read.GenBank<- function (access.nb, seq.names = access.nb, species.names = TRUE, gene.names = FALSE, as.character = FALSE, attributes= c("origin")) 
{
	require(ape)
	N <- length(access.nb)
	nrequest <- N%/%400 + as.logical(N%%400)
	X <- character(0)
	for (i in 1:nrequest) {
			a <- (i - 1) * 400 + 1
			b <- 400 * i
			if (i == nrequest) 
				b <- N
			URL <- paste("http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=", 
					paste(access.nb[a:b], collapse = ","), "&rettype=gb&retmode=text", 
					sep = "")
			X <- c(X, scan(file = URL, what = "", sep = "\n", quiet = TRUE))
		}
		
	FI <- grep("^ {0,}ORIGIN", X) + 1
	LA <- which(X == "//") - 1
	obj <- vector("list", N)
	for (i in 1:N) {
		tmp <- gsub("[[:digit:] ]", "", X[FI[i]:LA[i]])
		obj[[i]] <- unlist(strsplit(tmp, NULL))
	}
	names(obj) <- seq.names
	if (!as.character) 
		obj <- as.DNAbin(obj)
	if(length(attributes)) 
	{
		attr.lines	<- lapply(attributes, function(attr) grep(paste("^ {0,}/",attr,sep=''), X) 	)
		attr.fields	<- lapply(seq_along(attr.lines),function(i)		gsub("\"$","",gsub(paste("^ {0,}/",attributes[i],"=\"",sep=''),"",X[attr.lines[[i]]]))		)						
		for(i in seq_along(attributes))
			attr(obj, attributes[[i]])<- attr.fields[[i]]
	}
	if (gene.names) {
		tmp <- character(N)
		sp <- grep(" +gene +<", X)
		for (i in 1:N) tmp[i] <- unlist(strsplit(X[sp[i + 1L]], 
							" +/gene=\""))[2]
		attr(obj, "gene") <- gsub("\"$", "", tmp)
	}
	obj
}
######################################################################################
hivc.db.resetNegTbyPoslRNA_T1<- function(df.all,verbose=1)
{
	if(verbose)	cat(paste("\nreset NegT to NA for the following seq"))
	if(verbose) print(subset( df.all,!is.na(PoslRNA_T1) & PoslRNA_T1<NegT ))
	tmp		<- which(df.all[, NegT_Acc=="Yes" & !is.na(PoslRNA_T1) & PoslRNA_T1<NegT])
	if(verbose) cat(paste("\nnumber of seq with NegT_Acc=='Yes' & !is.na(PoslRNA_T1) & NegT > time of VL test, n=",length(tmp)))			
	set(df.all, tmp, "NegT", NA)
	tmp2	<- as.character(df.all[,NegT_Acc])
	tmp2[tmp]<- NA
	set(df.all, NULL, "NegT_Acc", as.factor(tmp2))		
	tmp		<- which(df.all[, NegT_Acc=="No" & !is.na(PoslRNA_T1) & PoslRNA_T1<NegT])
	if(verbose) cat(paste("\nnumber of seq with NegT_Acc=='No' & !is.na(PoslRNA_T1) & NegT > time of VL test, n=",length(tmp)))			
	set(df.all, tmp, "NegT", NA)
	tmp2	<- as.character(df.all[,NegT_Acc])
	tmp2[tmp]<- NA
	set(df.all, NULL, "NegT_Acc", as.factor(tmp2))
	df.all
}
######################################################################################
hivc.db.resetNegTbyAnyPosT<- function(df.all,verbose=1)
{	
	tmp		<- which(df.all[, NegT_Acc=="Yes" & NegT>AnyPos_T1])
	if(verbose) cat(paste("\nnumber of seq with NegT_Acc=='Yes' & NegT> pos test or pos seq, n=",length(tmp)))
	if(verbose) print(as.character(unique(df.all[tmp,Patient])))
	set(df.all, tmp, "NegT", NA)
	tmp2	<- as.character(df.all[,NegT_Acc])
	tmp2[tmp]<- NA
	set(df.all, NULL, "NegT_Acc", as.factor(tmp2))		
	tmp		<- which(df.all[, NegT_Acc=="No" & NegT>AnyPos_T1])
	if(verbose) cat(paste("\nnumber of seq with NegT_Acc=='No' & NegT> pos test or pos seq, n=",length(tmp)))
	if(verbose) print(as.character(unique(df.all[tmp,Patient])))
	set(df.all, tmp, "NegT", NA)
	tmp2	<- as.character(df.all[,NegT_Acc])
	tmp2[tmp]<- NA
	set(df.all, NULL, "NegT_Acc", as.factor(tmp2))
	df.all
}		
######################################################################################
hivc.db.resetCD4byNegT<- function(df.cross, with.NegT_Acc.No=0, verbose=1)
{
	tmp			<- which(df.cross[,NegT_Acc=="Yes" & 	NegT>PosCD4])
	if(verbose)		cat(paste("\nnumber of entries with NegT_Acc=='Yes' & 	NegT>PosCD4, n=", length(tmp),"SETTING CD4 to NA"))
	set(df.cross, tmp, "PosCD4", NA)
	set(df.cross, tmp, "CD4", NA)
	df.cross	<- subset(df.cross,!is.na(CD4))
	if(!with.NegT_Acc.No)
		return(df.cross)
	tmp			<- which(df.cross[,NegT_Acc=="No" & 	NegT>PosCD4])
	if(verbose)		cat(paste("\nnumber of entries with NegT_Acc=='No' & 	NegT>PosCD4, n=", length(tmp),"SETTING CD4 to NA"))
	set(df.cross, tmp, "PosCD4", NA)
	set(df.cross, tmp, "CD4", NA)
	df.cross	<- subset(df.cross,!is.na(CD4))
	df.cross
}	
######################################################################################
hivc.db.getTrIMo<- function(df.cross, verbose=1)
{
	ans		<- df.cross[, 	{ 
								x		<- data.table( PosSeqT, StartTime, StopTime, TrI )
								x		<- subset(x, !is.na(StartTime))		#discard first time period if unknown StartTime					
								if(!nrow(x) || is.na(x[1,PosSeqT]))
								{
									TrImo_bTS		<- TrImo_aTS	<- NA_real_
								}
								else if(x[1,PosSeqT]<=x[1,StartTime])
								{
									TrImo_bTS		<- 0
									TrImo_aTS		<- sum(as.numeric(subset(x, TrI=="Yes")[,difftime(StopTime,StartTime,units="days")/30]))
								}
								else
								{
									z				<- which( as.numeric( x[,difftime(PosSeqT,StartTime,units="days")] )>0 )	#all rows in which PosSeqT>=StartTime						
									z				<- z[length(z)]
									z2				<- seq.int(1,z)
									xb				<- x[z2,]
									set(xb,z,"StopTime",x[z,PosSeqT])
									z2				<- seq.int(z,nrow(x))
									xa				<- x[z2,]
									set(xa,1L,"StartTime",x[z,PosSeqT])		
									z				<- as.numeric(xb[,difftime(StopTime,StartTime,units="days")/30])
									z2				<- as.numeric(xa[,difftime(StopTime,StartTime,units="days")/30])
									TrImo_bTS		<- sum(z[ which(xb[, TrI=="Yes"]) ])
									TrImo_aTS		<- sum(z2[ which(xa[, TrI=="Yes"]) ])
								}					
								list( 	TrImo_bTS	= TrImo_bTS, 
										TrImo_aTS	= TrImo_aTS			)																					
							}, by=FASTASampleCode]
	set(ans,NULL,"TrImo_bTS",round(ans[,TrImo_bTS],d=1))
	set(ans,NULL,"TrImo_aTS",round(ans[,TrImo_aTS],d=1))
	if(verbose)	cat(paste("\nnumber of entries, n=",nrow(ans)))
	ans		
}
######################################################################################
hivc.db.getlRNA.T1andTS<- function(df.cross, lRNA.bTS.quantile= 0.75, lRNA.aTS.quantile= 0.25, lRNAi.min= log10(1e4), db.locktime= HIVC.db.locktime, verbose=1)
{						
	if(verbose)	cat(paste("\nlRNA.aTS.quantile is",lRNA.aTS.quantile))
	if(verbose)	cat(paste("\nlRNA.bTS.quantile is",lRNA.bTS.quantile))
	if(verbose)	cat(paste("\nlRNAi.min is",lRNAi.min))
	if(verbose)	cat(paste("\nnumber of entries in cross product is, n=",nrow(df.cross)))
	df.cross[, lRNAi:=lRNA>=lRNAi.min]			
	ans		<-	df.cross[, 	{
								z	<- which.min(PosRNA)
								if(!is.na(PosSeqT[1]))
								{
									z2			<- which.min(abs(difftime(PosSeqT, PosRNA, units="weeks")))
									PoslRNA_TS	<- PosRNA[z2]
									lRNA_TS		<- round( mean( lRNA[seq.int(z2-1,z2+1)], na.rm=1 ), d=1 )				
									lRNA_bTS	<- round( quantile(lRNA[seq_len(z2+1)], probs = lRNA.bTS.quantile, na.rm=T, names = F), d=1 )
									lRNA_aTS	<- round( quantile(lRNA[seq.int(z2-1,length(lRNA))], probs = lRNA.aTS.quantile, na.rm=T, names = F), d=1 )							
									lRNAit		<- as.numeric( difftime(c(PosRNA[-1],db.locktime), PosRNA,units="days") )
									#print(data.table(Patient, FASTASampleCode, AnyPos_T1, PosSeqT, PosRNA, lRNA, lRNAi, lRNA.hb4tr_LT, lRNA.early)); print(lRNAit); print(z2); print(lRNAit[ seq.int(z2,length(lRNA)) ]); print(lRNAit[ seq_len(z2) ][  lRNAi[ seq_len(z2) ]  ])
									lRNAi_bTS	<- sum( lRNAit[ seq_len(z2) ][  lRNAi[ seq_len(z2) ]  ] ) / sum( lRNAit[ seq_len(z2) ] ) 							
									lRNAi_aTS	<- sum( lRNAit[seq.int(z2,length(lRNAit))][  lRNAi[ seq.int(z2,length(lRNAi)) ]	] ) /	sum( lRNAit[seq.int(z2,length(lRNAit))] )
								}
								else
								{
									PoslRNA_TS	<- lRNAi_bTS<- lRNAi_aTS<- as.Date(NA)
									lRNA_TS 	<- lRNA_bTS <- lRNA_aTS <- NA_real_
								}
								list(	PoslRNA_T1		= PosRNA[z], 
										lRNA_T1			= lRNA[z], 
										PoslRNA_TS		= PoslRNA_TS, 
										lRNA_TS			= lRNA_TS, 
										lRNA_bTS		= lRNA_bTS, 
										lRNA_aTS		= lRNA_aTS,
										lRNAi_bTS		= lRNAi_bTS,
										lRNAi_aTS		= lRNAi_aTS,
										lRNA.hb4tr_LT	= lRNA.hb4tr_LT[1], 
										lRNA.early		= lRNA.early[1]				) 							 	
							},by=FASTASampleCode]
	if(verbose)	cat(paste("\nnumber of seq with PosCD4_T1 CD4_T1  PosCD4_TS CD4_TS is n=",nrow(ans)))
	df.cross[,lRNAi:=NULL]
	ans	
}
######################################################################################
hivc.db.getcoverage<- function(df)
{
	ans	<- paste( c("df[,list(", paste( 	sapply(colnames(df), function(x) paste(x,"= length(which(!is.na(",x,")))",sep='')), collapse=",", sep='' ),")]") , collapse='',sep='')
	ans	<- eval(parse(text=ans))	/ nrow(df)
	ans
}
######################################################################################
hivc.db.getCD4.T1andTS<- function(df.cross, verbose=1, CD4.HIVNeg.min= 500, CD4.bTS.quantile= 0.75, CD4.aTS.quantile= 0.25)
{
	if(verbose)	cat(paste("\nCD4.HIVNeg.min is",CD4.HIVNeg.min))
	if(verbose)	cat(paste("\nCD4.aTS.quantile is",CD4.aTS.quantile))
	if(verbose)	cat(paste("\nCD4.bTS.quantile is",CD4.bTS.quantile))
	if(verbose)	cat(paste("\nnumber of entries in cross product is, n=",nrow(df.cross)))
	tmp	<- which( df.cross[, AnyPos_T1>PosCD4 & CD4>CD4.HIVNeg.min] )
	if(verbose)	cat(paste("\nnumber of patients with AnyPos_T1>PosCD4 & CD4>CD4.HIVNeg.min, n=",length(unique(df.cross[tmp,Patient])),"SETTING PosCD4 to NA since Patients could be healthy for corresponding CD4"))
	set(df.cross, tmp, "CD4", NA)		
	df.cross		<- subset(df.cross, !is.na(CD4))
	if(verbose)	cat(paste("\nnumber of entries in cross product is, n=",nrow(df.cross)))
	
	ans		<-	df.cross[, 	{
				z	<- which.min(PosCD4)
				if(!is.na(PosSeqT[1]))
				{
					z2			<- which.min(abs(difftime(PosSeqT, PosCD4, units="weeks")))
					PosCD4_TS	<- PosCD4[z2]
					CD4_TS		<- round( mean( CD4[seq.int(z2-1,z2+1)], na.rm=1 ), d=0 )				
					CD4_bTS		<- round( quantile(CD4[seq_len(z2+1)], probs = CD4.bTS.quantile, na.rm=T, names = F), d=0 )
					CD4_aTS		<- round( quantile(CD4[seq.int(z2-1,length(CD4))], probs = CD4.aTS.quantile, na.rm=T, names = F), d=0 )
				}
				else
				{
					PosCD4_TS	<- as.Date(NA)
					CD4_TS <- CD4_bTS <- CD4_aTS <- NA_real_
				}
				list(PosCD4_T1=PosCD4[z], CD4_T1=CD4[z], PosCD4_TS=PosCD4_TS, CD4_TS=CD4_TS, CD4_bTS=CD4_bTS, CD4_aTS=CD4_aTS	 ) 	
			},by=FASTASampleCode]
	
	if(verbose)	cat(paste("\nnumber of seq with PosCD4_T1 CD4_T1  PosCD4_TS CD4_TS is n=",nrow(ans)))
	ans
}
######################################################################################
hivc.db.reset.inaccuratePosT<- function(df, nacc.dy.dy= 30, nacc.mody.mo= 11, nacc.mody.dy= 31, verbose=1)
{	
	nacc.dy		<- which( df[, PosT_Acc=="No" & !is.na(PosT) & as.POSIXlt(PosT)$mday==15] )
	nacc.mody	<- which( df[, PosT_Acc=="No" & !is.na(PosT) & as.POSIXlt(PosT)$mon==6 & as.POSIXlt(PosT)$mday==1] )
	if(verbose) cat(paste("\nnumber of uncertain PosT, day only, n=",length(nacc.dy)))
	if(verbose) cat(paste("\nnumber of uncertain PosT, month & day, n=",length(nacc.mody)))
	# reset nacc.dy
	tmp							<- as.POSIXlt(df[nacc.dy,PosT] )
	tmp$mday					<- nacc.dy.dy
	set(df, nacc.dy, "PosT", as.Date(tmp))
	# reset nacc.mody
	tmp							<- as.POSIXlt(df[nacc.mody,PosT] )
	tmp$mday					<- nacc.mody.dy
	tmp$mon						<- nacc.mody.mo
	set(df, nacc.mody, "PosT", as.Date(tmp))
	df
}
######################################################################################
hivc.db.reset.inaccurateNegT<- function(df, nacc.dy.dy= 1, nacc.mody.mo= 0, nacc.mody.dy= 1, verbose=1)
{	
	serocon.nacc.dy		<- which( df[, NegT_Acc=="No" & !is.na(NegT) & as.POSIXlt(NegT)$mday==15] )
	serocon.nacc.mody	<- which( df[, NegT_Acc=="No" & !is.na(NegT) & as.POSIXlt(NegT)$mon==6 & as.POSIXlt(NegT)$mday==1] )
	if(verbose) cat(paste("\nnumber of uncertain NegT, day only, n=",length(serocon.nacc.dy)))
	if(verbose) cat(paste("\nnumber of uncertain NegT, month & day, n=",length(serocon.nacc.mody)))
	# reset serocon.nacc.dy
	tmp							<- as.POSIXlt(df[serocon.nacc.dy,NegT] )
	tmp$mday					<- nacc.dy.dy
	set(df, serocon.nacc.dy, "NegT", as.Date(tmp))
	#set(df, serocon.nacc.dy, "NegT_Acc", "Yes")
	# reset serocon.nacc.mody
	tmp							<- as.POSIXlt(df[serocon.nacc.mody,NegT] )
	tmp$mday					<- nacc.mody.dy
	tmp$mon						<- nacc.mody.mo
	set(df, serocon.nacc.mody, "NegT", as.Date(tmp))
	#set(df, serocon.nacc.mody, "NegT_Acc", "Yes")
	df
}
######################################################################################
hivc.db.getplot.newdiagnosesbyCD4<- function(df, plot.file=NULL, plot.file.p=NULL, plot.ylab=NULL, verbose=1)
{
	set(df, NULL, "AnyPos_T1", hivc.db.Date2numeric(df[,AnyPos_T1]))
	df[,AnyPos_yr:=	floor(df[,AnyPos_T1])]
	df[,CD4_T1bin:= my.aggregate(df[,CD4_T1], c(-Inf, 200, 350, 500, Inf))]
	#
	t.newdiag			<- table( df[,AnyPos_yr,CD4_T1bin] )
	rownames(t.newdiag)	<- c("<200","200-349","350-499",">=500")	
	if(!is.null(plot.file))
	{
		cols				<- brewer.pal(nrow(t.newdiag),"Set1")
		if(verbose)	cat(paste("\nplot file prop to ",plot.file))
		pdf(file=plot.file, width=5, height=4)
		par(mar=c(4,6,0.5,0.5))
		my.barplot.table(t.newdiag, 0.4, "year", plot.ylab, cols, x.baseline=0, xax=as.numeric(colnames(t.newdiag)), legend.loc="topleft")
		dev.off()
	}
	#
	p.newdiag			<- t.newdiag / matrix(rep(apply(t.newdiag,2,sum), each=nrow(t.newdiag)), nrow(t.newdiag), ncol(t.newdiag))
	if(!is.null(plot.file.p))
	{
		if(verbose)	cat(paste("\nplot file prop to ",plot.file.p))
		pdf(file=plot.file.p, width=5, height=4)
		par(mar=c(4,6,0.5,0.5))
		cols				<- brewer.pal(nrow(p.newdiag),"Set1")
		my.barplot.table(p.newdiag, 0.4, "year", paste("%",plot.ylab), cols, x.baseline=0, xax=as.numeric(colnames(p.newdiag)), legend.loc=NULL)
		dev.off()
	}
	list( t.newdiag=t.newdiag, p.newdiag=p.newdiag)
}
######################################################################################
hivc.db.getplot.livingbyCD4<- function(df, df.immu, plot.file, plot.file.p, plot.ylab, db.endtime=2013.3, db.diff.lastcontact2died=0.5, db.diff.lastcontact2now= 2.3, verbose=1)
{
	#
	set(df, NULL, "AnyPos_T1", hivc.db.Date2numeric(df[,AnyPos_T1]))
	set(df, NULL, "DateDied", hivc.db.Date2numeric(df[,DateDied]))
	set(df, NULL, "DateLastContact", hivc.db.Date2numeric(df[,DateLastContact]))
	set(df, NULL, "AnyT_T1", hivc.db.Date2numeric(df[,AnyT_T1]))
	df[,DateEnd:= DateDied]
	#compute DateEnd, either death or lost contact
	tmp				<- which( df[, (DateDied-DateLastContact)>db.diff.lastcontact2died] )
	if(verbose)	cat(paste("\n(DateDied-DateLastContact)>db.diff.lastcontact2died for n=",length(tmp)))
	set(df, tmp, "DateEnd", df[tmp,DateLastContact])
	tmp				<- which(df[, is.na(DateDied) & (db.endtime-DateLastContact)>db.diff.lastcontact2now])
	if(verbose)	cat(paste("\n(db.endtime-DateLastContact)>db.diff.lastcontact2now for n=",length(tmp)))
	set(df, tmp, "DateEnd", df[tmp,DateLastContact])
	tmp				<- which(is.na(df[,DateLastContact]))
	if(verbose)	cat(paste("\nnumber patients with missing DateLastContact n=",length(tmp)))
	tmp				<- which(is.na(df[,DateEnd]))
	if(verbose)	cat(paste("\nnumber patients still alive n=",length(tmp)))
	set(df, tmp,"DateEnd",db.endtime)
	df[, AnyPos_yr:= floor(AnyPos_T1)]
	df[, DateEnd_yr:= floor(DateEnd)]
	#
	set(df, which(is.na(df[,AnyT_T1])), "AnyT_T1", max(df[,DateEnd_yr])+1)
	#
	df				<- subset(df, select=c(Patient, AnyPos_T1, CD4_T1, isAcute, AnyT_T1, DateEnd, AnyPos_yr, DateEnd_yr))	
	#compute patients alive in year Db_yr 
	tmp				<- seq.int(min(df[, AnyPos_yr]),max(df[, DateEnd_yr]))
	tmp				<- lapply(tmp, function(x)
			{				
				tmp	<- subset(df, AnyPos_yr<=x & DateEnd_yr>=x)
				cbind(tmp,data.table(Db_yr=rep(x,nrow(tmp))))
			})
	df				<- rbindlist(tmp)
	setkey(df, Patient)	
	#compute proportions of patients in recent, undiagnosed by CD4 count, treated per year	
	df	<- df[,	{
				x	<- data.table(Patient, AnyPos_T1, isAcute,  AnyT_T1, Db_yr, DateEnd_yr)
				#x	<- subset(df, Patient=="M10833", select=c(Patient, AnyPos_T1, isAcute,  AnyT_T1, Db_yr, DateEnd_yr))
				z	<- merge(x, df.immu, by="Patient")
				z	<- subset(z, floor(PosCD4)==Db_yr)
				tmp	<- setdiff( seq.int(x[,floor(AnyPos_T1)[1]], x[,DateEnd_yr[1]]), unique(z[,Db_yr]) )	#years for which no CD4 count available
				if(length(tmp))
					x	<- rbind(z,data.table(Patient=x[,Patient[1]], AnyPos_T1= x[,AnyPos_T1[1]], isAcute=x[,isAcute[1]], AnyT_T1=x[,AnyT_T1[1]], Db_yr=tmp, DateEnd_yr=x[,DateEnd_yr[1]],   PosCD4= tmp+.5, CD4=NA)) 
				else
					x	<- z
				
				tmp<- x[, 	{									
							tmp		<- min(Db_yr[1]+1,AnyT_T1[1])-max(Db_yr[1],AnyPos_T1[1])
							Acute_p	<- ifelse(!is.na(isAcute[1]) && isAcute[1]=="Yes" && floor(AnyPos_T1[1])==Db_yr, tmp, 0)
							AnyT_p	<- min(1,max(0,Db_yr[1]+1-AnyT_T1[1]))
							CD4_p	<- ifelse((is.na(isAcute[1]) || isAcute[1]!="Yes") && floor(AnyT_T1[1])>=Db_yr, tmp, 0)
							#print(Acute_p); print(AnyT_p); print(CD4_p)
							tmp2	<- which(PosCD4<=AnyT_T1)
							CD4_med	<- ifelse(length(tmp2),median(CD4[tmp2]), NA_real_)
							#print(CD4_yr)
							list(Patient=Patient[1], Acute_p=Acute_p, CD4_med=CD4_med, AnyT_p=AnyT_p, NotYetT_p=CD4_p )
						},by=Db_yr]
				tmp
			},by=Patient]
	df[,CD4_medbin:= my.aggregate(df[,CD4_med], c(-Inf, 200, 350, 500, Inf))]	
	set(df, which(is.na(df[,CD4_medbin])), "CD4_medbin", "unknown")
	#take sum for each stratification of interest
	t.living	<- df[,	{				
				CD4_b200_n<- which(CD4_medbin=="-Inf,200")
				CD4_200_n<- which(CD4_medbin=="200,350")
				CD4_350_n<- which(CD4_medbin=="350,500")
				CD4_500_n<- which(CD4_medbin=="500,Inf")
				CD4_NA_n<- which(CD4_medbin=="unknown")
				list(Acute_n=sum(Acute_p), T_n=sum(AnyT_p), CD4_b200_n=sum(NotYetT_p[CD4_b200_n]), CD4_200_n=sum(NotYetT_p[CD4_200_n]), CD4_350_n=sum(NotYetT_p[CD4_350_n]), CD4_500_n=sum(NotYetT_p[CD4_500_n]), CD4_NA_n=sum(NotYetT_p[CD4_NA_n]))	
			},by=Db_yr]
	setkey(t.living, Db_yr)
	tmp					<- t.living[,Db_yr]
	t.living			<- t( as.matrix(subset(t.living, select=c(Acute_n, CD4_b200_n, CD4_200_n, CD4_350_n, CD4_500_n, CD4_NA_n, T_n))) )
	colnames(t.living)	<- tmp
	rownames(t.living)	<- c("Acute","<200","200-349","350-499",">=500","unknown","Treated")
	
	if(!is.null(plot.file))
	{
		cols				<- brewer.pal(nrow(t.living),"Set1")
		if(verbose)	cat(paste("\nplot file prop to ",plot.file))
		pdf(file=plot.file, width=5, height=4)
		par(mar=c(4,6,0.5,0.5))
		my.barplot.table(t.living, 0.4, "year", plot.ylab, cols, x.baseline=0, xax=as.numeric(colnames(t.living)), legend.loc="topleft")
		dev.off()
	}
	p.living			<- t.living / matrix(rep(apply(t.living,2,sum), each=nrow(t.living)), nrow(t.living), ncol(t.living))
	if(!is.null(plot.file.p))
	{
		if(verbose)	cat(paste("\nplot file prop to ",plot.file.p))
		pdf(file=plot.file.p, width=5, height=4)
		par(mar=c(4,6,0.5,0.5))
		cols				<- brewer.pal(nrow(p.living),"Set1")
		my.barplot.table(p.living, 0.4, "year", paste("%",plot.ylab), cols, x.baseline=0, xax=as.numeric(colnames(p.living)), legend.loc=NULL)
		dev.off()
	}
	list(t.living=t.living,p.living=p.living)
}
######################################################################################
hivc.db.getplot.livingbyexposure<- function(df, plot.file, plot.file.p, plot.ylab, db.endtime=2013.3, db.diff.lastcontact2died=0.5, db.diff.lastcontact2now= 2.3, verbose=1)
{
	#simplify risk group 
	set(df,which(df[,Trm=="SXCH"]),"Trm","OTH")
	set(df,which(df[,Trm=="PREG"]),"Trm","OTH")
	set(df,which(df[,Trm=="NEEACC"]),"Trm","OTH")
	set(df,which(df[,Trm=="BLOOD"]),"Trm","OTH")
	set(df,which(df[,Trm=="HETfa"]),"Trm","HET")	
	set(df,which(is.na(df[,Trm])),"Trm","unknown")	
	set(df, NULL, "Trm", factor(df[,Trm]))
	#
	set(df, NULL, "AnyPos_T1", hivc.db.Date2numeric(df[,AnyPos_T1]))
	set(df, NULL, "DateDied", hivc.db.Date2numeric(df[,DateDied]))
	set(df, NULL, "DateLastContact", hivc.db.Date2numeric(df[,DateLastContact]))
	df[,DateEnd:= DateDied]
	#compute DateEnd, either death or lost contact
	tmp				<- which( df[, (DateDied-DateLastContact)>db.diff.lastcontact2died] )
	if(verbose)	cat(paste("\n(DateDied-DateLastContact)>db.diff.lastcontact2died for n=",length(tmp)))
	set(df, tmp, "DateEnd", df[tmp,DateLastContact])
	tmp				<- which(df[, is.na(DateDied) & (db.endtime-DateLastContact)>db.diff.lastcontact2now])
	if(verbose)	cat(paste("\n(db.endtime-DateLastContact)>db.diff.lastcontact2now for n=",length(tmp)))
	set(df, tmp, "DateEnd", df[tmp,DateLastContact])
	tmp				<- which(is.na(df[,DateLastContact]))
	if(verbose)	cat(paste("\nnumber patients with missing DateLastContact n=",length(tmp)))
	tmp				<- which(is.na(df[,DateEnd]))
	if(verbose)	cat(paste("\nnumber patients still alive n=",length(tmp)))
	set(df, tmp,"DateEnd",db.endtime)
	df				<- subset(df,select=c(Patient, AnyPos_T1, Trm, DateEnd))	
	df[, AnyPos_yr:= floor(AnyPos_T1)]
	df[, DateEnd_yr:= floor(DateEnd)]
	#compute patients alive in year Db_yr 
	tmp				<- seq.int(min(df[, AnyPos_yr]),max(df[, DateEnd_yr]))
	tmp				<- lapply(tmp, function(x)
			{				
				tmp	<- subset(df, AnyPos_yr<=x & DateEnd_yr>=x)
				cbind(tmp,data.table(Db_yr=rep(x,nrow(tmp))))
			})
	df				<- rbindlist(tmp)
	#
	t.living			<- table(df[,Db_yr,Trm])
	rownames(t.living)[rownames(t.living)=="BI"]<- "Bi"
	rownames(t.living)[rownames(t.living)=="HET"]<- "Het"
	rownames(t.living)[rownames(t.living)=="IDU"]<- "DU"
	rownames(t.living)[rownames(t.living)=="OTH"]<- "other"	
	t.living	<- t.living[c("MSM","Bi",setdiff(rownames(t.living),c("MSM","Bi","unknown")),"unknown"),]
	#	
	if(!is.null(plot.file))
	{
		cols				<- brewer.pal(nrow(t.living),"Set2")
		if(verbose)	cat(paste("\nplot file prop to ",plot.file))
		pdf(file=plot.file, width=5, height=4)
		par(mar=c(4,6,0.5,0.5))
		my.barplot.table(t.living, 0.4, "year", plot.ylab, cols, x.baseline=0, xax=as.numeric(colnames(t.living)), legend.loc="topleft")
		dev.off()
	}		
	p.living			<- t.living / matrix(rep(apply(t.living,2,sum), each=nrow(t.living)), nrow(t.living), ncol(t.living))
	if(!is.null(plot.file.p))
	{
		if(verbose)	cat(paste("\nplot file prop to ",plot.file.p))
		pdf(file=plot.file.p, width=5, height=4)
		par(mar=c(4,6,0.5,0.5))
		cols				<- brewer.pal(nrow(p.living),"Set2")
		my.barplot.table(p.living, 0.4, "year", paste("%",plot.ylab), cols, x.baseline=0, xax=as.numeric(colnames(p.living)), legend.loc=NULL)
		dev.off()
	}
	list(t.living=t.living,p.living=p.living)
}
######################################################################################
hivc.db.Date2numeric<- function( x )
{
	x	<- as.POSIXlt(x)
	tmp	<- x$year + 1900
	x	<- tmp + round( x$yday / ifelse((tmp%%4==0 & tmp%%100!=0) | tmp%%400==0,366,365), d=3 )
	x	
}
######################################################################################
#	get BEAST taxon labels:		cluster	FASTASampleCode	NegT	AnyPosT	SeqT  -> turn date into numerical format
hivc.beast.addBEASTLabel<- function( df, df.resetTipDate=NA )
{
	if(is.na(df.resetTipDate))
		df[,dummy:=NA]
	else if(df.resetTipDate=="LdTd")
		df	<- merge(df, df[,list(dummy=max(AnyPos_T1)),by="cluster"])
	else if(df.resetTipDate=="LsTd")
		df	<- merge(df, df[,list(dummy=max(PosSeqT)),by="cluster"])
	else if(df.resetTipDate=="UmTd")
		df[,dummy:=max(df[,PosSeqT])]
	else
		stop("Unexpected df.resetTipDate")
	df	<- merge(df, df[,list(PosSeqTadj= max(PosSeqT, dummy, na.rm=T)),by="FASTASampleCode"], by="FASTASampleCode")
	
	tmp	<- df[,		{
						z	<- as.POSIXlt(c(NegT, AnyPos_T1, PosSeqT, PosSeqTadj))
						tmp	<- z$year + 1900
						z	<- tmp + round( z$yday / ifelse((tmp%%4==0 & tmp%%100!=0) | tmp%%400==0,366,365), d=3 )												
						list(BEASTlabel= paste(c(cluster, z, FASTASampleCode), collapse='_', sep=''))
					}, by="FASTASampleCode"]
	df	<- merge(df, tmp, by="FASTASampleCode")
	subset(df,select=-which(colnames(df)=="dummy"))	
}
######################################################################################
#	pool clusters into sets containing roughly 'pool.ntip' sequences
hivc.beast.poolclusters<- function(cluphy.df, pool.ntip= 130, pool.includealwaysbeforeyear=NA, verbose=1)
{	
	df			<- cluphy.df[, list(clu.ntip=clu.ntip[1], clu.AnyPos_T1=clu.AnyPos_T1[1]), by="cluster"]	
	if(!is.na(pool.includealwaysbeforeyear))
	{
		if(verbose) cat(paste("\nalways include clusters starting before ",pool.includealwaysbeforeyear,"and then pool evenly across clu.AnyPos_T1"))
		df.always	<- subset(df,as.POSIXlt(clu.AnyPos_T1)$year<(pool.includealwaysbeforeyear-1900))
		df			<- subset(df,as.POSIXlt(clu.AnyPos_T1)$year>=(pool.includealwaysbeforeyear-1900))
		pool.ntip	<- pool.ntip - sum(df.always[,clu.ntip])		
		setkey(df, clu.AnyPos_T1)
		pool.n		<- ceiling( sum( df[,clu.ntip] ) / pool.ntip )
		tmp			<- lapply( seq_len(pool.n), function(x)	seq.int(x,nrow(df),by=pool.n) )		
		pool.df		<- lapply(seq_along(tmp), function(i) merge(subset(rbind(df.always, df[tmp[[i]],]), select=cluster), cluphy.df, by="cluster") )		
	}
	else
	{
		if(verbose) cat(paste("\npool evenly across clu.AnyPos_T1"))
		setkey(df, clu.AnyPos_T1)
		pool.n	<- ceiling( sum( df[,clu.ntip] ) / pool.ntip )
		tmp		<- lapply( seq_len(pool.n), function(x)	seq.int(x,nrow(df),by=pool.n) )
		pool.df	<- lapply(seq_along(tmp), function(i) merge(subset(df[tmp[[i]],], select=cluster), cluphy.df, by="cluster") )		
	}
	if(verbose) cat(paste("\nnumber of pools is n=",pool.n))		
	if(verbose) cat(paste("\nnumber of seq in pools is n=",paste( sapply(pool.df, nrow), sep='', collapse=', ' )))
	list(pool.df=pool.df, pool.ntip=pool.ntip)
}	
######################################################################################
#	add taxa and alignment in bxml from BEASTlabels in df and alignment in seq.PROT.RT
#	beast.label.datepos= 4; beast.label.sep= '_'; beast.date.direction= "forwards"; beast.date.units= "years"; beast.alignment.dataType= "nucleotide"; xml.resetTipDate2LastDiag=1
hivc.beast.add.seq<- function(bxml, df, seq.PROT.RT, beast.label.datepos= 4, beast.label.sep= '_', beast.date.direction= "forwards", beast.date.units= "years", beast.alignment.dataType= "nucleotide", verbose=1)
{			
	bxml.beast	<- getNodeSet(bxml, "//beast")[[1]]
	
	#	add sequence taxa
	tmp.label	<- df[,BEASTlabel]	
	if(verbose)	cat(paste("\nsetting tip date to time at pos x in label, x=", beast.label.datepos))
	tmp.date	<- sapply( strsplit(tmp.label, beast.label.sep, fixed=1), function(x) x[beast.label.datepos] )
	dummy		<- newXMLCommentNode(text="The list of taxa to be analysed (can also include dates/ages).", parent=bxml.beast, doc=bxml, addFinalizer=T)
	dummy		<- newXMLCommentNode(text=paste("ntax=",nrow(df),sep=''), parent=bxml.beast, doc=bxml, addFinalizer=T)	
	seqtaxa		<- newXMLNode("taxa", attrs= list(id="taxa"), parent=bxml.beast, doc=bxml, addFinalizer=T)
	dummy		<- lapply(seq_along(tmp.label), function(i)
			{
				taxon	<- newXMLNode("taxon", attrs= list(id=tmp.label[i]), parent=seqtaxa, doc=bxml, addFinalizer=T )
				dummy	<- newXMLNode("date", attrs= list(value=tmp.date[i], direction=beast.date.direction, units=beast.date.units), parent=taxon, doc=bxml, addFinalizer=T )
				taxon
			})	
	if(verbose)	cat(paste("\nadded new seq taxa, n=", xmlSize(seqtaxa)))
	#	add alignment	
	dummy		<- newXMLCommentNode(text="The sequence alignment (each sequence refers to a taxon above).", parent=bxml.beast, doc=bxml, addFinalizer=T)
	dummy		<- newXMLCommentNode(text=paste("ntax=",nrow(df)," nchar=",ncol(seq.PROT.RT),sep=''), parent=bxml.beast, doc=bxml, addFinalizer=T)	
	seqalign	<- newXMLNode("alignment", attrs= list(id="alignment", dataType=beast.alignment.dataType), parent=bxml.beast, doc=bxml, addFinalizer=T)
	dummy		<- lapply( seq_len(nrow(df)), function(i)
			{
				seq		<- newXMLNode("sequence", parent=seqalign, doc=bxml, addFinalizer=T)
				dummy	<- newXMLNode("taxon", attrs= list(idref= df[i, BEASTlabel]), parent=seq, doc=bxml, addFinalizer=T)
				tmp		<- which( rownames(seq.PROT.RT)==df[i, FASTASampleCode] )
				if(length(tmp)!=1)	stop("unexpected row in seq.PROT.RT selected")
				tmp		<- paste(as.character(seq.PROT.RT[tmp,])[1,],collapse='',sep='')						
				dummy	<- newXMLTextNode(text=tmp, parent=seq, doc=bxml,addFinalizer=T)
				seq
			})	
	if(verbose)	cat(paste("\nadded new alignments, n=", xmlSize(seqalign)))	
	bxml
}	
######################################################################################
#	Read a BEAST treeannotator file. The node.label is set to a data.table that contains the SIMMAP annotation for the interior nodes in the newick tree.
hivc.treeannotator.read<- function(file, add.to.tiplabel=NA, rate.multiplier=NA, round.digit=NA, verbose=1) 
{
	require(data.table)
	require(ape)
	ph 			<- read.nexus(file)
	if(verbose)	cat(paste("\nReading BEAST file",file,sep=''))
	X 			<- scan(file = file, what = "", sep = "\n", quiet = TRUE)	
	#	read annotations of node in order as they appear in 'file'
	tab			<- X[grep("tree TREE1[[:space:]]+=", X)]
	tab 		<- gsub("tree TREE1[[:space:]]+= \\[&R\\] ", "", tab)
	tab 		<- unlist(strsplit(tab, "\\["))[-1]
	tab 		<- gsub("&|;|\\]", "", tab)
	tab 		<- gsub(":.+$", "", tab)
	tab 		<- lapply(tab, function(x) unlist(strsplit(x, ","))	)	
	tab			<- lapply(tab, function(x)
			{
				x			<- gsub('%','',x,fixed=1)
				ind 		<- grep("[{]", x)
				names 		<- gsub("=.+$", "", x[ind])
				x[ind] 		<- gsub("[{]", "", x[ind])
				x[ind] 		<- gsub("=", "_MIN=", x[ind])
				x[ind + 1] 	<- gsub("[}]", "", x[ind + 1])
				x[ind + 1] 	<- paste(paste(names, "MAX=", sep = "_"), x[ind + 1],sep='')
				x
			})
	colnames	<- unique(gsub("=.+$", "", unlist(tab)))
	if(verbose)	cat(paste("\nFound BEAST variables ",paste(colnames,collapse=', '),sep=''))
	tab			<- c(list(paste(colnames,-1,sep='=')), tab)		#rbindlist bug fix
	tab			<- lapply(tab, function(x)
			{
				ans									<- rep(NA, length(colnames))
				names(ans)							<- colnames				
				x									<- strsplit(x,'=',fixed=1)
				ans[ sapply(x, function(z) z[1]) ]	<- sapply(x, function(z) z[2])
				ans									<- paste(apply( rbind( names(ans), ans ), 2, function(z) paste(z,collapse='=',sep='')),collapse=',')
				eval(parse(text=paste("data.table(",ans,")",sep='')))				
			})
	df.beast	<- rbindlist(tab)[-1,]
	if(!is.na(rate.multiplier))
	{
		tmp	<- colnames(df.beast)[grepl("rate",colnames(df.beast))]
		sapply(tmp, function(x)		set(df.beast, NULL, x, as.numeric(unlist(df.beast[,x,with=F]))*rate.multiplier)			)		
	}
	if(!any(is.na(round.digit)))
	{
		if(length(round.digit)!=ncol(df.beast))
			round.digit<- rep(round.digit[1], ncol(df.beast))
		tmp	<- colnames(df.beast)
		sapply(seq_along(tmp), function(i)  
				{
					if(class(df.beast[[i]])=="numeric")			
						set(df.beast, NULL, tmp[i], round(as.numeric(unlist(df.beast[,tmp[i],with=F])), d=round.digit[i]))				
				})
	}
	
	tmp			<- length(which(df.beast[,!is.na(posterior)]))
	if(verbose)	cat(paste("\nFound annotated nodes, n=", tmp))
	if(verbose)	cat(paste("\nFound annotated tips, n=", nrow(df.beast)-tmp))	
	#	determine node index for 'df.beast':
	#
	#	- delete SIMMAP information from 'X'	
	LEFT 		<- grep("\\[", X)
	RIGHT 		<- grep("\\]", X)
	if (length(LEFT)) 
	{
		w 	<- LEFT == RIGHT
		if (any(w)) 
		{
			s 		<- LEFT[w]
			X[s] 	<- gsub("\\[[^]]*\\]", "", X[s])
		}
		w <- !w
		if(any(w)) 
		{
			s 		<- LEFT[w]
			X[s] 	<- gsub("\\[.*", "", X[s])
			sb	 	<- RIGHT[w]
			X[sb] 	<- gsub(".*\\]", "", X[sb])
			if(any(s < sb - 1)) 
				X 	<- X[-unlist(mapply(":", (s + 1), (sb - 1)))]
		}
	}	
	#	- read tree block
	endblock 			<- grep("END;|ENDBLOCK;", X, ignore.case = TRUE)
	semico 				<- grep(";", X)
	i1 					<- grep("BEGIN TREES;", X, ignore.case = TRUE)
	i2 					<- grep("TRANSLATE", X, ignore.case = TRUE)
	tree 				<- X[(semico[semico > i2][1] + 1):(endblock[endblock > i1][1] - 1)]
	tree 				<- gsub("^.*= *", "", tree)
	tree				<- substr(tree, 1, nchar(tree)-1)	 
	# 	- For each node, add a dummy node label that is the index in 'df.beast'	
	tmp					<- unlist(strsplit(tree, ":"))
	interiorm1			<- which( df.beast[,!is.na(posterior)] )
	tmp					<- sapply(seq_along(tmp), function(i) 			ifelse(i %in% interiorm1, 	paste(tmp[i],i,sep=''),	tmp[i])				)
	tmp					<- paste(paste(tmp, collapse=':'),';',sep='')
	#	- read this newick string and determine the node index in 'df.beast'
	tmp					<- read.tree(text=tmp)
	ph[["node.label"]]	<- cbind(data.table(node=Ntip(ph) + seq_len(Nnode(ph))), df.beast[as.numeric( tmp$node.label ),])
	setkey(ph[["node.label"]], node)
		
	if(!any(is.na(add.to.tiplabel)))
	{
		if(length(intersect(add.to.tiplabel,colnames(df.beast)))!=length(add.to.tiplabel))
			stop("Cannot find add.to.tiplabel")
		tmp				<- as.matrix(subset( df.beast[as.numeric( tmp$tip.label ),], select=add.to.tiplabel, with=F))
		tmp				<- cbind(ph$tip.label, tmp)
		ph$tip.label	<- apply(tmp, 1, function(x) paste(x,collapse='_'))		
	}	
	
	ph
}
######################################################################################
#	write nexus file for all sequences specified in df. assumes df has BEASTlabel. assumes seq.DNAbin.matrix and ph contain FASTASampleCode in df.
hivc.beast.writeNexus4Beauti<- function( seq.DNAbin.matrix, df, ph=NULL, file=NULL )
{
	# 	select sequences and relabel					
	seq.DNAbin.matrix			<- seq.DNAbin.matrix[ df[,FASTASampleCode], ]	
	rownames(seq.DNAbin.matrix)	<- df[,BEASTlabel]
	#	generate starting tree and relabel
	if(!is.null(ph))
	{			
		tmp					<- setdiff( ph$tip.label, df[,FASTASampleCode] )
		tmp					<- match( tmp, ph$tip.label)
		ph.start			<- drop.tip(ph, tmp)
		ph.start$node.label	<- NULL
		setkey(df, FASTASampleCode)
		ph.start$tip.label	<- df[ph.start$tip.label,][,BEASTlabel]		
	}
	#	produce nexus text			
	ans<- hivc.seq.write.dna.nexus(seq.DNAbin.matrix, ph=ph.start, file=file)
	ans
}
######################################################################################
#	read tip stem samples after burn in from log file and return upper left points of a histogram of the Monte Carlo sample
hivc.beast.read.log2tstem<- function(file.log, file.xml, beastlabel.idx.samplecode=6, burn.in= 5e6, breaks.n= 30, verbose=0)
{
	require(data.table)
	require(XML)	
	if(verbose)	cat(paste("\nReading file ",file.log))
	df.log	<- read.delim2(file.log, header = TRUE, sep = "\t", quote="\"", dec=".", fill = TRUE, comment.char="#")
	df.log	<- as.data.table(df.log)		
	if(verbose)	cat(paste("\nRead log file with ncol=", ncol(df.log)))
	tmp		<- c("state", colnames(df.log)[ grepl("tstem", colnames(df.log)) ] )
	if(length(tmp)==1)
	{
		if(verbose)	cat(paste("\nNo tstem found, return NA"))
		return( data.table( FASTASampleCode=NA, tstem= NA, density= NA ) )
	}
	df.log	<- subset( df.log, state>burn.in, select=tmp )
	if(verbose)	cat(paste("\nFound tstem data for n=", ncol(df.log)-1))
	#	translate tips back into FASTASampleCode
	if(verbose)	cat(paste("\nReading file ",file.xml))
	bxml				<- xmlTreeParse(file.xml, useInternalNodes=TRUE, addFinalizer = TRUE)
	tmp					<- regexpr("tip[0-9]+",colnames(df.log))
	if(!length(which(tmp>0)))
	{
		if(verbose)	cat(paste("\nNo tip tstem found, return NA"))
		return( data.table( FASTASampleCode=NA, tstem= NA, density= NA ) )
	}	
	log.tips			<- sapply(seq_along(tmp)[tmp>0], function(i)	substr(colnames(df.log)[i],tmp[i],tmp[i]+attr(tmp,"match.length")[i]-1)		)
	log.FASTASampleCode	<- sapply(log.tips, function(x) unlist( xpathApply(bxml, paste("//taxa[@id='",x,"']/taxon",sep=''), xmlGetAttr, "idref") )	)
	log.FASTASampleCode	<- sapply( strsplit(log.FASTASampleCode, '_'), function(x) x[beastlabel.idx.samplecode])		
	setnames(df.log, colnames(df.log)[tmp>0], log.FASTASampleCode)
	#	compute histograms for each tip stem
	ans		<- lapply(log.FASTASampleCode, function(x)
			{
				#x<- "R11-11357"
				tmp	<- hist( as.numeric(unlist( df.log[, x, with=0] )), breaks=breaks.n, plot=0 )
				data.table( FASTASampleCode=x, tstem= tmp$breaks, density= c(tmp$density,0) )										
			})
	if(verbose)	cat(paste("\nProcessed tstem data for tips n=", length(log.FASTASampleCode)))
	ans	<- rbindlist(ans)
}
######################################################################################
hivc.beast2.add.alignment<- function(bxml, seq.PROT.RT, df, beast2.spec, verbose=1)
{
	bxml.beast	<- getNodeSet(bxml, "//beast")[[1]]
	dummy		<- newXMLCommentNode(text="The sequence alignment.", parent=bxml.beast, doc=bxml, addFinalizer=T)
	dummy		<- newXMLCommentNode(text=paste("ntax=",nrow(df)," nchar=",ncol(seq.PROT.RT),sep=''), parent=bxml.beast, doc=bxml, addFinalizer=T)	
	seqalign	<- newXMLNode("data", attrs= list(id=beast2.spec$alignment.id, name="alignment", dataType=beast2.spec$alignment.dataType, missing=beast2.spec$alignment.missing), parent=bxml.beast, doc=bxml, addFinalizer=T)
	dummy		<- lapply( seq_len(nrow(df)), function(i)
			{
				tmp		<- which( rownames(seq.PROT.RT)==df[i, FASTASampleCode] )
				if(length(tmp)!=1)	stop("unexpected row in seq.PROT.RT selected")
				tmp		<- paste(as.character(seq.PROT.RT[tmp,])[1,],collapse='',sep='')										
				seq		<- newXMLNode("sequence", attrs=list(id= df[i, BEASTlabel], taxon=df[i, BEASTlabel], value=tmp), parent=seqalign, doc=bxml, addFinalizer=T)
				seq
			})	
	if(verbose)	cat(paste("\nadded new alignment with sequences, n=", xmlSize(seqalign)))
	bxml
}
######################################################################################
hivc.beast2.add.datetrait<- function(bxml, df, beast2.spec, verbose=1)	
{			
	tmp			<- df[,BEASTlabel]	
	if(verbose)	cat(paste("\nsetting tip date to time at pos x in label, x=", beast2.spec$beast.label.datepos))
	tmp			<- sapply( strsplit(tmp, beast2.spec$beast.label.sep, fixed=1), function(x) paste(paste(x,collapse=beast2.spec$beast.label.sep,sep=''),'=',x[beast2.spec$beast.label.datepos],sep='') )
	
	bxml.beast	<- getNodeSet(bxml, "//beast")[[1]]
	dummy		<- newXMLCommentNode(text="The tip dates.", parent=bxml.beast, doc=bxml, addFinalizer=T)
	dummy		<- newXMLCommentNode(text=paste("ntax=",nrow(df),sep=''), parent=bxml.beast, doc=bxml, addFinalizer=T)
	
	seqtrait	<- newXMLNode("trait", attrs= list(id=paste("dateTrait.t",beast2.spec$alignment.id,sep=':'), 
					spec= beast2.spec$datetrait.spec,
					units= beast2.spec$datetrait.units,
					traitname= beast2.spec$datetrait.traitname, 
					value=paste(tmp, collapse=", ",sep='')), parent=bxml.beast, doc=bxml, addFinalizer=T)
	if(verbose)	cat(paste("\nadded new date trait for sequences, n=", length(tmp)))
	#define taxonset in here
	seqtaxa		<- newXMLNode("taxa", attrs= list(id=paste("TaxonSet.t",beast2.spec$alignment.id,sep=':'), spec=beast2.spec$datetrait.taxa.spec), parent=seqtrait, doc=bxml, addFinalizer=T)
	dummy		<- newXMLNode("alignment", attrs= list(idref=beast2.spec$alignment.id), parent=seqtaxa, doc=bxml, addFinalizer=T)
	bxml
}
######################################################################################
hivc.beast2.add.satree<- function(bxml, beast2.spec, verbose=verbose)
{
	bxml.beast	<- getNodeSet(bxml, "//beast")[[1]]
	beast.tree	<- newXMLNode("tree", attrs= list(	id=beast2.spec$tree.id, 													 
													trait=paste('@',"dateTrait.t:",beast2.spec$alignment.id,sep=''),
													nodetype=beast2.spec$sasky.tree.nodetype,
													clusterType=beast2.spec$sasky.tree.clusterType,
													spec=beast2.spec$sasky.tree.spec,
													taxa=paste('@',beast2.spec$alignment.id,sep='')), parent=bxml.beast, doc=bxml, addFinalizer=T)	
	if(verbose)	cat(paste("\nadded SA trees for taxonsets, n=", length(beast2.spec$tree.taxonset)))
	bxml
}
######################################################################################
hivc.beast2.add.tree<- function(bxml, beast2.spec, verbose=verbose)
{
	bxml.beast	<- getNodeSet(bxml, "//beast")[[1]]
	beast.tree	<- newXMLNode("tree", attrs= list(id=beast2.spec$tree.id, name=beast2.spec$tree.id, trait=paste('@',"dateTrait.t:",beast2.spec$alignment.id,sep='')), parent=bxml.beast, doc=bxml, addFinalizer=T)
	dummy		<- newXMLNode("taxonset", attrs= list(idref=paste("TaxonSet",beast2.spec$tree.taxonset,sep='.t:')), parent=beast.tree, doc=bxml, addFinalizer=T)
	if(verbose)	cat(paste("\nadded trees for taxonsets, n=", length(beast2.spec$tree.taxonset)))
	bxml
}
######################################################################################
hivc.beast2.add.treemodel.bdsky<- function(bxml, beast2.spec, verbose=1)
{
	bxml.beast		<- getNodeSet(bxml, "//beast")[[1]]
	beast.treemodel	<- newXMLNode("BirthDeathSkylineModel", attrs= list(	id=beast2.spec$treemodel.id, 
					name=beast2.spec$treemodel.id, 
					tree=paste('@',beast2.spec$tree.id,sep=''),
					spec=beast2.spec$bdsky.spec,
					intervalNumber=as.character(beast2.spec$bdsky.intervalNumber)), parent=bxml.beast, doc=bxml, addFinalizer=T)
	dummy			<- newXMLNode("parameter", attrs= list(	id=beast2.spec$bdsky.sprop.id, 
					name="samplingProportion",
					dimension=beast2.spec$bdsky.intervalNumber,
					lower=as.character(beast2.spec$bdsky.sprop.lower),
					upper=as.character(beast2.spec$bdsky.sprop.upper),
					value=paste(beast2.spec$bdsky.sprop.value, collapse=' ')), parent=beast.treemodel, doc=bxml, addFinalizer=T)
	dummy			<- newXMLNode("parameter", attrs= list(	id=beast2.spec$bdsky.R0.id, 
					name="R0",
					dimension=beast2.spec$bdsky.intervalNumber,
					lower=as.character(beast2.spec$bdsky.R0.lower),
					upper=as.character(beast2.spec$bdsky.R0.upper),
					value=paste(beast2.spec$bdsky.R0.value, collapse=' ')), parent=beast.treemodel, doc=bxml, addFinalizer=T)											
	dummy			<- newXMLNode("parameter", attrs= list(	id=beast2.spec$bdsky.notInf.id, 
					name="becomeUninfectiousRate",
					dimension=beast2.spec$bdsky.intervalNumber,
					lower=as.character(beast2.spec$bdsky.notInf.lower),
					upper=as.character(beast2.spec$bdsky.notInf.upper),
					value=paste(beast2.spec$bdsky.notInf.value, collapse=' ')), parent=beast.treemodel, doc=bxml, addFinalizer=T)
	dummy			<- newXMLNode("parameter", attrs= list(	id=beast2.spec$bdsky.origin.id, 
					name="origin",
					lower=as.character(beast2.spec$bdsky.origin.lower),
					upper=as.character(beast2.spec$bdsky.origin.upper),
					value=paste(beast2.spec$bdsky.origin.value)), parent=beast.treemodel, doc=bxml, addFinalizer=T)
	tmp				<- rep("false",4)
	tmp[which(!sapply(c(beast2.spec$bdsky.R0.changepoint.id[1],beast2.spec$bdsky.notInf.changepoint.id[1],beast2.spec$bdsky.sprop.changepoint.id[1]),is.null))]		<- "true"
	dummy			<- newXMLNode("reverseTimeArrays", attrs= list(	id=beast2.spec$bdsky.reverseTimeArrays.id, 
																	spec=beast2.spec$bdsky.reverseTimeArrays.spec,
																	value=paste(tmp, collapse=' ')), parent=beast.treemodel, doc=bxml, addFinalizer=T)
	if(!is.null(beast2.spec$bdsky.R0.changepoint.id[1]))												
		dummy		<- newXMLNode("parameter", attrs= list(	id=beast2.spec$bdsky.R0.changepoint.id, 
						name="birthRateChangeTimes",
						value=paste(beast2.spec$bdsky.R0.changepoint.value, collapse=' ')), parent=beast.treemodel, doc=bxml, addFinalizer=T)
	
	if(!is.null(beast2.spec$bdsky.notInf.changepoint.id[1]))												
		dummy		<- newXMLNode("parameter", attrs= list(	id=beast2.spec$bdsky.notInf.changepoint.id, 
						name="deathRateChangeTimes",
						value=paste(beast2.spec$bdsky.notInf.changepoint.value, collapse=' ')), parent=beast.treemodel, doc=bxml, addFinalizer=T)
	if(!is.null(beast2.spec$bdsky.sprop.changepoint.id[1]))												
		dummy		<- newXMLNode("parameter", attrs= list(	id=beast2.spec$bdsky.sprop.changepoint.id, 
						name="samplingRateChangeTimes",
						value=paste(beast2.spec$bdsky.sprop.changepoint.value, collapse=' ')), parent=beast.treemodel, doc=bxml, addFinalizer=T)
	if(verbose)	cat(paste("\nadded BDSKY tree models for taxonsets, n=", length(beast2.spec$tree.taxonset)))
	bxml	
}
######################################################################################
hivc.beast2.add.treemodel.sasky<- function(bxml, beast2.spec, verbose=1)
{
	bxml.beast		<- getNodeSet(bxml, "//beast")[[1]]
	beast.treemodel	<- newXMLNode("BirthDeathSkylineModel", attrs= list(	id=beast2.spec$treemodel.id, 
					name=beast2.spec$treemodel.id, 
					tree=paste('@',beast2.spec$tree.id,sep=''),
					spec=beast2.spec$sasky.spec,
					intervalNumber=as.character(beast2.spec$bdsky.intervalNumber)), parent=bxml.beast, doc=bxml, addFinalizer=T)
	dummy			<- newXMLNode("parameter", attrs= list(	id=beast2.spec$bdsky.sprop.id, 
					name="samplingProportion",
					dimension=beast2.spec$bdsky.intervalNumber,
					lower=as.character(beast2.spec$bdsky.sprop.lower),
					upper=as.character(beast2.spec$bdsky.sprop.upper),
					value=paste(beast2.spec$bdsky.sprop.value, collapse=' ')), parent=beast.treemodel, doc=bxml, addFinalizer=T)
	dummy			<- newXMLNode("parameter", attrs= list(	id=beast2.spec$bdsky.R0.id, 
					name="R0",
					dimension=beast2.spec$bdsky.intervalNumber,
					lower=as.character(beast2.spec$bdsky.R0.lower),
					upper=as.character(beast2.spec$bdsky.R0.upper),
					value=paste(beast2.spec$bdsky.R0.value, collapse=' ')), parent=beast.treemodel, doc=bxml, addFinalizer=T)											
	dummy			<- newXMLNode("parameter", attrs= list(	id=beast2.spec$bdsky.notInf.id, 
					name="becomeUninfectiousRate",
					dimension=beast2.spec$bdsky.intervalNumber,
					lower=as.character(beast2.spec$bdsky.notInf.lower),
					upper=as.character(beast2.spec$bdsky.notInf.upper),
					value=paste(beast2.spec$bdsky.notInf.value, collapse=' ')), parent=beast.treemodel, doc=bxml, addFinalizer=T)
	dummy			<- newXMLNode("parameter", attrs= list(	id=beast2.spec$sasky.r.id, 
					name="becomeNoninfectiousAfterSamplingProbability",
					dimension=beast2.spec$bdsky.intervalNumber,
					lower=as.character(beast2.spec$sasky.r.lower),
					upper=as.character(beast2.spec$sasky.r.upper),
					value=paste(beast2.spec$sasky.r.value, collapse=' ')), parent=beast.treemodel, doc=bxml, addFinalizer=T)	
	dummy			<- newXMLNode("parameter", attrs= list(	id=beast2.spec$bdsky.origin.id, 
					name="origin",
					lower=as.character(beast2.spec$bdsky.origin.lower),
					upper=as.character(beast2.spec$bdsky.origin.upper),
					value=paste(beast2.spec$bdsky.origin.value)), parent=beast.treemodel, doc=bxml, addFinalizer=T)
	tmp				<- rep("false",4)
	tmp[which(!sapply(c(beast2.spec$bdsky.R0.changepoint.id[1],beast2.spec$bdsky.notInf.changepoint.id[1],beast2.spec$bdsky.sprop.changepoint.id[1]),is.null))]		<- "true"
	dummy			<- newXMLNode("reverseTimeArrays", attrs= list(	id=beast2.spec$bdsky.reverseTimeArrays.id, 
					spec=beast2.spec$bdsky.reverseTimeArrays.spec,
					value=paste(tmp, collapse=' ')), parent=beast.treemodel, doc=bxml, addFinalizer=T)
	if(!is.null(beast2.spec$bdsky.R0.changepoint.id[1]))												
		dummy		<- newXMLNode("parameter", attrs= list(	id=beast2.spec$bdsky.R0.changepoint.id, 
						name="birthRateChangeTimes",
						value=paste(beast2.spec$bdsky.R0.changepoint.value, collapse=' ')), parent=beast.treemodel, doc=bxml, addFinalizer=T)	
	if(!is.null(beast2.spec$bdsky.notInf.changepoint.id[1]))												
		dummy		<- newXMLNode("parameter", attrs= list(	id=beast2.spec$bdsky.notInf.changepoint.id, 
						name="deathRateChangeTimes",
						value=paste(beast2.spec$bdsky.notInf.changepoint.value, collapse=' ')), parent=beast.treemodel, doc=bxml, addFinalizer=T)
	if(!is.null(beast2.spec$bdsky.sprop.changepoint.id[1]))												
		dummy		<- newXMLNode("parameter", attrs= list(	id=beast2.spec$bdsky.sprop.changepoint.id, 
						name="samplingRateChangeTimes",
						value=paste(beast2.spec$bdsky.sprop.changepoint.value, collapse=' ')), parent=beast.treemodel, doc=bxml, addFinalizer=T)
	if(verbose)	cat(paste("\nadded SASKY tree models for taxonsets, n=", length(beast2.spec$tree.taxonset)))
	bxml	
}
######################################################################################
hivc.beast2.init.xml<- function( beast2.spec=NULL, verbose=1)
{
	bxml		<- newXMLDoc(addFinalizer=T)
	bxml.beast	<- newXMLNode("beast", attrs=list(beautitemplate='HIVCLUST', beautistatus='', version="2.0", namespace=paste(beast2.spec$namespace,sep='',collapse=':')), doc=bxml, addFinalizer=T)
	dummy		<- newXMLNode("map", attrs= list(name="Beta"), parent=bxml.beast, doc=bxml, addFinalizer=T)
	dummy		<- newXMLTextNode(beast2.spec$map.Beta, parent=dummy, doc=bxml,addFinalizer=T)
	dummy		<- newXMLNode("map", attrs= list(name="ExcludablePrior"), parent=bxml.beast, doc=bxml, addFinalizer=T)
	dummy		<- newXMLTextNode(beast2.spec$map.ExcludablePrior, parent=dummy, doc=bxml,addFinalizer=T)
	dummy		<- newXMLNode("map", attrs= list(name="Exponential"), parent=bxml.beast, doc=bxml, addFinalizer=T)
	dummy		<- newXMLTextNode(beast2.spec$map.Exponential, parent=dummy, doc=bxml,addFinalizer=T)
	dummy		<- newXMLNode("map", attrs= list(name="InverseGamma"), parent=bxml.beast, doc=bxml, addFinalizer=T)
	dummy		<- newXMLTextNode(beast2.spec$map.InverseGamma, parent=dummy, doc=bxml,addFinalizer=T)
	dummy		<- newXMLNode("map", attrs= list(name="LogNormal"), parent=bxml.beast, doc=bxml, addFinalizer=T)
	dummy		<- newXMLTextNode(beast2.spec$map.LogNormal, parent=dummy, doc=bxml,addFinalizer=T)
	dummy		<- newXMLNode("map", attrs= list(name="Gamma"), parent=bxml.beast, doc=bxml, addFinalizer=T)
	dummy		<- newXMLTextNode(beast2.spec$map.Gamma, parent=dummy, doc=bxml,addFinalizer=T)
	dummy		<- newXMLNode("map", attrs= list(name="Uniform"), parent=bxml.beast, doc=bxml, addFinalizer=T)
	dummy		<- newXMLTextNode(beast2.spec$map.Uniform, parent=dummy, doc=bxml,addFinalizer=T)
	dummy		<- newXMLNode("map", attrs= list(name="prior"), parent=bxml.beast, doc=bxml, addFinalizer=T)
	dummy		<- newXMLTextNode(beast2.spec$map.prior, parent=dummy, doc=bxml,addFinalizer=T)
	dummy		<- newXMLNode("map", attrs= list(name="LaplaceDistribution"), parent=bxml.beast, doc=bxml, addFinalizer=T)
	dummy		<- newXMLTextNode(beast2.spec$map.Laplace, parent=dummy, doc=bxml,addFinalizer=T)
	dummy		<- newXMLNode("map", attrs= list(name="OneOnX"), parent=bxml.beast, doc=bxml, addFinalizer=T)
	dummy		<- newXMLTextNode(beast2.spec$map.OneOnX, parent=dummy, doc=bxml,addFinalizer=T)
	dummy		<- newXMLNode("map", attrs= list(name="Normal"), parent=bxml.beast, doc=bxml, addFinalizer=T)
	dummy		<- newXMLTextNode(beast2.spec$map.Normal, parent=dummy, doc=bxml,addFinalizer=T)
	bxml
}	
######################################################################################
hivc.beast2.get.prior<- function(args, parent, parentid, bxml)
{
	args	<- strsplit(args,'/')[[1]]
	if(args[1]=="Uniform")
		prior	<- newXMLNode("Uniform", attrs= list( id= paste('U',parentid,sep='-'), name="distr", lower=args[2], upper=args[3]), parent=parent, doc=bxml, addFinalizer=T)											
	else if(args[1]=="OneOnX")
		prior	<- newXMLNode("OneOnX", attrs= list( id= paste('OneOnX.',parentid,sep=''), name="distr"), parent=parent, doc=bxml, addFinalizer=T)
	else if(args[1]=="Exponential")
	{
		prior	<- newXMLNode("Exponential", attrs= list( id= paste('Exponential',parentid,sep='-'), name="distr", offset=args[3]), parent=parent, doc=bxml, addFinalizer=T)
		dummy	<- newXMLNode("parameter", attrs= list( id= paste('pExponential',parentid,sep='-'), name="mean", value=args[2], estimate="false"), parent=prior, doc=bxml, addFinalizer=T)
	}
	else if(args[1]=="Beta")
	{
		prior	<- newXMLNode("Beta", attrs= list( id= paste('Beta',parentid,sep='-'), name="distr", offset=args[4]), parent=parent, doc=bxml, addFinalizer=T)
		dummy	<- newXMLNode("parameter", attrs= list( id= paste('pBeta1',parentid,sep='-'), lower="0.0", upper="10.0", name="alpha", value=args[2], estimate="false"), parent=prior, doc=bxml, addFinalizer=T)
		dummy	<- newXMLNode("parameter", attrs= list( id= paste('pBeta2',parentid,sep='-'), lower="0.0", upper="10.0", name="beta", value=args[3], estimate="false"), parent=prior, doc=bxml, addFinalizer=T)
	}
	else if(args[1]=="Gamma")
	{
		prior	<- newXMLNode("Gamma", attrs= list( id= paste('Gamma',parentid,sep='-'), name="distr", offset=args[4]), parent=parent, doc=bxml, addFinalizer=T)
		dummy	<- newXMLNode("parameter", attrs= list( id= paste('pGamma1',parentid,sep='-'), name="alpha",  value=args[2], estimate="false"), parent=prior, doc=bxml, addFinalizer=T)
		dummy	<- newXMLNode("parameter", attrs= list( id= paste('pGamma2',parentid,sep='-'), name="beta",  value=args[3], estimate="false"), parent=prior, doc=bxml, addFinalizer=T)
	}
	else if(args[1]=="LogNormal")
	{
		prior	<- newXMLNode("LogNormal", attrs= list( id= paste('LogNormal',parentid,sep='-'), name="distr", offset=args[4], meanInRealSpace=args[5]), parent=parent, doc=bxml, addFinalizer=T)
		dummy	<- newXMLNode("parameter", attrs= list( id= paste('pLogNormal1',parentid,sep='-'), name="M",  value=args[2], estimate="false"), parent=prior, doc=bxml, addFinalizer=T)
		dummy	<- newXMLNode("parameter", attrs= list( id= paste('pLogNormal2',parentid,sep='-'), name="S",  value=args[3], estimate="false"), parent=prior, doc=bxml, addFinalizer=T)
	}
	else	stop("prior not implemented")	
	prior
}
######################################################################################
hivc.beast2.add.bdsky.serialpriors<- function(bxml, beast2.spec, verbose=1)
{
	bxml.prior	<- getNodeSet(bxml, "//*[@id='prior']")
	if(length(bxml.prior)!=1)	stop("unexpected length of bxml.prior")
	bxml.prior	<- bxml.prior[[1]]	
	dummy		<- newXMLCommentNode(text="start: serial BDSKY priors.", parent=bxml.prior, doc=bxml, addFinalizer=T)
	dummy		<- newXMLCommentNode(text="The BDSKY model prior.", parent=bxml.prior, doc=bxml, addFinalizer=T)
	tmp			<- newXMLNode("distribution", attrs= list(	id=paste("prior",beast2.spec$treemodel.id,sep='-'),spec=beast2.spec$compoundprior.spec), parent=bxml.prior, doc=bxml, addFinalizer=T)	
	dummy		<- newXMLNode("distribution", attrs= list(	idref=beast2.spec$treemodel.id), parent=tmp, doc=bxml, addFinalizer=T)
	dummy		<- newXMLCommentNode(text="The origin prior.", parent=bxml.prior, doc=bxml, addFinalizer=T)	
	tmp			<- newXMLNode("prior", attrs= list(	id=paste("sprior",beast2.spec$bdsky.origin.id,sep='-'),
													name="distribution",
													x= paste('@',beast2.spec$bdsky.origin.id,sep='')), parent=bxml.prior, doc=bxml, addFinalizer=T)	
	dummy		<- hivc.beast2.get.prior(beast2.spec$bdsky.origin.prior, tmp, xmlAttrs(tmp)["id"], bxml)
	dummy		<- newXMLCommentNode(text="The serial samplingProb priors.", parent=bxml.prior, doc=bxml, addFinalizer=T)
	dummy		<- lapply(seq_len(beast2.spec$bdsky.intervalNumber), function(serial.id)
			{
				tmp				<- rep(0,beast2.spec$bdsky.intervalNumber)
				tmp[serial.id]	<- 1
				serialwrapper	<- newXMLNode("distribution", attrs= list(	id=paste("sprior",serial.id,beast2.spec$bdsky.sprop.id,sep='-'),
								name="distribution",
								x= paste('@',beast2.spec$bdsky.sprop.id,sep=''),
								xInclude=paste(tmp, collapse=' ',sep=''), 					
								spec= beast2.spec$bdsky.prior.spec), parent=bxml.prior, doc=bxml, addFinalizer=T)
				dummy			<- hivc.beast2.get.prior(beast2.spec$bdsky.sprop.prior[serial.id], serialwrapper, xmlAttrs(serialwrapper)["id"], bxml)				
			})
	if(verbose)	cat(paste("\nadded serial samplingProb priors, n=", length(dummy)))
	dummy		<- newXMLCommentNode(text="The serial R0 priors.", parent=bxml.prior, doc=bxml, addFinalizer=T)
	dummy		<- lapply(seq_len(beast2.spec$bdsky.intervalNumber), function(serial.id)
			{
				tmp				<- rep(0,beast2.spec$bdsky.intervalNumber)
				tmp[serial.id]	<- 1
				serialwrapper	<- newXMLNode("distribution", attrs= list(	id=paste("sprior",serial.id,beast2.spec$bdsky.R0.id,sep='-'),
								name="distribution",
								x= paste('@',beast2.spec$bdsky.R0.id,sep=''),
								xInclude=paste(tmp, collapse=' ',sep=''), 					
								spec= beast2.spec$bdsky.prior.spec), parent=bxml.prior, doc=bxml, addFinalizer=T)
				dummy			<- hivc.beast2.get.prior(beast2.spec$bdsky.R0.prior[serial.id], serialwrapper, xmlAttrs(serialwrapper)["id"], bxml)				
			})												
	if(verbose)	cat(paste("\nadded serial R0 priors, n=", length(dummy)))
	dummy		<- newXMLCommentNode(text="The serial becomingUninfectious priors.", parent=bxml.prior, doc=bxml, addFinalizer=T)
	dummy		<- lapply(seq_len(beast2.spec$bdsky.intervalNumber), function(serial.id)
			{
				tmp				<- rep(0,beast2.spec$bdsky.intervalNumber)
				tmp[serial.id]	<- 1
				serialwrapper	<- newXMLNode("distribution", attrs= list(	id=paste("sprior",serial.id,beast2.spec$bdsky.notInf.id,sep='-'),
								name="distribution",
								x= paste('@',beast2.spec$bdsky.notInf.id,sep=''),
								xInclude=paste(tmp, collapse=' ',sep=''), 					
								spec= beast2.spec$bdsky.prior.spec), parent=bxml.prior, doc=bxml, addFinalizer=T)
				dummy			<- hivc.beast2.get.prior(beast2.spec$bdsky.notInf.prior[serial.id], serialwrapper, xmlAttrs(serialwrapper)["id"], bxml)				
			})	
	if(verbose)	cat(paste("\nadded serial becomingUninfectious priors, n=", length(dummy)))
	dummy		<- newXMLCommentNode(text="end: serial BDSKY priors.", parent=bxml.prior, doc=bxml, addFinalizer=T)
	bxml.prior
}
######################################################################################
hivc.beast2.add.sasky.serialpriors<- function(bxml, beast2.spec, verbose=1)
{
	bxml.prior	<- getNodeSet(bxml, "//*[@id='prior']")
	if(length(bxml.prior)!=1)	stop("unexpected length of bxml.prior")
	bxml.prior	<- bxml.prior[[1]]	
	dummy		<- newXMLCommentNode(text="start: serial SASKY priors.", parent=bxml.prior, doc=bxml, addFinalizer=T)
	dummy		<- newXMLCommentNode(text="The SASKY model prior.", parent=bxml.prior, doc=bxml, addFinalizer=T)
	tmp			<- newXMLNode("distribution", attrs= list(	id=paste("prior",beast2.spec$treemodel.id,sep='-'),spec=beast2.spec$compoundprior.spec), parent=bxml.prior, doc=bxml, addFinalizer=T)	
	dummy		<- newXMLNode("distribution", attrs= list(	idref=beast2.spec$treemodel.id), parent=tmp, doc=bxml, addFinalizer=T)
	dummy		<- newXMLCommentNode(text="The origin prior.", parent=bxml.prior, doc=bxml, addFinalizer=T)	
	tmp			<- newXMLNode("prior", attrs= list(	id=paste("sprior",beast2.spec$bdsky.origin.id,sep='-'),
					name="distribution",
					x= paste('@',beast2.spec$bdsky.origin.id,sep='')), parent=bxml.prior, doc=bxml, addFinalizer=T)	
	dummy		<- hivc.beast2.get.prior(beast2.spec$bdsky.origin.prior, tmp, xmlAttrs(tmp)["id"], bxml)
	dummy		<- newXMLCommentNode(text="The serial samplingProb priors.", parent=bxml.prior, doc=bxml, addFinalizer=T)
	dummy		<- lapply(seq_len(beast2.spec$bdsky.intervalNumber), function(serial.id)
			{
				tmp				<- rep(0,beast2.spec$bdsky.intervalNumber)
				tmp[serial.id]	<- 1
				serialwrapper	<- newXMLNode("distribution", attrs= list(	id=paste("sprior",serial.id,beast2.spec$bdsky.sprop.id,sep='-'),
								name="distribution",
								x= paste('@',beast2.spec$bdsky.sprop.id,sep=''),
								xInclude=paste(tmp, collapse=' ',sep=''), 					
								spec= beast2.spec$bdsky.prior.spec), parent=bxml.prior, doc=bxml, addFinalizer=T)
				dummy			<- hivc.beast2.get.prior(beast2.spec$bdsky.sprop.prior[serial.id], serialwrapper, xmlAttrs(serialwrapper)["id"], bxml)				
			})
	if(verbose)	cat(paste("\nadded serial samplingProb priors, n=", length(dummy)))
	dummy		<- newXMLCommentNode(text="The serial R0 priors.", parent=bxml.prior, doc=bxml, addFinalizer=T)
	dummy		<- lapply(seq_len(beast2.spec$bdsky.intervalNumber), function(serial.id)
			{
				tmp				<- rep(0,beast2.spec$bdsky.intervalNumber)
				tmp[serial.id]	<- 1
				serialwrapper	<- newXMLNode("distribution", attrs= list(	id=paste("sprior",serial.id,beast2.spec$bdsky.R0.id,sep='-'),
								name="distribution",
								x= paste('@',beast2.spec$bdsky.R0.id,sep=''),
								xInclude=paste(tmp, collapse=' ',sep=''), 					
								spec= beast2.spec$bdsky.prior.spec), parent=bxml.prior, doc=bxml, addFinalizer=T)
				dummy			<- hivc.beast2.get.prior(beast2.spec$bdsky.R0.prior[serial.id], serialwrapper, xmlAttrs(serialwrapper)["id"], bxml)				
			})												
	if(verbose)	cat(paste("\nadded serial R0 priors, n=", length(dummy)))
	dummy		<- newXMLCommentNode(text="The serial becomingUninfectious priors.", parent=bxml.prior, doc=bxml, addFinalizer=T)
	dummy		<- lapply(seq_len(beast2.spec$bdsky.intervalNumber), function(serial.id)
			{
				tmp				<- rep(0,beast2.spec$bdsky.intervalNumber)
				tmp[serial.id]	<- 1
				serialwrapper	<- newXMLNode("distribution", attrs= list(	id=paste("sprior",serial.id,beast2.spec$bdsky.notInf.id,sep='-'),
								name="distribution",
								x= paste('@',beast2.spec$bdsky.notInf.id,sep=''),
								xInclude=paste(tmp, collapse=' ',sep=''), 					
								spec= beast2.spec$bdsky.prior.spec), parent=bxml.prior, doc=bxml, addFinalizer=T)
				dummy			<- hivc.beast2.get.prior(beast2.spec$bdsky.notInf.prior[serial.id], serialwrapper, xmlAttrs(serialwrapper)["id"], bxml)				
			})	
	if(verbose)	cat(paste("\nadded serial becomingUninfectious priors, n=", length(dummy)))
	dummy		<- newXMLCommentNode(text="The serial becomeNoninfectiousAfterSamplingProbability priors.", parent=bxml.prior, doc=bxml, addFinalizer=T)
	dummy		<- lapply(seq_len(beast2.spec$bdsky.intervalNumber), function(serial.id)
			{
				tmp				<- rep(0,beast2.spec$bdsky.intervalNumber)
				tmp[serial.id]	<- 1
				serialwrapper	<- newXMLNode("distribution", attrs= list(	id=paste("sprior",serial.id,beast2.spec$sasky.r.id,sep='-'),
								name="distribution",
								x= paste('@',beast2.spec$sasky.r.id,sep=''),
								xInclude=paste(tmp, collapse=' ',sep=''), 					
								spec= beast2.spec$bdsky.prior.spec), parent=bxml.prior, doc=bxml, addFinalizer=T)
				dummy			<- hivc.beast2.get.prior(beast2.spec$sasky.r.prior[serial.id], serialwrapper, xmlAttrs(serialwrapper)["id"], bxml)				
			})	
	if(verbose)	cat(paste("\nadded serial becomeNoninfectiousAfterSamplingProbability priors, n=", length(dummy)))	
	dummy		<- newXMLCommentNode(text="end: serial SASKY priors.", parent=bxml.prior, doc=bxml, addFinalizer=T)
	bxml.prior
}
######################################################################################
hivc.beast2.get.specifications	<- function(xml.dir=NA, xml.filename=NA, mcmc.length=20e6, bdsky.intervalNumber=4)
{
	beast2.spec<- list()	
	beast2.spec$namespace			<- c("beast.core","beast.evolution.alignment","beast.evolution.tree.coalescent","beast.core.util","beast.evolution.nuc","beast.evolution.operators","beast.evolution.sitemodel","beast.evolution.substitutionmodel","beast.evolution.likelihood","beast.evolution.speciation","beast.core.parameter")
	beast2.spec$xml.dir				<- xml.dir		
	beast2.spec$xml.filename		<- xml.filename
	beast2.spec$map.Beta			<- "beast.math.distributions.Beta"
	beast2.spec$map.Exponential		<- "beast.math.distributions.Exponential"
	beast2.spec$map.ExcludablePrior	<- "beast.math.distributions.ExcludablePrior"
	beast2.spec$map.InverseGamma	<- "beast.math.distributions.InverseGamma"
	beast2.spec$map.LogNormal		<- "beast.math.distributions.LogNormalDistributionModel"
	beast2.spec$map.Gamma			<- "beast.math.distributions.Gamma"
	beast2.spec$map.Uniform			<- "beast.math.distributions.Uniform"
	beast2.spec$map.prior			<- "beast.math.distributions.Prior"
	beast2.spec$map.Laplace			<- "beast.math.distributions.LaplaceDistribution"
	beast2.spec$map.OneOnX			<- "beast.math.distributions.OneOnX"
	beast2.spec$map.Normal			<- "beast.math.distributions.Normal"
	beast2.spec$mcmc.length			<- mcmc.length
	beast2.spec$beast.label.datepos	<- 4
	beast2.spec$beast.label.sep		<- '_'
	beast2.spec$alignment.dataType	<- 'nucleotide'
	beast2.spec$alignment.id		<- 'ds'
	beast2.spec$alignment.missing	<- "-?"
	beast2.spec$tree.taxonset		<- beast2.spec$alignment.id
	beast2.spec$tree.id				<- paste('Tree',beast2.spec$tree.taxonset,sep='.t:')
	beast2.spec$sequence.totalcount	<- 4
	beast2.spec$datetrait.spec		<- "beast.evolution.tree.TraitSet"
	beast2.spec$datetrait.taxa.spec	<- "TaxonSet"
	beast2.spec$datetrait.units		<- "year"
	beast2.spec$datetrait.traitname	<- "date"
	beast2.spec$treemodel			<- "BirthDeathSkylineModel"
	beast2.spec$treemodel.id		<- paste("birthDeath",beast2.spec$tree.taxonset,sep='.t:')	
	beast2.spec$bdsky.spec			<- "beast.evolution.speciation.BirthDeathSkylineModel"
	beast2.spec$bdsky.prior.spec	<- "beast.math.distributions.ExcludablePrior"
	beast2.spec$bdsky.intervalNumber<- bdsky.intervalNumber
	beast2.spec$bdsky.origin.id		<- paste('originS',beast2.spec$tree.taxonset,sep='.t:')
	beast2.spec$bdsky.origin.value	<- rep(35, length(beast2.spec$bdsky.origin.id))
	beast2.spec$bdsky.origin.lower	<- 0.0
	beast2.spec$bdsky.origin.upper	<- 1000.0
	beast2.spec$bdsky.origin.prior	<- "Uniform/20.0/40.0"
	beast2.spec$bdsky.sprop.id					<- paste('samplingProportionS',beast2.spec$tree.taxonset,sep='.t:')
	beast2.spec$bdsky.sprop.value				<- rep(0.4, beast2.spec$bdsky.intervalNumber)
	beast2.spec$bdsky.sprop.lower				<- 0.0
	beast2.spec$bdsky.sprop.upper				<- 1.0
	beast2.spec$bdsky.sprop.prior				<- rep("Uniform/0.2/1.0",beast2.spec$bdsky.intervalNumber)
	beast2.spec$bdsky.sprop.changepoint.id		<- paste('samplingRateChangeTimes',beast2.spec$tree.taxonset,sep='.t:')
	beast2.spec$bdsky.sprop.changepoint.value	<- c(9.596, 5.596, 1.596, 0.)	
	beast2.spec$bdsky.R0.id						<- paste('R0S',beast2.spec$tree.taxonset,sep='.t:')
	beast2.spec$bdsky.R0.value					<- rep(1.2, beast2.spec$bdsky.intervalNumber)
	beast2.spec$bdsky.R0.lower					<- 0.0
	beast2.spec$bdsky.R0.upper					<- 10.0
	beast2.spec$bdsky.R0.prior					<- rep("Gamma/1.5/1.5/0",beast2.spec$bdsky.intervalNumber)
	beast2.spec$bdsky.R0.changepoint.id			<- paste('birthRateChangeTimes',beast2.spec$tree.taxonset,sep='.t:')
	beast2.spec$bdsky.R0.changepoint.value		<- c(9.596, 5.596, 1.596, 0.)	
	beast2.spec$bdsky.notInf.id					<- paste('becomeUninfectiousRateS',beast2.spec$tree.taxonset,sep='.t:')
	beast2.spec$bdsky.notInf.value				<- rep(0.1, beast2.spec$bdsky.intervalNumber)
	beast2.spec$bdsky.notInf.lower				<- 0.0
	beast2.spec$bdsky.notInf.upper				<- 10.0
	beast2.spec$bdsky.notInf.prior				<- rep("OneOnX/0",beast2.spec$bdsky.intervalNumber)
	beast2.spec$bdsky.notInf.changepoint.id		<- paste('deathRateChangeTimes',beast2.spec$tree.taxonset,sep='.t:')
	beast2.spec$bdsky.notInf.changepoint.value	<- c(9.596, 5.596, 1.596, 0.)
	beast2.spec$bdsky.reverseTimeArrays.id		<- paste('reverseTimeArrays',beast2.spec$tree.taxonset,sep='.t:')
	beast2.spec$bdsky.reverseTimeArrays.spec	<- "parameter.BooleanParameter"	
	beast2.spec$sasky.spec						<- "beast.evolution.speciation.SABDSkylineModel"
	beast2.spec$sasky.tree.spec					<- "beast.util.ClusterZBSATree"
	beast2.spec$sasky.tree.nodetype				<- "beast.evolution.tree.ZeroBranchSANode" 
	beast2.spec$sasky.tree.clusterType			<- "upgma"
	beast2.spec$sasky.r.id						<- paste('r',beast2.spec$tree.taxonset,sep='.t:')
	beast2.spec$sasky.r.value					<- rep(0.5, beast2.spec$bdsky.intervalNumber)
	beast2.spec$sasky.r.lower					<- 0.0
	beast2.spec$sasky.r.upper					<- 1.0
	beast2.spec$sasky.r.prior					<- rep("Uniform/0.0/1.0",beast2.spec$bdsky.intervalNumber)
	beast2.spec$sasky.r.changepoint.id			<- paste('rChangeTimes',beast2.spec$tree.taxonset,sep='.t:')
	beast2.spec$sasky.r.changepoint.value		<- c(9.596, 5.596, 1.596, 0.)		
	beast2.spec$compoundprior.spec	<- "util.CompoundDistribution"
	beast2.spec
}
######################################################################################
hivc.beast2.get.xml<- function(	bxml.template, seq.PROT.RT, df, beast2.spec, ph=NULL, verbose=1)
{	
	require(XML)
	#	init XML
	bxml		<- hivc.beast2.init.xml( beast2.spec, verbose=verbose)
	bxml.beast	<- getNodeSet(bxml, "//beast")[[1]]	
	#	add alignment	
	dummy		<- hivc.beast2.add.alignment(bxml, seq.PROT.RT, df, beast2.spec, verbose=verbose)
	#	add tip dates
	dummy		<- hivc.beast2.add.datetrait(bxml, df, beast2.spec, verbose=verbose)
	#	add tree for alignment
	if(beast2.spec$treemodel=="BirthDeathSkylineModel")
		dummy	<- hivc.beast2.add.tree(bxml, beast2.spec, verbose=verbose)
	if(beast2.spec$treemodel=="SampledAncestorSkylineModel")
		dummy	<- hivc.beast2.add.satree(bxml, beast2.spec, verbose=verbose)	
	#	add tree model
	if(beast2.spec$treemodel=="BirthDeathSkylineModel")
		dummy	<- hivc.beast2.add.treemodel.bdsky(bxml, beast2.spec, verbose=verbose)
	if(beast2.spec$treemodel=="SampledAncestorSkylineModel")
		dummy	<- hivc.beast2.add.treemodel.sasky(bxml, beast2.spec, verbose=verbose)	
	#	copy branchRateModel from template
	tmp			<- getNodeSet(bxml.template, "//branchRateModel")
	if(length(tmp)!=1) stop("unexpected number of //branchRateModel")
	dummy<- addChildren( bxml.beast, xmlClone( tmp[[1]], addFinalizer=T, doc=bxml ) )
	if(verbose)	cat(paste("\nadded branchRateModel from template, size=", xmlSize(tmp[[1]])))	
	#	copy siteModel from template
	tmp			<- getNodeSet(bxml.template, "//siteModel")
	if(length(tmp)!=1) stop("unexpected number of //siteModel")
	dummy<- addChildren( bxml.beast, xmlClone( tmp[[1]], addFinalizer=T, doc=bxml ) )
	if(verbose)	cat(paste("\nadded siteModel from template, size=", xmlSize(tmp[[1]])))
	#	copy run from template
	tmp			<- getNodeSet(bxml.template, "//run")
	if(length(tmp)!=1) stop("unexpected number of //run")
	dummy		<- addChildren( bxml.beast, xmlClone( tmp[[1]], addFinalizer=T, doc=bxml ) )
	if(verbose)	cat(paste("\nadded run from template, size=", xmlSize(tmp[[1]])))
	#	add tree model prior
	if(beast2.spec$treemodel=="BirthDeathSkylineModel")
		dummy	<- hivc.beast2.add.bdsky.serialpriors(bxml, beast2.spec, verbose=verbose)
	if(beast2.spec$treemodel=="SampledAncestorSkylineModel")
		dummy	<- hivc.beast2.add.sasky.serialpriors(bxml, beast2.spec, verbose=verbose)
	
	#	reset output fileNames
	bxml.onodes	<- getNodeSet(bxml, "//*[@fileName]")
	tmp			<- sapply(bxml.onodes, function(x) xmlGetAttr(x,"fileName"))	
	tmp			<- sapply(strsplit(tmp,'.',fixed=1), function(x)	paste(beast2.spec$xml.filename, '.', rev(x)[1], sep=''))		
	dummy		<- sapply(seq_along(bxml.onodes), function(i){		xmlAttrs(bxml.onodes[[i]])["fileName"]<- tmp[i]		})
	if(verbose)	cat(paste("\nchanged trunk filename to", beast2.spec$xml.filename))
	#	reset chain length and logEvery
	bxml.onodes	<- getNodeSet(bxml, "//*[@chainLength]")
	dummy		<- sapply(seq_along(bxml.onodes), function(i){		xmlAttrs(bxml.onodes[[i]])["chainLength"]<- sprintf("%d",beast2.spec$mcmc.length)		})
	if(verbose)	cat(paste("\nchanged chain length to", beast2.spec$mcmc.length))
	bxml.onodes	<- getNodeSet(bxml, "//*[@logEvery]")
	dummy		<- sapply(seq_along(bxml.onodes), function(i){		xmlAttrs(bxml.onodes[[i]])["logEvery"]<- sprintf("%d",round(beast2.spec$mcmc.length/1e4))		})
	if(verbose)	cat(paste("\nchanged logEvery to", round(beast2.spec$mcmc.length/1e4)))
	
	bxml	
}
######################################################################################
#	create xml file from btemplate and seq.PROT.RT, using seq in df 
# 	beast.label.datepos= 4; beast.label.sep= '_'; beast.date.direction= "forwards"; beast.date.units= "years"; verbose=1; xml.prior4tipstem="uniform"; xml.resetTipDate2LastDiag=1
hivc.beast.get.xml<- function(	btemplate, seq.PROT.RT, df, file, ph=NULL, xml.monophyly4clusters=0, xml.taxon4tipstem=0, xml.prior4tipstem=NA, 
								beast.label.datepos= 4, beast.label.sep= '_', beast.date.direction= "forwards", beast.usingDates="true", beast.date.units= "years", beast.mcmc.chainLength=50000000, verbose=1)
{
	
	bxml		<- newXMLDoc(addFinalizer=T)
	bxml.beast	<- newXMLNode("beast", doc=bxml, addFinalizer=T)
	newXMLCommentNode(text=paste("Generated by HIVCLUST from template",file), parent=bxml.beast, doc=bxml, addFinalizer=T)
	#	add new set of sequences
	dummy		<- hivc.beast.add.seq(bxml, df, seq.PROT.RT, beast.label.datepos=beast.label.datepos, beast.label.sep=beast.label.sep, beast.date.direction=beast.date.direction, beast.date.units=beast.date.units, verbose=verbose)
	#	copy everything after alignment up to but not including constantSize from template
	bt.beast	<- getNodeSet(btemplate, "//beast")[[1]]
	dummy		<- sapply(seq.int( which( xmlSApply(bt.beast, xmlName)=="alignment" )+1, which( xmlSApply(bt.beast, xmlName)=="patterns" ) ), function(i)
			{
				if( class(bt.beast[[i]])[1]=="XMLInternalCommentNode" )
					dummy<- newXMLCommentNode(text=xmlValue(bt.beast[[i]]), parent=bxml.beast, doc=bxml, addFinalizer=T)
				else
					dummy<- addChildren( bxml.beast, xmlClone( bt.beast[[i]], addFinalizer=T, doc=bxml ) )
			})
	#	add startingTree
	if(!is.null(ph))		#	add startingTree if provided
		dummy	<- hivc.beast.add.startingtree(bxml, ph, df, beast.usingDates=beast.usingDates)		
	else					#	otherwise copy upgmaTree
	{
		dummy		<- sapply(seq.int( which( xmlSApply(bt.beast, xmlName)=="constantSize" )-2, which( xmlSApply(bt.beast, xmlName)=="upgmaTree" ) ), function(i)
				{
					if( class(bt.beast[[i]])[1]=="XMLInternalCommentNode" )
						dummy<- newXMLCommentNode(text=xmlValue(bt.beast[[i]]), parent=bxml.beast, doc=bxml, addFinalizer=T)
					else
						dummy<- addChildren( bxml.beast, xmlClone( bt.beast[[i]], addFinalizer=T, doc=bxml ) )
				})
	}
	#	copy everything from 'treeModel'-1 from template
	dummy		<- sapply(seq.int( which( xmlSApply(bt.beast, xmlName)=="treeModel" )-1, xmlSize(bt.beast) ), function(i)
			{
				if( class(bt.beast[[i]])[1]=="XMLInternalCommentNode" )
					dummy<- newXMLCommentNode(text=xmlValue(bt.beast[[i]]), parent=bxml.beast, doc=bxml, addFinalizer=T)
				else
					dummy<- addChildren( bxml.beast, xmlClone( bt.beast[[i]], addFinalizer=T, doc=bxml ) )
			})
	#	if user-provided startingTree, reset <upgmaTree idref="startingTree"/> to <newick idref="startingTree"/>
	if(!is.null(ph))
	{
		tmp					<- getNodeSet(bxml, "//*[@idref='startingTree']")
		if(length(tmp)!=1) stop("unexpected number of //*[@idref='startingTree']")
		xmlName(tmp[[1]])	<- "newick"
	}
	# 	reset dimension of GMRF likelihood
	tmp			<- getNodeSet(bxml, "//*[@id='skyride.logPopSize']")
	if(length(tmp)!=1)	stop("unexpected number of *[@id='skyride.logPopSize'")
	tmp			<- tmp[[1]]
	xmlAttrs(tmp)["dimension"]	<-	nrow(df)-1  
	tmp			<- getNodeSet(bxml, "//*[@id='skyride.groupSize']")
	if(length(tmp)!=1)	stop("unexpected number of *[@id='skyride.groupSize'")
	tmp			<- tmp[[1]]
	xmlAttrs(tmp)["dimension"]	<-	nrow(df)-1
	#	if uniform prior for rootheight, set minimum to earliest sample time in data set
	dummy		<- hivc.beast.adjust.rootheightprior(bxml, df, verbose=verbose)
	dummy		<- hivc.beast.adjust.mcmc(bxml, beast.mcmc.chainLength=beast.mcmc.chainLength, verbose=verbose)
	#	for tips, add taxon sets and tmrcaStatistics
	if(xml.taxon4tipstem)
		dummy	<- hivc.beast.add.taxonsets4tips(bxml, df, log=1, verbose=1)
	if(!is.na(xml.prior4tipstem))
		dummy	<- hivc.beast.add.prior4tips(bxml, df, xml.prior4stem=xml.prior4tipstem, beast.label.datepos=beast.label.datepos, verbose=verbose)

	#	for clusters, add taxon sets and tmrcaStatistics 
	dummy		<- hivc.beast.add.taxonsets4clusters(bxml, df, xml.monophyly4clusters=xml.monophyly4clusters, verbose=verbose)		
	#	reset output fileNames
	bxml.onodes	<- getNodeSet(bxml, "//*[@fileName]")
	tmp			<- sapply(bxml.onodes, function(x) xmlGetAttr(x,"fileName"))
	tmp			<- gsub("(time).","time",tmp,fixed=1)
	tmp			<- gsub("(subst).","subst",tmp,fixed=1)
	tmp			<- sapply(strsplit(tmp,'.',fixed=1), function(x)	paste(file, '.', x[2], sep=''))		
	dummy		<- sapply(seq_along(bxml.onodes), function(i){		xmlAttrs(bxml.onodes[[i]])["fileName"]<- tmp[i]		})
	#
	bxml
}	
######################################################################################
#	adjust mcmc BEAST XML element
hivc.beast.adjust.mcmc<- function(bxml, beast.mcmc.chainLength=50000000, verbose=1)
{
	tmp			<- getNodeSet(bxml, "//*[@id='mcmc']")
	if(length(tmp)!=1)	stop("unexpected number of *[@id='mcmc']")
	tmp			<- tmp[[1]]
	if(verbose)	cat(paste("\nSet MCMC chainLength to",beast.mcmc.chainLength))
	xmlAttrs(tmp)["chainLength"]	<-	sprintf("%d",beast.mcmc.chainLength)
	bxml
}
######################################################################################
#	if rootheight prior uniform, sets lower bound to earliest sampling time in data set
hivc.beast.adjust.rootheightprior<- function(bxml, df, verbose=1)
{
	bxml.prior	<- getNodeSet(bxml, "//uniformPrior[descendant::parameter[@idref='treeModel.rootHeight']]")
	if(length(bxml.prior)>1)	stop("unexpected length of //uniformPrior[descendant::parameter[@idref='treeModel.rootHeight']]")
	if(length(bxml.prior)==1)
	{
		tmp								<- df[,range(AnyPos_T1)]
		if(verbose)	cat(paste("\nfound uniformPrior for treeModel.rootHeight. Range of tip dates is",tmp[1],tmp[2]))
		tmp								<- difftime(tmp[2],tmp[1], units="days")	
		bxml.prior						<- bxml.prior[[1]]
		if(verbose)	cat(paste("\nset lower bound of uniformPrior for treeModel.rootHeight to", floor( tmp/365 )))
		xmlAttrs(bxml.prior)["lower"]	<- floor( tmp/365 )	
	}
	bxml
}
######################################################################################
#	For each cluster, create a taxonset. Assumes df has BEASTlabel and cluster
hivc.beast.get.taxonsets4clusters	<- function(bxml, df)
{	
	ans	<- lapply( unique( df[,cluster] ), function(clu)
			{							
				taxonset	<- newXMLNode("taxa", attrs= list(id=paste("c",clu,sep='')), doc=bxml, addFinalizer=T )
				tmp<- df[cluster==clu,][,BEASTlabel]
				tmp<- lapply(tmp, function(x)		newXMLNode("taxon", attrs= list(idref=x), parent=taxonset, doc=bxml, addFinalizer=T )	)
				taxonset
			})
	ans		
}
######################################################################################
#	For each tip, create a taxonset. Assumes df has BEASTlabel
hivc.beast.get.taxonsets4tips	<- function(bxml, df)
{	
	ans	<- lapply( seq_len(nrow(df)), function(i)
			{							
				taxonset	<- newXMLNode("taxa", attrs= list(id=paste("tip",i,sep='')), doc=bxml, addFinalizer=T )
				newXMLNode("taxon", attrs= list(idref=df[i,BEASTlabel]), parent=taxonset, doc=bxml, addFinalizer=T )
				taxonset
			})
	ans		
}
######################################################################################
# 	extract starting tree from 'ph' by tips in 'df'. Only keeps tree topology and resets branch lengths so that the maximum root distance is 'beast.rootHeight'
#	beast.rootHeight= 35; beast.usingDates= "false"; beast.newickid= "startingTree"
hivc.beast.add.startingtree<- function(bxml, ph, df, beast.rootHeight= 35, beast.usingDates="true", beast.newickid= "startingTree", beast.brlunits="years", verbose=1)
{
	require(adephylo)
	if(verbose) cat(paste("\ncreate startingTree with root height=",beast.rootHeight))
	tmp					<- setdiff( ph$tip.label, df[,FASTASampleCode] )
	tmp					<- match( tmp, ph$tip.label)
	ph.start			<- drop.tip(ph, tmp)		
	ph.start$node.label	<- NULL
	setkey(df, FASTASampleCode)
	ph.start$tip.label	<- df[ph.start$tip.label,][,BEASTlabel]
	if(verbose) cat(paste("\nselected tips for startingTree, n=",Ntip(ph.start)))
	#	adjust rootHeight to 'beast.rootHeight'
	tmp					<- beast.rootHeight / max(distRoot(ph.start))
	ph.start$edge.length<- ph.start$edge.length*tmp
	#	compute adjusted branch lengths for each tip: midpoint within NegT and AnyPos_T1
	df.length			<- t( sapply( strsplit(df[,BEASTlabel],'_',fixed=1), function(x)  as.numeric( x[2:4]) ) )
	df.length			<- data.table(BEASTlabel=df[,BEASTlabel], NegT=df.length[,1], AnyPos_T1=df.length[,2], PosSeqT=df.length[,3])
	tmp					<- max( df.length[, PosSeqT])		#TODO should this be height or length ?
	df.length			<- df.length[, list(BEASTlabel=BEASTlabel, brl=(AnyPos_T1-NegT)/2+PosSeqT-AnyPos_T1)]
	setkey(df.length, BEASTlabel)
	#	adjust stem of each tip to be within NegT and AnyPos_T1
	tmp							<- sapply(seq_len(Ntip(ph.start)), function(x) which( ph.start$edge[,2]==x ) )
	ph.start$edge.length[ tmp ]	<- df.length[ph.start$tip.label,][,brl]
	if(verbose) cat(paste("\nadjusted branch lengths of tips to be within NegT and AnyPos_T1. New root height is",max(distRoot(ph.start))))
	#	write ph.start as newick tree to bxml
	tmp					<- write.tree( ph.start )	
	bxml.beast			<- getNodeSet(bxml, "//beast")[[1]]
	dummy				<- newXMLCommentNode(text="The user-specified starting tree in a newick tree format", parent=bxml.beast, doc=bxml, addFinalizer=T)
	bxml.startingTree	<- newXMLNode("newick", attrs= list(id=beast.newickid, usingDates=beast.usingDates, units=beast.brlunits), parent= bxml.beast, doc=bxml, addFinalizer=T)
	dummy				<- newXMLTextNode(text=tmp, parent=bxml.startingTree, doc=bxml, addFinalizer=T) 
	bxml
}
######################################################################################
#	get a list of monopylyStatistics. assumes each taxonset in 'btaxonsets' is monophyletic
hivc.beast.get.monophylyStatistic<- function(bxml, btaxonsets, treeModel.id) 
{		
	btaxonset.id	<- sapply(btaxonsets, function(x)	xmlGetAttr(x,"id")	)
	ans				<- lapply(btaxonset.id, function(x)
			{
				monophylyStatistic	<- newXMLNode("monophylyStatistic", attrs= list(id=paste("monophyly(",x,")",sep='')), doc=bxml, addFinalizer=T )
				mrca			<- newXMLNode("mrca", parent=monophylyStatistic, doc=bxml, addFinalizer=T)
				newXMLNode("taxa", attrs= list(idref=x), parent=mrca, doc=bxml, addFinalizer=T )
				newXMLNode("treeModel", attrs= list(idref=treeModel.id), parent=monophylyStatistic, doc=bxml, addFinalizer=T )
				monophylyStatistic					
			})
	ans
}	
######################################################################################
#	add a list of monopylyStatistics to the BEAST prior. assumes all monophylyStatistics are referenced in the list 'monophylyStatistics'
hivc.beast.add.monophylylkl<- function(bxml, monophylyStatistics)
{
	bxml.prior				<- getNodeSet(bxml, "//*[@id='prior']")
	if(length(bxml.prior)!=1)	stop("unexpected length of bxml.prior")
	bxml.prior				<- bxml.prior[[1]]
	#see if there is a booleanLikelihood, and if yes use it, otherwise add a new XML node
	bxml.bool.lkl			<- getNodeSet(bxml.prior,"//booleanLikelihood")
	if(length(bxml.bool.lkl)>1)	stop("unexpected length of bxml.bool.lkl")
	if(length(bxml.bool.lkl)<1)
		bxml.bool.lkl		<- newXMLNode("booleanLikelihood", parent=bxml.prior, doc=bxml, addFinalizer=T )
	else
		bxml.bool.lkl		<- bxml.bool.lkl[[1]]
	
	monophylyStatistics.id	<- sapply(monophylyStatistics, function(x)	xmlGetAttr(x,"id")	)
	dummy					<- lapply(monophylyStatistics.id, function(x)
			{
				dummy		<- newXMLNode("monophylyStatistic", attrs= list(idref=x), parent=bxml.bool.lkl, doc=bxml, addFinalizer=T )
			})
	bxml
}
######################################################################################
#	For each taxonset, get a list of tmrcaStatistics. Assumes all tmrcaStatistics share a treeModel.id and the same includeStem attribute
hivc.beast.get.tmrcaStatistic<- function(bxml, btaxonsets, treeModel.id, includeStem="false") 
{		
	btaxonset.id	<- sapply(btaxonsets, function(x)	xmlGetAttr(x,"id")	)
	prefix.id		<- ifelse(includeStem=="false","tmrca","tstem")
	ans				<- lapply(btaxonset.id, function(x)
			{
				tmrcaStatistic	<- newXMLNode("tmrcaStatistic", attrs= list(id=paste(prefix.id,"(",x,")",sep=''), includeStem=includeStem), doc=bxml, addFinalizer=T )
				mrca			<- newXMLNode("mrca", parent=tmrcaStatistic, doc=bxml, addFinalizer=T)
				newXMLNode("taxa", attrs= list(idref=x), parent=mrca, doc=bxml, addFinalizer=T )
				newXMLNode("treeModel", attrs= list(idref=treeModel.id), parent=tmrcaStatistic, doc=bxml, addFinalizer=T )
				tmrcaStatistic					
			})
	ans
}	
######################################################################################
#	For each tip: construct a prior for the corresponding tmrcaStatistics
hivc.beast.get.tipPrior<- function(bxml, df, btmrcaStatistics.tips, xml.prior4stem="uniform", beast.label.negpos=2, beast.label.diagpos=3, beast.label.datepos=4, verbose=1)
{
	if(xml.prior4stem!="uniform")	stop("unexpected xml.tipprior")
	#
	if(verbose)	cat(paste("\nuse tip date found at pos x in label, x=", beast.label.datepos))
	df.height	<- t( sapply( strsplit(df[,BEASTlabel],'_',fixed=1), function(x)  as.numeric( x[c(beast.label.negpos, beast.label.diagpos, beast.label.datepos)]) ) )
	df.height	<- data.table(BEASTlabel=df[,BEASTlabel], NegT=df.height[,1], AnyPos_T1=df.height[,2], PosSeqT=df.height[,3])
	tmp			<- max( df.height[, PosSeqT])		#TODO should this be height or length ?
	df.height	<- df.height[, list(BEASTlabel=BEASTlabel, NegT=tmp-NegT, AnyPos_T1=tmp-AnyPos_T1, PosSeqT=tmp-PosSeqT)]
	#add uniform prior according to last NegT and first PosT
	ans			<- lapply(btmrcaStatistics.tips, function(x)
					{						
						tmrcaStatistics.id	<- xmlGetAttr(x,"id")
						tmp					<- xpathApply(x, "mrca/taxa", xmlGetAttr, "idref" )
						if(length(tmp)!=1)	stop("unexpected length of idref for mrca/taxa") 
						tip					<- as.numeric( substr(tmp[[1]],4,nchar(tmp[[1]])) )												
						bxml.tipprior		<- newXMLNode("uniformPrior", attrs= list(lower=df.height[tip,AnyPos_T1], upper=df.height[tip,NegT]), doc=bxml, addFinalizer=T )
						dummy				<- newXMLNode("statistic", attrs= list(idref=tmrcaStatistics.id), parent=bxml.tipprior, doc=bxml, addFinalizer=T )
						bxml.tipprior
					})
	ans			
}
######################################################################################
#	For each tip: add a taxonset, tmrcaStatistic, and potentially a reference to fileLog to 'bxml'
#	The aim of this as standalone is to log the length of tip stems. 
hivc.beast.add.taxonsets4tips<- function(bxml, df, log=1, verbose=1)
{
	bxml.treeModel.id			<- unlist(xpathApply(bxml, "//treeModel[@id]", xmlGetAttr, "id"))
	
	#get taxon sets for each tip
	btaxonsets.tips				<- hivc.beast.get.taxonsets4tips(bxml, df)	
	#get tmrcaStatistic for each tip
	btmrcaStatistics.tips		<- hivc.beast.get.tmrcaStatistic(bxml, btaxonsets.tips, bxml.treeModel.id, includeStem="true") 			
	#	modify from template
	bxml.beast					<- getNodeSet(bxml, "//beast")[[1]]
	#	add 'btaxonsets.tips' after last taxa
	bxml.idx					<- which(xmlSApply(bxml.beast, xmlName)=="taxa")		
	addChildren(bxml.beast, btaxonsets.tips, at=bxml.idx[length(bxml.idx)] )
	if(verbose) cat(paste("\nadded taxon sets for each tip, n=",length(btaxonsets.tips)))
	#	add 'btmrcaStatistics.tips' after last treeModel or tmrcaStatistic
	bxml.idx					<- which(xmlSApply(bxml.beast, xmlName)%in%c("treeModel","tmrcaStatistic"))
	addChildren(bxml.beast, btmrcaStatistics.tips, at=tail(bxml.idx,1) )
	if(verbose) cat(paste("\nadded tmrcaStatistics for stem of each tip, n=",length(btmrcaStatistics.tips)))
	#add tmrcaStatistics to log 
	if(log)
	{
		bxml.fileLog				<- getNodeSet(bxml, "//log[@id='fileLog']")
		tmrcaStatistics.id			<- sapply(btmrcaStatistics.tips, function(x)	xmlGetAttr(x,"id")	)
		tmp							<- lapply(tmrcaStatistics.id, function(x)	newXMLNode("tmrcaStatistic", attrs= list(idref=x), parent=bxml.fileLog, doc=bxml, addFinalizer=T )	)
		if(verbose) cat(paste("\nadded tmrcaStatistics for tip stems to fileLog, n=",length(tmrcaStatistics.id)))
	}
	
	bxml
}
######################################################################################
#	For each tip: add a taxonset, tmrcaStatistic, prior for the tmrcaStatistics and potentially a reference to fileLog to 'bxml'
hivc.beast.add.prior4tips<- function(bxml, df, xml.prior4stem="uniform", beast.label.datepos=4, verbose=1)
{
	#find list of tmrcaStatistics with id containing 'tip'
	btmrcaStatistics.tips		<- getNodeSet(bxml, "//tmrcaStatistic[starts-with(@id,'tstem(tip')]")
	print(btmrcaStatistics.tips)
	#get prior for each tip stem
	bprior.tips					<- hivc.beast.get.tipPrior(bxml, df, btmrcaStatistics.tips, xml.prior4stem=xml.prior4stem, beast.label.datepos=beast.label.datepos, verbose=verbose)		
	#	add 'bprior.tips' to prior
	bxml.prior					<- getNodeSet(bxml, "//prior[@id='prior']")
	if(length(bxml.prior)!=1)	stop("unexpected number of //prior[@id='prior']")
	bxml.prior					<- bxml.prior[[1]]
	addChildren(bxml.prior, bprior.tips)
	if(verbose) cat(paste("\nadded priors for stem of each tip, n=",length(bprior.tips)))	
	bxml
}
######################################################################################
hivc.beast.add.taxonsets4clusters<- function(bxml, df, xml.monophyly4clusters=1, verbose=1)
{
	bxml.treeModel.id			<- unlist(xpathApply(bxml, "//treeModel[@id]", xmlGetAttr, "id"))
	
	#get taxon sets for each cluster
	btaxonsets.clusters			<- hivc.beast.get.taxonsets4clusters(bxml, df)	
	#get tmrcaStatistic for each cluster
	btmrcaStatistics.clusters	<- hivc.beast.get.tmrcaStatistic(bxml, btaxonsets.clusters, bxml.treeModel.id, includeStem="false") 		
	
	#	modify from template
	bxml.beast					<- getNodeSet(bxml, "//beast")[[1]]
	#	add 'btaxonsets.clusters' after last taxa
	bxml.idx					<- which(xmlSApply(bxml.beast, xmlName)=="taxa")		
	addChildren(bxml.beast, btaxonsets.clusters, at=bxml.idx[length(bxml.idx)] )
	if(verbose) cat(paste("\nadded taxon sets comprising each cluster, n=",length(btaxonsets.clusters)))
	#	add 'btmrcaStatistics.clusters' after last treeModel
	bxml.idx					<- which(xmlSApply(bxml.beast, xmlName)=="treeModel")
	addChildren(bxml.beast, btmrcaStatistics.clusters, at=bxml.idx[length(bxml.idx)] )
	if(verbose) cat(paste("\nadded tmrcaStatistics for mrca of each cluster, n=",length(btmrcaStatistics.clusters)))
	#add tmrcaStatistics to log 
	bxml.fileLog				<- getNodeSet(bxml, "//log[@id='fileLog']")
	tmrcaStatistics.id			<- sapply(btmrcaStatistics.clusters, function(x)	xmlGetAttr(x,"id")	)
	tmp							<- lapply(tmrcaStatistics.id, function(x)	newXMLNode("tmrcaStatistic", attrs= list(idref=x), parent=bxml.fileLog, doc=bxml, addFinalizer=T )	)
	if(verbose) cat(paste("\nadded tmrcaStatistics to fileLog"))
	#
	#	get monophylyStatistic and add after last 'tmrcaStatistic'
	#	populate booleanLikelihood with monophylyStatistic constraints
	#
	if(xml.monophyly4clusters)
	{
		bmStatistics.clusters	<- hivc.beast.get.monophylyStatistic(bxml, btaxonsets.clusters, bxml.treeModel.id)
		bxml.idx				<- which(xmlSApply(bxml.beast, xmlName)=="tmrcaStatistic")
		addChildren(bxml.beast, bmStatistics.clusters, at=bxml.idx[length(bxml.idx)] )
		if(verbose) cat(paste("\nadded monophylyStatistics for mrca of each cluster, n=",length(bmStatistics.clusters)))
		hivc.beast.add.monophylylkl(bxml, bmStatistics.clusters)		
		if(verbose) cat(paste("\nadded monophylyStatistics to booleanLikelihood, n=",length(bmStatistics.clusters)))
	}
	
	bxml
}
######################################################################################
hivc.phy.plotupon<- function (x, type = "phylogram", use.edge.length = TRUE, node.pos = NULL, 
		show.tip.label = TRUE, show.node.label = FALSE, edge.color = "black", 
		edge.width = 1, edge.lty = 1, font = 3, cex = par("cex"), 
		adj = NULL, srt = 0, no.margin = FALSE, root.edge = FALSE, 
		label.offset = 0, underscore = FALSE, y.lim = NULL, lab4ut = "horizontal", 
		tip.color = "black", ...) 
{
	lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
	direction <- lastPP$direction
	type <- lastPP$type
	x.lim <- lastPP$x.lim
	Ntip <- length(x$tip.label)
	if (Ntip == 1) {
		warning("found only one tip in the tree")
		return(NULL)
	}
	if (any(tabulate(x$edge[, 1]) == 1)) 
		stop("there are single (non-splitting) nodes in your tree; you may need to use collapse.singles()")
	.nodeHeight <- function(Ntip, Nnode, edge, Nedge, yy) .C("node_height", 
				as.integer(Ntip), as.integer(Nnode), as.integer(edge[, 
								1]), as.integer(edge[, 2]), as.integer(Nedge), as.double(yy), 
				DUP = FALSE, PACKAGE = "ape")[[6]]
	.nodeDepth <- function(Ntip, Nnode, edge, Nedge) .C("node_depth", 
				as.integer(Ntip), as.integer(Nnode), as.integer(edge[, 
								1]), as.integer(edge[, 2]), as.integer(Nedge), double(Ntip + 
								Nnode), DUP = FALSE, PACKAGE = "ape")[[6]]
	.nodeDepthEdgelength <- function(Ntip, Nnode, edge, Nedge, 
			edge.length) .C("node_depth_edgelength", as.integer(Ntip), 
				as.integer(Nnode), as.integer(edge[, 1]), as.integer(edge[, 
								2]), as.integer(Nedge), as.double(edge.length), double(Ntip + 
								Nnode), DUP = FALSE, PACKAGE = "ape")[[7]]
	Nedge <- dim(x$edge)[1]
	Nnode <- x$Nnode
	ROOT <- Ntip + 1
	type <- match.arg(type, c("phylogram", "cladogram", "fan", 
					"unrooted", "radial"))
	direction <- match.arg(direction, c("rightwards", "leftwards", 
					"upwards", "downwards"))
	if (is.null(x$edge.length)) 
		use.edge.length <- FALSE
	if (type %in% c("unrooted", "radial") || !use.edge.length || 
			is.null(x$root.edge) || !x$root.edge) 
		root.edge <- FALSE
	if (type == "fan" && root.edge) {
		warning("drawing root edge with type = 'fan' is not yet supported")
		root.edge <- FALSE
	}
	phyloORclado <- type %in% c("phylogram", "cladogram")
	horizontal <- direction %in% c("rightwards", "leftwards")
	xe <- x$edge
	if (phyloORclado) {
		phyOrder <- attr(x, "order")
		if (is.null(phyOrder) || phyOrder != "cladewise") {
			x <- reorder(x)
			if (!identical(x$edge, xe)) {
				ereorder <- match(x$edge[, 2], xe[, 2])
				if (length(edge.color) > 1) {
					edge.color <- rep(edge.color, length.out = Nedge)
					edge.color <- edge.color[ereorder]
				}
				if (length(edge.width) > 1) {
					edge.width <- rep(edge.width, length.out = Nedge)
					edge.width <- edge.width[ereorder]
				}
				if (length(edge.lty) > 1) {
					edge.lty <- rep(edge.lty, length.out = Nedge)
					edge.lty <- edge.lty[ereorder]
				}
			}
		}
		yy <- numeric(Ntip + Nnode)
		TIPS <- x$edge[x$edge[, 2] <= Ntip, 2]
		yy[TIPS] <- 1:Ntip
	}
	z <- reorder(x, order = "pruningwise")
	if (phyloORclado) {
		if (is.null(node.pos)) {
			node.pos <- 1
			if (type == "cladogram" && !use.edge.length) 
				node.pos <- 2
		}
		if (node.pos == 1) 
			yy <- .nodeHeight(Ntip, Nnode, z$edge, Nedge, yy)
		else {
			ans <- .C("node_height_clado", as.integer(Ntip), 
					as.integer(Nnode), as.integer(z$edge[, 1]), as.integer(z$edge[, 
									2]), as.integer(Nedge), double(Ntip + Nnode), 
					as.double(yy), DUP = FALSE, PACKAGE = "ape")
			xx <- ans[[6]] - 1
			yy <- ans[[7]]
		}
		if (!use.edge.length) {
			if (node.pos != 2) 
				xx <- .nodeDepth(Ntip, Nnode, z$edge, Nedge) - 
						1
			xx <- max(xx) - xx
		}
		else {
			xx <- .nodeDepthEdgelength(Ntip, Nnode, z$edge, Nedge, 
					z$edge.length)
		}
	}
	else switch(type, fan = {
					TIPS <- x$edge[which(x$edge[, 2] <= Ntip), 2]
					xx <- seq(0, 2 * pi * (1 - 1/Ntip), 2 * pi/Ntip)
					theta <- double(Ntip)
					theta[TIPS] <- xx
					theta <- c(theta, numeric(Nnode))
					theta <- .nodeHeight(Ntip, Nnode, z$edge, Nedge, theta)
					if (use.edge.length) {
						r <- .nodeDepthEdgelength(Ntip, Nnode, z$edge, Nedge, 
								z$edge.length)
					} else {
						r <- .nodeDepth(Ntip, Nnode, z$edge, Nedge)
						r <- 1/r
					}
					xx <- r * cos(theta)
					yy <- r * sin(theta)
				}, unrooted = {
					nb.sp <- .nodeDepth(Ntip, Nnode, z$edge, Nedge)
					XY <- if (use.edge.length) unrooted.xy(Ntip, Nnode, z$edge, 
										z$edge.length, nb.sp) else unrooted.xy(Ntip, Nnode, 
										z$edge, rep(1, Nedge), nb.sp)
					xx <- XY$M[, 1] - min(XY$M[, 1])
					yy <- XY$M[, 2] - min(XY$M[, 2])
				}, radial = {
					X <- .nodeDepth(Ntip, Nnode, z$edge, Nedge)
					X[X == 1] <- 0
					X <- 1 - X/Ntip
					yy <- c((1:Ntip) * 2 * pi/Ntip, rep(0, Nnode))
					Y <- .nodeHeight(Ntip, Nnode, z$edge, Nedge, yy)
					xx <- X * cos(Y)
					yy <- X * sin(Y)
				})
	if (phyloORclado) {
		if (!horizontal) {
			tmp <- yy
			yy <- xx
			xx <- tmp - min(tmp) + 1
		}
		if (root.edge) {
			if (direction == "rightwards") 
				xx <- xx + x$root.edge
			if (direction == "upwards") 
				yy <- yy + x$root.edge
		}
	}
	if (no.margin) 
		par(mai = rep(0, 4))
	if (is.null(x.lim)) {
		if (phyloORclado) {
			if (horizontal) {
				x.lim <- c(0, NA)
				pin1 <- par("pin")[1]
				strWi <- strwidth(x$tip.label, "inches")
				xx.tips <- xx[1:Ntip] * 1.04
				alp <- try(uniroot(function(a) max(a * xx.tips + 
													strWi) - pin1, c(0, 1e+06))$root, silent = TRUE)
				if (is.character(alp)) 
					tmp <- max(xx.tips) * 1.5
				else {
					tmp <- if (show.tip.label) 
								max(xx.tips + strWi/alp)
							else max(xx.tips)
				}
				x.lim[2] <- tmp
			}
			else x.lim <- c(1, Ntip)
		}
		else switch(type, fan = {
						if (show.tip.label) {
							offset <- max(nchar(x$tip.label) * 0.018 * max(yy) * 
											cex)
							x.lim <- c(min(xx) - offset, max(xx) + offset)
						} else x.lim <- c(min(xx), max(xx))
					}, unrooted = {
						if (show.tip.label) {
							offset <- max(nchar(x$tip.label) * 0.018 * max(yy) * 
											cex)
							x.lim <- c(0 - offset, max(xx) + offset)
						} else x.lim <- c(0, max(xx))
					}, radial = {
						if (show.tip.label) {
							offset <- max(nchar(x$tip.label) * 0.03 * cex)
							x.lim <- c(-1 - offset, 1 + offset)
						} else x.lim <- c(-1, 1)
					})
	}
	else if (length(x.lim) == 1) {
		x.lim <- c(0, x.lim)
		if (phyloORclado && !horizontal) 
			x.lim[1] <- 1
		if (type %in% c("fan", "unrooted") && show.tip.label) 
			x.lim[1] <- -max(nchar(x$tip.label) * 0.018 * max(yy) * 
							cex)
		if (type == "radial") 
			x.lim[1] <- if (show.tip.label) 
						-1 - max(nchar(x$tip.label) * 0.03 * cex)
					else -1
	}
	if (phyloORclado && direction == "leftwards") 
		xx <- x.lim[2] - xx
	if (is.null(y.lim)) {
		if (phyloORclado) {
			if (horizontal) 
				y.lim <- c(1, Ntip)
			else {
				y.lim <- c(0, NA)
				pin2 <- par("pin")[2]
				strWi <- strwidth(x$tip.label, "inches")
				yy.tips <- yy[1:Ntip] * 1.04
				alp <- try(uniroot(function(a) max(a * yy.tips + 
													strWi) - pin2, c(0, 1e+06))$root, silent = TRUE)
				if (is.character(alp)) 
					tmp <- max(yy.tips) * 1.5
				else {
					tmp <- if (show.tip.label) 
								max(yy.tips + strWi/alp)
							else max(yy.tips)
				}
				y.lim[2] <- tmp
			}
		}
		else switch(type, fan = {
						if (show.tip.label) {
							offset <- max(nchar(x$tip.label) * 0.018 * max(yy) * 
											cex)
							y.lim <- c(min(yy) - offset, max(yy) + offset)
						} else y.lim <- c(min(yy), max(yy))
					}, unrooted = {
						if (show.tip.label) {
							offset <- max(nchar(x$tip.label) * 0.018 * max(yy) * 
											cex)
							y.lim <- c(0 - offset, max(yy) + offset)
						} else y.lim <- c(0, max(yy))
					}, radial = {
						if (show.tip.label) {
							offset <- max(nchar(x$tip.label) * 0.03 * cex)
							y.lim <- c(-1 - offset, 1 + offset)
						} else y.lim <- c(-1, 1)
					})
	}
	else if (length(y.lim) == 1) {
		y.lim <- c(0, y.lim)
		if (phyloORclado && horizontal) 
			y.lim[1] <- 1
		if (type %in% c("fan", "unrooted") && show.tip.label) 
			y.lim[1] <- -max(nchar(x$tip.label) * 0.018 * max(yy) * 
							cex)
		if (type == "radial") 
			y.lim[1] <- if (show.tip.label) 
						-1 - max(nchar(x$tip.label) * 0.018 * max(yy) * 
										cex)
					else -1
	}
	if (phyloORclado && direction == "downwards") 
		yy <- y.lim[2] - yy
	if (phyloORclado && root.edge) {
		if (direction == "leftwards") 
			x.lim[2] <- x.lim[2] + x$root.edge
		if (direction == "downwards") 
			y.lim[2] <- y.lim[2] + x$root.edge
	}
	asp <- if (type %in% c("fan", "radial")) 
				1
			else NA
	if (is.null(adj)) 
		adj <- if (phyloORclado && direction == "leftwards") 
					1
				else 0
	if (phyloORclado && show.tip.label) {
		MAXSTRING <- max(strwidth(x$tip.label, cex = cex))
		loy <- 0
		if (direction == "rightwards") {
			lox <- label.offset + MAXSTRING * 1.05 * adj
		}
		if (direction == "leftwards") {
			lox <- -label.offset - MAXSTRING * 1.05 * (1 - adj)
		}
		if (!horizontal) {
			psr <- par("usr")
			MAXSTRING <- MAXSTRING * 1.09 * (psr[4] - psr[3])/(psr[2] - 
						psr[1])
			loy <- label.offset + MAXSTRING * 1.05 * adj
			lox <- 0
			srt <- 90 + srt
			if (direction == "downwards") {
				loy <- -loy
				srt <- 180 + srt
			}
		}
	}
	if (type == "phylogram") {
		phylogram.plot(x$edge, Ntip, Nnode, xx, yy, horizontal, 
				edge.color, edge.width, edge.lty)
	}
	else {
		if (type == "fan") {
			ereorder <- match(z$edge[, 2], x$edge[, 2])
			if (length(edge.color) > 1) {
				edge.color <- rep(edge.color, length.out = Nedge)
				edge.color <- edge.color[ereorder]
			}
			if (length(edge.width) > 1) {
				edge.width <- rep(edge.width, length.out = Nedge)
				edge.width <- edge.width[ereorder]
			}
			if (length(edge.lty) > 1) {
				edge.lty <- rep(edge.lty, length.out = Nedge)
				edge.lty <- edge.lty[ereorder]
			}
			circular.plot(z$edge, Ntip, Nnode, xx, yy, theta, 
					r, edge.color, edge.width, edge.lty)
		}
		else cladogram.plot(x$edge, xx, yy, edge.color, edge.width, 
					edge.lty)
	}
	if (root.edge) 
		switch(direction, rightwards = segments(0, yy[ROOT], 
						x$root.edge, yy[ROOT]), leftwards = segments(xx[ROOT], 
						yy[ROOT], xx[ROOT] + x$root.edge, yy[ROOT]), upwards = segments(xx[ROOT], 
						0, xx[ROOT], x$root.edge), downwards = segments(xx[ROOT], 
						yy[ROOT], xx[ROOT], yy[ROOT] + x$root.edge))
	if (show.tip.label) {
		if (is.expression(x$tip.label)) 
			underscore <- TRUE
		if (!underscore) 
			x$tip.label <- gsub("_", " ", x$tip.label)
		if (phyloORclado) 
			text(xx[1:Ntip] + lox, yy[1:Ntip] + loy, x$tip.label, 
					adj = adj, font = font, srt = srt, cex = cex, 
					col = tip.color)
		if (type == "unrooted") {
			if (lab4ut == "horizontal") {
				y.adj <- x.adj <- numeric(Ntip)
				sel <- abs(XY$axe) > 0.75 * pi
				x.adj[sel] <- -strwidth(x$tip.label)[sel] * 1.05
				sel <- abs(XY$axe) > pi/4 & abs(XY$axe) < 0.75 * 
						pi
				x.adj[sel] <- -strwidth(x$tip.label)[sel] * (2 * 
							abs(XY$axe)[sel]/pi - 0.5)
				sel <- XY$axe > pi/4 & XY$axe < 0.75 * pi
				y.adj[sel] <- strheight(x$tip.label)[sel]/2
				sel <- XY$axe < -pi/4 & XY$axe > -0.75 * pi
				y.adj[sel] <- -strheight(x$tip.label)[sel] * 
						0.75
				text(xx[1:Ntip] + x.adj * cex, yy[1:Ntip] + y.adj * 
								cex, x$tip.label, adj = c(adj, 0), font = font, 
						srt = srt, cex = cex, col = tip.color)
			}
			else {
				adj <- as.numeric(abs(XY$axe) > pi/2)
				srt <- 180 * XY$axe/pi
				srt[as.logical(adj)] <- srt[as.logical(adj)] - 
						180
				for (i in 1:Ntip) text(xx[i], yy[i], cex = cex, 
							x$tip.label[i], adj = adj[i], font = font, 
							srt = srt[i], col = tip.color[i])
			}
		}
		if (type %in% c("fan", "radial")) {
			xx.tips <- xx[1:Ntip]
			angle <- atan2(yy[1:Ntip], xx.tips) * 180/pi
			s <- xx.tips < 0
			angle[s] <- angle[s] + 180
			adj <- numeric(Ntip)
			adj[xx.tips < 0] <- 1
			for (i in 1:Ntip) text(xx[i], yy[i], x$tip.label[i], 
						font = font, cex = cex, srt = angle[i], adj = adj[i], 
						col = tip.color[i])
		}
	}
	if (show.node.label) 
		text(xx[ROOT:length(xx)] + label.offset, yy[ROOT:length(yy)], 
				x$node.label, adj = adj, font = font, srt = srt, 
				cex = cex)
}
######################################################################################
hivc.treeannotator.get.rates<- function(ph, tip.df, nodelabel.idx.edgewidth=5)
{
	tmp				<- as.numeric( sapply( strsplit(ph$node.label,'_'),function(x)	x[nodelabel.idx.edgewidth] ) )	
	rates.df		<- data.table(node=seq_len(Nnode(ph))+Ntip(ph), rate=tmp)
	ans				<- rbind(rates.df, subset(tip.df, select=c(node, rate)))
	#set(ans, which(is.na(ans[,rate])),"rate",mean(ans[,rate],na.rm=1))
	setkey(ans,node)
	ans
}
######################################################################################
hivc.treeannotator.get.phy<- function(ph.beast, beastlabel.idx.clu=1, beastlabel.idx.hivs=4, beastlabel.idx.samplecode=5, beastlabel.idx.rate=6, verbose=1, debug=0)
{
	#	get root height for final tree in calendar time
	ph.tip.ctime	<- sapply(ph.beast, function(x) max( as.numeric( sapply(strsplit(x$tip.label,'_'), function(x)	x[beastlabel.idx.hivs] ) ) ))			
	ph.root.ctime	<- min( sapply(seq_along(ph.beast), function(i)	ph.tip.ctime[i]-max(node.depth.edgelength(ph.beast[[i]]))	) )
	#	extract tip label information
	tip.df			<- lapply(ph.beast, function(x) hivc.treeannotator.tiplabel2df(x, beastlabel.idx.clu=beastlabel.idx.clu, beastlabel.idx.hivs=beastlabel.idx.hivs, beastlabel.idx.samplecode=beastlabel.idx.samplecode, beastlabel.idx.rate=beastlabel.idx.rate) )
	#	prepare cluster subtrees
	clu.subtrees	<- lapply( seq_along(ph.beast), function(i)
			{				
				x<- ph.beast[[i]]
				#	convert heights into calendar time and collapse node.label
				tmp					<- ph.tip.ctime[i]								
				#subset(x$node.label, node==147, select=c(node, height_median, height_95_HPD_MIN, height_95_HPD_MAX, posterior))
				tmp					<- x$node.label[, list(	node=node, 
															height_median=ph.tip.ctime[i]-height_median, 
															height_95_HPD_MIN=ph.tip.ctime[i]-height_95_HPD_MAX, 
															height_95_HPD_MAX=ph.tip.ctime[i]-height_95_HPD_MIN, 
															rate_median=rate_median,
															posterior=posterior)]			
				tmp					<- tmp[, list(node.label= paste(posterior,height_median, height_95_HPD_MIN, height_95_HPD_MAX, rate_median, sep='_')),by="node"]
				x$node.label		<- tmp[, node.label]
				x$node.label.format	<- "posterior height_median height_95_HPD_MIN height_95_HPD_MAX rate_median"				
				#
				#	extract rooted ExaML clusters
				#								
				x$tip.label	<- tip.df[[i]][, FASTASampleCode]
				tmp			<- tip.df[[i]][, list(node=getMRCA(x,FASTASampleCode)),by=cluster]
				clu.subtrees<- lapply(tmp[,node], function(z)
						{	
							ans						<- extract.clade(x, z, root.edge= 1, interactive = FALSE)
							ans$root.edge			<- as.numeric(strsplit(ans$node.label[1],'_')[[1]][2])-ph.root.ctime		#reset root edge against root of all runs combined
							ans$node.label.format	<- x$node.label.format
							ans
						})	
				names(clu.subtrees)<- tmp[,cluster]	
				clu.subtrees
			})
	clu.subtrees	<- eval(parse(text= paste("c(",paste('clu.subtrees[[',seq_along(clu.subtrees),']]', sep='',collapse=','),")",sep='') ))
	if(verbose)	cat(paste("\nFound ExaML clusters in treeannotator files, number of clusters is n=", length(clu.subtrees) ))
	if(debug)
		clu.subtrees	<- lapply(1:3, function(i) clu.subtrees[[i]] )
	#	join all clusters 
	cluphy				<- eval(parse(text=paste('clu.subtrees[[',seq_along(clu.subtrees),']]', sep='',collapse='+')))	
	if(verbose)	cat(paste("\nFound ExaML clusters in treeannotator files, number of sequences is n=", Ntip(cluphy) ))
	#	retain tip info for those tip labels in cluphy
	tip.df			<- rbindlist(tip.df)	
	tip.df			<- merge(data.table(FASTASampleCode=cluphy$tip.label), tip.df, by="FASTASampleCode")
	tip.df			<- cbind(tip.df, node=seq_len(Ntip(cluphy)))

	list(cluphy=cluphy, ph.tip.df=tip.df, ph.tip.ctime=ph.tip.ctime, ph.root.ctime=ph.root.ctime)
}
######################################################################################
hivc.treeannotator.get.tmrcas<- function(ph.beast, beastlabel.idx.hivs=4)
{
	#	for each tree, return 	height_median height_95_HPD_MIN height_95_HPD_MAX for the parent of each tip 		
	ph.trmca	<- lapply(ph.beast, function(x)
			{
				#select heights and convert heights into calendar time
				tip.latest			<- max( as.numeric( sapply(strsplit(x$tip.label,'_'), function(x)	x[beastlabel.idx.hivs] ) ) )
				tip.parents			<- unique( x$edge[x$edge[,2]<=Ntip(x),1] )		
				ans					<- x$node.label[J(tip.parents)][, list(node=node, height_median=tip.latest-height_median, height_95_HPD_MIN=tip.latest-height_95_HPD_MAX, height_95_HPD_MAX=tip.latest-height_95_HPD_MIN)]
				tmp					<- sapply(ans[,node], function(z)	x$edge[ x$edge[,1]==z, 2 ] )
				tmp[tmp>Ntip(x)]	<- NA
				tmp					<- t( apply(tmp,2,function(z) x$tip.label[ sort(z,na.last=1) ]) )
				ans[,height_95_diff:= height_95_HPD_MAX-height_95_HPD_MIN]
				ans[,tip1:= tmp[,1]]		
				ans[,tip2:= tmp[,2]]					
				ans
			})
	ph.trmca	<- rbindlist( ph.trmca )
}
######################################################################################
hivc.treeannotator.tiplabel2df<- function(x, beastlabel.idx.clu=1, beastlabel.idx.hivn=2, beastlabel.idx.hivd=3, beastlabel.idx.hivs=4, beastlabel.idx.samplecode=5, beastlabel.idx.rate=6)
{
	tmp		<- t( sapply(strsplit(x$tip.label,'_'), function(z)	z[c(beastlabel.idx.clu,beastlabel.idx.hivn, beastlabel.idx.hivd, beastlabel.idx.hivs, beastlabel.idx.samplecode, beastlabel.idx.rate)] ) )
	data.table(cluster=tmp[,1], NegT= tmp[,2], AnyPos_T1= tmp[,3], TipT= tmp[,4], FASTASampleCode=tmp[,5], rate=tmp[,6] )				
}	
######################################################################################
hivc.treeannotator.get.clusterprob<- function(ph.beast, beastlabel.idx.clu=1, beastlabel.idx.samplecode=5, verbose=1)
{
	#	for each of the clusters, compute the posterior probability of a common MRCA
	clu.df	<- lapply(ph.beast, function(x)
			{							
				tmp		<- data.table(	tip=seq_len(Ntip(x)), 
										cluster=as.numeric( sapply( strsplit(x$tip.label,'_'),function(z)  z[beastlabel.idx.clu] ) ),
										FASTASampleCode=sapply( strsplit(x$tip.label,'_'),function(z)  z[beastlabel.idx.samplecode] )										
										)
				ans		<- merge( tmp[, list(node=hivc.clu.mrca(x, x.tip=tip), FASTASampleCode=FASTASampleCode), by=cluster], subset(x$node.label, select=c(node, posterior)), by="node" )
				subset(ans, select=c(cluster, FASTASampleCode, posterior))
			})
	clu.df	<- rbindlist(clu.df)		
	if(verbose)	cat(paste("\nRange of posterior probabilities that each of the putative clusters each has a common MRCA, min=",min(clu.df[,posterior])," max=",max(clu.df[,posterior]) ))
	clu.df
}
######################################################################################
hivc.treeannotator.get.edgewidth<- function(ph, rates.df, scale.edgewidth= 12)
{
	if(is.null(rates.df))
		rates.df<- data.table(node=Nnode(ph,internal.only=0), rate=NA)
	
	edge.width	<- merge(rates.df, data.table(node=ph$edge[,2], edge=seq_len(nrow(ph$edge))), all.y=1, by="node")
	edge.width[, width:=edge.width[,rate]/mean(edge.width[,rate], na.rm=1)]
	set(edge.width, NULL, "width", (edge.width[,width]-1)*scale.edgewidth + 1)
	set(edge.width, which(edge.width[,is.na(width)]), "width", 1)
	set(edge.width, which(edge.width[,width<0.1]), "width", 0.1)
	setkey(edge.width, edge)
	edge.width
}
######################################################################################
#	ph<- cluphy; end.ctime=2013.3; cex.nodelabel=0.5; cex.tiplabel=0.5; file=NULL; pdf.width=7; pdf.height=20
hivc.treeannotator.plot<- function(ph, ph.root.ctime, youngest.tip.ctime, df.all, df.viro, df.immu, df.treatment=NULL, df.tstem=NULL, df.rates=NULL, end.ctime=2013.3, cex.nodelabel=0.5, cex.tiplabel=0.5, file=NULL, pdf.width=7, pdf.height=20)
{		
	require(RColorBrewer)
	if(class(file)=="character")
		pdf(file, width=pdf.width, height=pdf.height)
	par(mar=c(0,0,0,0))
	
	cols			<- brewer.pal(12,"Paired")
	ph.xlim			<- end.ctime-ph.root.ctime+ c(-22,6)
	ph.ylim			<- c(1,Ntip(ph)) + c(-1,1)			
	
	plot(ph, x.lim=ph.xlim, y.lim= ph.ylim, show.tip.label=0, edge.color = 0, tip.color = 0)
	# add calendar timeline			
	hivc.treeannotator.plot.ctimeline(ph, youngest.tip.ctime, end.ctime, add.yinch= 0.5)
	# add BEAST TMRCA 95% credibility interval
	ph.nodeheighthpds	<- hivc.treeannotator.nodelabels.getnodeheightHPD(ph, youngest.tip.ctime)
	#do not plot the BEAST TMRCA 95% credibility interval for those nodes for which more data than the 95% interval is available
	if(!is.null(df.tstem))	
		ph.nodeheighthpds	<- subset( ph.nodeheighthpds, !node%in%unique(df.tstem[,mrca]) )
	hivc.treeannotator.plot.hpdbars(ph, ph.nodeheighthpds, col=cols[1], lwd=4)
	# add NegT and AnyPos_T1
	ph.seronodeheight	<- hivc.treeannotator.sero.getnodeheight.range(ph, df.all, youngest.tip.ctime)
	hivc.treeannotator.plot.seronodeheightrange(ph, ph.seronodeheight, add.yinch= -0.03, width.yinch= 0.03, width.yinch.past.AnyPos_T1= 0, col=cols[2])
	# add lRNA timeline
	ph.viro.timeline	<- hivc.treeannotator.get.viro.timeline(ph, df.all, df.viro, youngest.tip.ctime, df.treatment=df.treatment)
	hivc.treeannotator.plot.viro.timeline(ph, ph.viro.timeline, viro.min= log10(300), width.yinch= 0.15, add.yinch= 0.005, col.bg= cols[c(5,10,12)], col.legend= cols[6], cex.txt= 0.2)
	# add BEAST posterior density of TMRCAs where available
	if(!is.null(df.tstem))
		hivc.treeannotator.plot.tipstem.timeline(ph, youngest.tip.ctime, df.tstem, width.yinch=0.1, add.yinch=0.005, col.bg=cols[1] )		
	# add CD4 timeline
	ph.immu.timeline	<- hivc.treeannotator.get.immu.timeline(ph, df.all, df.immu, youngest.tip.ctime, end.ctime=2013.3)
	hivc.treeannotator.plot.immu.timeline(ph, ph.immu.timeline, immu.min= 150, immu.max= 800, width.yinch= 0.15, add.yinch= -0.005, col.bg= cols[3], col.legend= cols[4], cex.txt= 0.2)
	# re-plot phylogeny
	ph$node.label		<- as.numeric(sapply( strsplit( ph$node.label, '_' ), function(x)	x[1] ))
	edge.width			<- hivc.treeannotator.get.edgewidth(ph, df.rates, scale.edgewidth= 8)
	hivc.phy.plotupon(ph, show.tip.label=0, show.node.label=1, cex=cex.nodelabel, edge.width=edge.width[,width],)
	# add rate labels			
	if(!is.null(df.rates))	
		hivc.treeannotator.plot.rates(ph, edge.width, add.xinch=-0.1, cex.rate=0.3)
	# add tip labels			
	ph.tiplabel			<- hivc.clu.get.tiplabels(ph, 	df.all, col.notmsm="#4EB3D3", col.Early="#EF9708", col.highVL="#FEE391", col.AfterTreat="#D4B9DA", col.green="#D9F0A3", col.latePres="#FA9FB5", select=c("CountryInfection","Trm","Sex","isAcute","lRNA.early","Patient","RegionHospital") )				
	tmp					<- rep( max(node.depth.edgelength(ph)) - (youngest.tip.ctime-ceiling(end.ctime)), Ntip(ph))
	hivc.clu.plot.tiplabels(seq_len(Ntip(ph)), ph.tiplabel$text, ph.tiplabel$col, xx=tmp, adj = c(-0.05, 0.5), cex=cex.tiplabel, add.xinch= 0.03, add.yinch= 0.02)
	# add legend	
	legend("topright", fill= cols[c(1,2,3,5,10,12)], legend=c("BEAST 95% TMRCA", "interval [last HIV-, diagnosis]", "CD4 timeline", "VL timeline", "VL timeline under treatment", "VL timeline under treatment"), bty='n', border=NA, cex=cex.tiplabel)
	
	if(class(file)=="character")
		dev.off()				
}
######################################################################################
hivc.treeannotator.plot.rates<- function(ph, edge.width, add.xinch=-0.1, cex.rate=0.5)
{
	lastPP 	<- get("last_plot.phylo", envir = .PlotPhyloEnv)	
	tmp	<- edge.width[, list(xx= lastPP$xx[ ph$edge[edge,1] ]+xinch(add.xinch), yy= mean( lastPP$yy[ ph$edge[edge,] ] ), rate=rate), by="edge"]
	text(tmp[,xx], tmp[,yy], tmp[,rate], cex=cex.rate)
}
######################################################################################
hivc.treeannotator.plot.tipstem.timeline<- function(ph, youngest.tip.ctime, df.tstem, width.yinch=0.15, add.yinch=0, col.bg="grey75", density=30, lwd=0.5)
{
	lastPP 	<- get("last_plot.phylo", envir = .PlotPhyloEnv)
	set(df.tstem, NULL, "tstem", youngest.tip.ctime - df.tstem[,tstem])
	set(df.tstem, NULL, "tstem", max(node.depth.edgelength(ph)) - df.tstem[,tstem])
	df.tstem[, yyi:= density]
	setkey(df.tstem, tip)
	df.tstem<- merge(df.tstem, df.tstem[,	list(scale=yinch(width.yinch)/max(yyi))	,by="tip"], by="tip")	
	set(df.tstem, NULL, "yyi",  df.tstem[,yyi*scale])
	df.tstem[, yy:= yinch(add.yinch)+lastPP$yy[df.tstem[,tip]]]
	df.tstem[, col:=col.bg]
	setkey(df.tstem, tip)
	dummy<- sapply( unique(df.tstem[,tip]), function(x)
			{
				z		<- df.tstem[J(x)]
				setkey(z, tstem)
				polygon( c( z[,tstem],z[nrow(z),tstem] ), z[1,yy]+c( z[,yyi],0 ), border=z[1,col], col=z[1,col], density=density, lwd=lwd )  				
			})
}
######################################################################################
hivc.treeannotator.plot.ctimeline<- function(ph, youngest.tip.ctime, end.ctime, col.bg= c(my.fade.col("black",0.15),"transparent"), col.txt= c(my.fade.col("black",1),"transparent"), cex.txt= 0.5, add.yinch= 0.5)
{
	lastPP 	<- get("last_plot.phylo", envir = .PlotPhyloEnv)	
	
	tmp 	<- seq( max(lastPP$xx) - ( youngest.tip.ctime-floor(youngest.tip.ctime) )+1, min(lastPP$xx), -1 )
	df.time	<- data.table(ctime= seq( floor(youngest.tip.ctime), by=-1, len= length(tmp) ), xx.u=tmp) 
	if(1+floor(youngest.tip.ctime)<floor(end.ctime))
	{
		tmp		<- seq( 1+floor(youngest.tip.ctime),floor(end.ctime),by=1 )
		tmp		<- data.table( ctime=rev(tmp), xx.u=rev(seq(df.time[1,xx.u]+1, len=length(tmp), by=1)))
		df.time	<- rbind(tmp, df.time)
	}
	df.time	<- cbind( 	df.time[-nrow(df.time), ], 
						data.table(	xx.l	= df.time[-1,xx.u], 
									col.bg	= rep(col.bg,ceiling(nrow(df.time)/length(col.bg)))[seq_len(nrow(df.time)-1)],
									col.txt = rep(col.txt,ceiling(nrow(df.time)/length(col.txt)))[seq_len(nrow(df.time)-1)]		) )
	
	rect(df.time[,xx.l], lastPP$y.lim[1]-yinch(add.yinch), df.time[,xx.u], lastPP$y.lim[2]+yinch(add.yinch), col=df.time[,col.bg], border=NA)				
	text(df.time[,xx.l+(xx.u-xx.l)/2.1], lastPP$y.lim[1]-yinch(add.yinch)/4, df.time[,ctime], cex=cex.txt, col=df.time[,col.txt], offset=0)
	text(df.time[,xx.l+(xx.u-xx.l)/2.1], lastPP$y.lim[2]+yinch(add.yinch)/4, df.time[,ctime], cex=cex.txt, col=df.time[,col.txt], offset=0)
}
######################################################################################
hivc.treeannotator.plot.immu.timeline<- function(ph, ph.immu.timeline, immu.min= 150, immu.max= 800, immu.legend= c(200, 350, 500, immu.max), width.yinch= 0.15, add.yinch= -0.005, col.bg= cols[3], col.legend= cols[4], cex.txt= 0.2)
{
	lastPP 	<- get("last_plot.phylo", envir = .PlotPhyloEnv)
	set(ph.immu.timeline, NULL, "PosCD4", max(node.depth.edgelength(ph)) - ph.immu.timeline[,PosCD4])
	ph.immu.timeline[, yyCD4:= CD4]
	set(ph.immu.timeline, which(ph.immu.timeline[,yyCD4]>immu.max), "yyCD4", immu.max)
	set(ph.immu.timeline, NULL, "yyCD4", ph.immu.timeline[,yyCD4]-immu.min)
	set(ph.immu.timeline, which(ph.immu.timeline[,yyCD4]<0), "yyCD4", 0.)			
	scale	<- yinch(width.yinch) / max( ph.immu.timeline[,yyCD4])
	set(ph.immu.timeline, NULL, "yyCD4", ph.immu.timeline[,yyCD4] * scale )
	ph.immu.timeline[, yy:= yinch(add.yinch)+lastPP$yy[ph.immu.timeline[,tip]]]
	
	dummy<- sapply( unique(ph.immu.timeline[,tip]), function(x)
			{
				z<- ph.immu.timeline[J(x)]
				polygon( c( z[,PosCD4], z[nrow(z),PosCD4], z[1,PosCD4] ), c( z[,yy-yyCD4], z[nrow(z),yy], z[1,yy] ), border=NA, col=col.bg	)
				sapply(z[1,yy]-(immu.legend-immu.min)*scale,function(i)		lines(z[c(1,nrow(z)),PosCD4], rep(i,2), col=col.legend, lty=3, lwd=0.2)		)		
				text(rep(z[nrow(z),PosCD4],3),z[1,yy]-(immu.legend-immu.min)*scale,immu.legend,cex=cex.txt, col=col.legend)
			})
}
######################################################################################
hivc.treeannotator.get.immu.timeline<- function(ph, df, df.immu, youngest.tip.ctime, end.ctime=2013.3)
{				
	setkey(df, FASTASampleCode)
	tmp				<- df[J(ph$tip.label)][, list(tip=seq_along(ph$tip.label), FASTASampleCode=FASTASampleCode, Patient=Patient, DateDied=DateDied)]		
	ans				<- merge(subset(df.immu, select=c(Patient, PosCD4, CD4)), tmp, all.y=1, by="Patient")
	set(ans,NULL,"PosCD4",	youngest.tip.ctime - hivc.db.Date2numeric(ans[,PosCD4]))
	set(ans,NULL,"DateDied",	youngest.tip.ctime - hivc.db.Date2numeric(ans[,DateDied]))				
	set(ans,which(is.na(ans[,DateDied])), "DateDied", youngest.tip.ctime - end.ctime)				
	ans				<- ans[,list(PosCD4=c(PosCD4,DateDied[1]), CD4=c(CD4,tail(CD4,1))),by="tip"]
	setkey(ans,tip)
	ans
}
######################################################################################
#	viro.min= log10(300); width.yinch= 0.15; add.yinch= 0.005; col.bg= cols[c(5,9,10)]; col.legend= cols[6]; cex.txt= 0.2
hivc.treeannotator.plot.viro.timeline<- function(ph, ph.viro.timeline, viro.min= log10(300), viro.legend=c(3,4,5,6), width.yinch= 0.2, add.yinch= 0.005, col.bg= "red", col.legend="red", cex.txt= 0.2)
{	
	lastPP 	<- get("last_plot.phylo", envir = .PlotPhyloEnv)
	set(ph.viro.timeline, NULL, "PosRNA", max(node.depth.edgelength(ph)) - ph.viro.timeline[,PosRNA])
	ph.viro.timeline[, yylRNA:= lRNA]
	set(ph.viro.timeline, NULL, "yylRNA", ph.viro.timeline[,yylRNA]-viro.min)
	scale	<- yinch(width.yinch) / max( ph.viro.timeline[,yylRNA])
	set(ph.viro.timeline, NULL, "yylRNA", ph.viro.timeline[,yylRNA] * scale )
	ph.viro.timeline[, yy:= yinch(add.yinch)+lastPP$yy[ph.viro.timeline[,tip]]]
	#define treatment periods if there
	if(!any(ph.viro.timeline[, NoDrug!=0]))
		ph.viro.timeline[, col:=col.bg[1]]
	else
	{
		if(length(col.bg)==1)	stop("for TPeriod != NA, expect more than one 'col.bg'")
		col.bg.nNA	<- rep(col.bg[-1], ceiling(max(ph.viro.timeline[, TPeriod])/(length(col.bg)-1)))
		ph.viro.timeline[, col:=""]
		set(ph.viro.timeline, which(ph.viro.timeline[,NoDrug==0]), "col", col.bg[1])
		set(ph.viro.timeline, which(ph.viro.timeline[,NoDrug!=0]), "col", col.bg.nNA[ subset(ph.viro.timeline, NoDrug!=0)[,TPeriod] ])		
	}	
	setkey(ph.viro.timeline, tip)
	dummy<- sapply( unique(ph.viro.timeline[,tip]), function(x)
			{
				z		<- ph.viro.timeline[J(x)]	
				#reset TPeriod because there can be multiple off treatment periods (ie TPeriod==0) and we cannot lump them together in the next line
				#NOT NEEDED ANY LONGER	set(z, NULL, "TPeriod",cumsum(c(0,as.numeric(abs(diff(z[,TPeriod]))>0))))
				dummy	<- z[,	{								
									polygon( c( PosRNA, PosRNA[length(PosRNA)], PosRNA[1] ), c( yylRNA+yy, yy[length(yy)], yy[1] ), border=NA, col=col[1]	)	
							}, by="TPeriod"]								
				sapply(z[1,yy]+(viro.legend-viro.min)*scale,function(i)		lines(z[c(1,nrow(z)),PosRNA], rep(i,2), col=col.legend, lty=3, lwd=0.2)		)		
				text(rep(z[nrow(z),PosRNA],3),z[1,yy]+(viro.legend-viro.min)*scale,paste("1e",viro.legend,sep=''),cex=cex.txt, col=col.legend)
				#stop()
			})
}
######################################################################################
hivc.treeannotator.get.viro.timeline<- function(ph, df, df.viro, youngest.tip.ctime, df.treatment=NULL, end.ctime=2013.3)
{
	setkey(df, FASTASampleCode)
	
	if(is.null(df.treatment))		#prepare a single viral load timeline without treatment periods
	{
		tmp				<- df[J(ph$tip.label)][, list(tip=seq_along(ph$tip.label), FASTASampleCode=FASTASampleCode, Patient=Patient, DateDied=DateDied)]		
		ans				<- merge(subset(df.viro, select=c(Patient, PosRNA, lRNA)), tmp, all.y=1, by="Patient")
		set(ans,NULL,"PosRNA",	youngest.tip.ctime - hivc.db.Date2numeric(ans[,PosRNA]))
		set(ans,NULL,"DateDied",	youngest.tip.ctime - hivc.db.Date2numeric(ans[,DateDied]))				
		set(ans,which(is.na(ans[,DateDied])), "DateDied", youngest.tip.ctime - end.ctime)				
		ans				<- ans[,list(PosRNA=c(PosRNA,DateDied[1]), lRNA=c(lRNA,tail(lRNA,1), TPeriod=0, NoDrug=0)),by="tip"]
		setkey(ans,tip)
	}
	else							#prepare a single viral load timeline with treatment periods
	{
		#as above except TPeriod=NA
		tmp				<- df[J(ph$tip.label)][, list(tip=seq_along(ph$tip.label), FASTASampleCode=FASTASampleCode, Patient=Patient, DateDied=DateDied)]		
		ans				<- merge(subset(df.viro, select=c(Patient, PosRNA, lRNA)), tmp, all.y=1, by="Patient")
		set(ans,NULL,"PosRNA",	youngest.tip.ctime - hivc.db.Date2numeric(ans[,PosRNA]))
		set(ans,NULL,"DateDied",	youngest.tip.ctime - hivc.db.Date2numeric(ans[,DateDied]))				
		set(ans,which(is.na(ans[,DateDied])), "DateDied", youngest.tip.ctime - end.ctime)				
		ans				<- ans[,list(Patient= rep(Patient[1], length(Patient)+1), PosRNA=c(PosRNA,DateDied[1]), lRNA=c(lRNA,tail(lRNA,1))),by="tip"]
		#prepare subset of df.treatment
		tmp				<- merge(subset(df.treatment, select=c(Patient, StartTime, StopTime, NoDrug)), unique(subset(ans,select=Patient)), all.y=1, by="Patient")	
		set(tmp,NULL,"StartTime",	youngest.tip.ctime - hivc.db.Date2numeric(tmp[,StartTime]))
		set(tmp,NULL,"StopTime",	youngest.tip.ctime - hivc.db.Date2numeric(tmp[,StopTime]))
		set(tmp, which( tmp[,StopTime==min(StopTime, na.rm=1)] ), "StopTime", youngest.tip.ctime - 2013.3)
		ans				<- merge(tmp, ans, by="Patient", allow.cartesian=1)
		tmp				<- which(ans[,is.na(StartTime)])
		set(ans,tmp,"StartTime",youngest.tip.ctime - 2013.3)
		set(ans,tmp,"StopTime",youngest.tip.ctime - 2013.3)
		set(ans,tmp,"NoDrug",0L)
		#now have all treatments and viral loads measures togethers
		ans					<- ans[,	{
											x		<- data.table(StartTime, StopTime, NoDrug, tip, PosRNA,  lRNA)											
											#select outside any treatment period
											tmp		<- subset(x, PosRNA>max(StartTime))								#handle no drug before first treatment
											tmp		<- rbind(tmp, subset(x, PosRNA<min(StopTime)) )					#handle no drug after stop treatment
											if(nrow(tmp))
											{
												setkey(tmp, PosRNA)
												tmp		<- unique(tmp)
												set(tmp,NULL,"NoDrug",0L)
												#select only those viral loads within a particular treatment period
												tmp		<- rbind(tmp, subset(x, StartTime>=PosRNA & PosRNA>=StopTime))												
											}
											else
												tmp		<- subset(x, StartTime>=PosRNA & PosRNA>=StopTime)
											#assign integer value to different treatment periods
											setkey(tmp, NoDrug, StartTime)
											tmp2	<- subset(unique(tmp), select=c(StartTime, NoDrug))[order(-StartTime, NoDrug)]
											tmp2	<- tmp2[, list(StartTime=StartTime, NoDrug=NoDrug, TPeriod=seq_along(NoDrug))]
											tmp		<- merge(tmp, tmp2, all.x=1, by=c("StartTime", "NoDrug"))											
											tmp		<- subset(tmp[order(-PosRNA)], select=c(PosRNA,  lRNA, NoDrug, TPeriod))
											#add endpoints for each treatment period
											tmp2	<- subset(tmp, TPeriod>min(TPeriod))[,  list(PosRNA=PosRNA[1], lRNA=lRNA[1] ),by=TPeriod]
											set(tmp2, NULL, "TPeriod", tmp2[,TPeriod]-1)
											tmp2	<- merge(tmp2, unique(subset(tmp, select=c(TPeriod, NoDrug))), by="TPeriod")
											tmp		<- rbind(tmp, subset(tmp2, select=c(PosRNA,  lRNA, NoDrug, TPeriod)))
											tmp
										},by="tip"]
		ans					<- ans[order(tip, -PosRNA, TPeriod)]		
	}
	ans
}
######################################################################################
hivc.treeannotator.sero.getnodeheight.range<- function(ph, df, youngest.tip.ctime)
{
	setkey(df, FASTASampleCode)
	ans					<- cbind( data.table(tip=seq_along(ph$tip.label)), subset(df[J(ph$tip.label)], select=c(NegT, AnyPos_T1, DateDied)) )
	set(ans,NULL,"NegT",		youngest.tip.ctime - hivc.db.Date2numeric(ans[,NegT]))
	set(ans,NULL,"AnyPos_T1",	youngest.tip.ctime - hivc.db.Date2numeric(ans[,AnyPos_T1]))				
	set(ans,NULL,"DateDied",	youngest.tip.ctime - hivc.db.Date2numeric(ans[,DateDied]))
	set(ans,which(is.na(ans[,DateDied])), "DateDied", 0.)
	ans
}
######################################################################################
hivc.treeannotator.plot.seronodeheightrange<- function(ph, ph.seronodeheight, add.yinch= -0.05, width.yinch= 0.1, width.yinch.past.AnyPos_T1= 0.02, col="red")
{						
	lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
	if (!lastPP$use.edge.length) 
		stop("function needs edge length information")
	if (lastPP$type != "phylogram") 
		stop("currently only 'type == phylogram' supported")
	if (lastPP$dir!="rightwards")
		stop("currently only rightwards supported")			
	
	
	tmp				<- max(node.depth.edgelength(ph)) - subset(ph.seronodeheight, select=c(NegT, AnyPos_T1, DateDied))		#from node heights to root heights, assuming root is plotted at 0
	yy.l 			<- lastPP$yy[ph.seronodeheight[,tip]] + yinch(add.yinch)			
	rect(tmp[,NegT], yy.l, tmp[,AnyPos_T1], yy.l+yinch(width.yinch), col = col, border=NA)			
	rect(tmp[,AnyPos_T1], yy.l, tmp[,DateDied], yy.l+yinch(width.yinch.past.AnyPos_T1), col=col, border=NA)	
}
######################################################################################
hivc.treeannotator.nodelabels.getnodeheightHPD<- function(ph, youngest.tip.ctime, node.label.hpd.l= 3, node.label.hpd.u= 4)
{
	nodes			<- which(!is.na(ph$node.label))
	node.hpd		<- t( sapply(strsplit(ph$node.label[nodes],'_'), function(x)	as.numeric(x[c(node.label.hpd.l, node.label.hpd.u)])) )
	node.hpd		<- youngest.tip.ctime - node.hpd				#from calendar time to node heights
	data.table(node=nodes, hpd.l=node.hpd[,1], hpd.u=node.hpd[,2])
}
######################################################################################
hivc.treeannotator.plot.hpdbars<- function(ph, ph.nodeheighthpds, col="grey75", lwd=4 ) 
{
	lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
	if (!lastPP$use.edge.length) 
		stop("function needs edge length information")
	if (lastPP$type != "phylogram") 
		stop("currently only 'type == phylogram' supported")
	if (lastPP$dir!="rightwards")
		stop("currently only rightwards supported")				
	tmp				<- max(node.depth.edgelength(ph)) - subset(ph.nodeheighthpds, select=c(hpd.l,hpd.u))		#from node heights to root heights, assuming root is plotted at 0
	yy 				<- lastPP$yy[Ntip(ph)+ph.nodeheighthpds[,node]]
	segments(tmp[,hpd.l], yy, tmp[,hpd.u], yy, col = col, lwd = lwd)					
}	
######################################################################################
#' @export
or.dist.dna<- function (x, model = "K80", variance = FALSE, gamma = FALSE, 
		pairwise.deletion = FALSE, base.freq = NULL, as.matrix = FALSE) 
{
	MODELS <- c("RAW", "JC69", "K80", "F81", "K81", "F84", "T92", 
			"TN93", "GG95", "LOGDET", "BH87", "PARALIN", "N", "TS", 
			"TV", "INDEL", "INDELBLOCK")
	imod <- pmatch(toupper(model), MODELS)
	if (is.na(imod)) 
		stop(paste("'model' must be one of:", paste("\"", MODELS, 
								"\"", sep = "", collapse = " ")))
	if (imod == 11 && variance) {
		warning("computing variance not available for model BH87")
		variance <- FALSE
	}
	if (gamma && imod %in% c(1, 5:7, 9:17)) {
		warning(paste("gamma-correction not available for model", 
						model))
		gamma <- FALSE
	}
	if (is.list(x)) 
		x <- as.matrix(x)
	nms <- dimnames(x)[[1]]
	n <- dim(x)
	s <- n[2]
	n <- n[1]
	if (imod %in% c(4, 6:8)) {
		BF <- if (is.null(base.freq)) 
					base.freq(x)
				else base.freq
	}
	else BF <- 0
	if (imod %in% 16:17) 
		pairwise.deletion <- TRUE
	if (!pairwise.deletion) {
		keep <- .C("GlobalDeletionDNA", x, n, s, rep(1L, s), 
				PACKAGE = "ape")[[4]]
		x <- x[, as.logical(keep)]
		s <- dim(x)[2]
	}
	Ndist <- if (imod == 11) 
				n * n
			else n * (n - 1)/2
	var <- if (variance) 
				double(Ndist)
			else 0
	if (!gamma) 
		gamma <- alpha <- 0
	else alpha <- gamma <- 1
print(x)	
	d <- .C("dist_dna", x, n, s, imod, double(Ndist), BF, as.integer(pairwise.deletion), 
			as.integer(variance), var, as.integer(gamma), alpha, 
			DUP = FALSE, NAOK = TRUE, PACKAGE = "ape")
	if (variance) 
		var <- d[[9]]
	d <- d[[5]]
	if (imod == 11) {
		dim(d) <- c(n, n)
		dimnames(d) <- list(nms, nms)
	}
	else {
		attr(d, "Size") <- n
		attr(d, "Labels") <- nms
		attr(d, "Diag") <- attr(d, "Upper") <- FALSE
		attr(d, "call") <- match.call()
		attr(d, "method") <- model
		class(d) <- "dist"
		if (as.matrix) 
			d <- as.matrix(d)
	}
	if (variance) 
		attr(d, "variance") <- var
	d
}

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

hivc.seq.length<- function(seq.DNAbin.mat, exclude=c('-','?'))
{
	counts	<- apply(seq.DNAbin.mat,1,function(x) base.freq(x, freq=1, all=1))
	apply(counts[ !rownames(counts)%in%exclude, ],2,sum)
}

hivc.seq.proportion.ambiguous<- function(seq.DNAbin.mat, exclude=c('-','?'))
{
	counts	<- apply(seq.DNAbin.mat,1,function(x) base.freq(x, freq=1, all=1))
	len		<- apply(counts[ !rownames(counts)%in%exclude, ],2,sum)
	pa		<- apply(counts[c("r", "m", "w", "s", "k", "y", "v", "h", "d", "b"),],2,sum)
	pa/len
}

hivc.seq.gc.content<- function(seq.DNAbin.mat)
{	
	rna.gc.fraction.n		<- c('a','c','g','t',	'r','m','w','s',	'k','y','v','h',		'd','b','n','-','?')
	rna.gc.fraction			<- c( 0, 1, 1, 0,		0.5, 0.5, 0, 1, 	1/2, 1/2, 2/3, 1/3,		1/3,2/3, 1/4, 0, 0)		#this fraction assumes that U is synonymous with T and that U does not occur in the code
	rna.gc.sum				<- c( T, T, T, T,       T, T, T, T,         T, T, T, T,				T, T, T, F, F )	
	counts					<- apply(seq.DNAbin.mat,1,function(x) base.freq(x, freq=1, all=1))
	apply(counts*rna.gc.fraction,2,sum) / apply(counts[rna.gc.sum,],2,sum)
}

#slight modification of blastSequences() in pkg annotate
hivc.seq.blast<- function (x, database = "nr", hitListSize = "10", filter = "L", expect = "10", program = "blastn", organism= "HIV-1") 
{
	baseUrl <- "http://www.ncbi.nlm.nih.gov/blast/Blast.cgi"
	query <- paste("QUERY=", as.character(x), "&DATABASE=", database, "&ORGANISM=", organism,
			"&HITLIST_SIZE=", hitListSize, "&FILTER=", filter, "&EXPECT=", 
			expect, "&PROGRAM=", program, sep = "")
	url0 <- sprintf("%s?%s&CMD=Put", baseUrl, query)
	results <- tempfile()
	Sys.sleep(5)
	require(XML)
	post <- htmlTreeParse(url0, useInternalNodes = TRUE)
	x <- post[["string(//comment()[contains(., \"QBlastInfoBegin\")])"]]
	rid <- sub(".*RID = ([[:alnum:]]+).*", "\\1", x)
	rtoe <- as.integer(sub(".*RTOE = ([[:digit:]]+).*", "\\1", 
					x))
	url1 <- sprintf("%s?RID=%s&FORMAT_TYPE=XML&CMD=Get", baseUrl, 
			rid)
	Sys.sleep(rtoe)
	result <- annotate:::.tryParseResult(url1)
	qseq <- xpathApply(result, "//Hsp_qseq", xmlValue)
	hseq <- xpathApply(result, "//Hsp_hseq", xmlValue)
	require(Biostrings)
	res <- list()
	for (i in seq_len(length(qseq))) {
		res[i] <- DNAMultipleAlignment(c(hseq[[i]], qseq[[i]]), 
				rowmask = as(IRanges(), "NormalIRanges"), colmask = as(IRanges(), 
						"NormalIRanges"))
	}
	res
}

#slight modification of read.blast() in pkg RFLPtools; expects blast was run with -outfmt 6
hivc.seq.blast.read<- function (file, sep = "\t") 
{
	require(data.table)
	x <- read.table(file = file, header = FALSE, sep = sep, quote = "\"", dec = ".", fill = TRUE, comment.char = "", stringsAsFactors = FALSE)
	if (ncol(x) != 12) 
		stop("Data in given file", basename(file), "has wrong dimension!")
	names(x) <- c("query.id", "subject.id", "identity", "alignment.length","mismatches", "gap.opens", "q.start", "q.end", "s.start","s.end", "evalue", "bit.score")
	data.table(x, key="query.id")
}

#' @export
hivc.seq.rm.drugresistance<- function(char.matrix, dr, verbose=1, rtn.DNAbin=1)
{
	if(verbose)	cat(paste("\nchecking for drug resistance mutations, n=",nrow(dr)))
	tmp	<- rep(0, nrow(dr))
	for(i in seq_len(nrow(dr)))
	{		
		query.yes	<- hivc.seq.find(char.matrix, dr[i,Alignment.nuc.pos], unlist(strsplit(unlist(dr[i,Mutant.NTs]),'')))
		if(length(query.yes))
		{
			if(verbose)	
				cat(paste("\nfound DR mutation",dr[i,DR.name],"NT code",dr[i,Mutant.NTs],"at pos",dr[i,Alignment.nuc.pos],"n=",length(query.yes), ". Replacing NT seq with nnn."))			
			char.matrix[query.yes,	seq.int(dr[i,Alignment.nuc.pos], length.out=3) ]<- matrix("n", nrow=length(query.yes), ncol=3)			
			#print( char.matrix[query.yes, seq.int(dr[i,Alignment.nuc.pos]-3, length.out=9)] ); stop()
		}	
		tmp[i]		<- length(query.yes)
	}
	if(verbose)
		cat(paste("\nremoved DR mutations, n=",sum(tmp), "n per patient=",sum(tmp)/nrow(char.matrix)))
	if(rtn.DNAbin)
		return( as.DNAbin(char.matrix) )
	else
		return( char.matrix )
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
hivc.seq.rmgaps<- function(seq.DNAbin.matrix, rm.only.col.gaps=1, verbose=0)
{
	seq.DNAbin.matrix		<- as.character(seq.DNAbin.matrix)		
	if(!rm.only.col.gaps)
	{	
		if(is.matrix(seq.DNAbin.matrix))
		{
			tmp					<- lapply(seq_len(nrow(seq.DNAbin.matrix)), function(i){	seq.DNAbin.matrix[i, seq.DNAbin.matrix[i,]!="-" & seq.DNAbin.matrix[i,]!="?"]	})
			names(tmp)			<- rownames(seq.DNAbin.matrix)
		}
		else
		{
			tmp					<- lapply(seq_along(seq.DNAbin.matrix), function(i){	seq.DNAbin.matrix[[i]][ seq.DNAbin.matrix[[i]]!="-" & seq.DNAbin.matrix[[i]]!="?"]	})
			names(tmp)			<- names(seq.DNAbin.matrix)
		}		
		seq.DNAbin.matrix	<- tmp
	}
	else
	{		
		nogap				<- which( !apply(seq.DNAbin.matrix,2,function(x) all(x=="-" || x=="?")) )
		if(verbose)	cat(paste("\nremove gaps, n=",ncol(seq.DNAbin.matrix)-length(nogap)))
		seq.DNAbin.matrix	<- seq.DNAbin.matrix[,nogap]	
	}
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
######################################################################################
hivc.phy.get.TP.and.TN.bootstrapvalues<- function(ph, bs.linked.bypatient, ph.mrca=NULL, df.seqinfo=NULL, bs.unlinkedpairs=NULL, bs.unlinked.byspace=NULL, dist.brl=NULL, thresh.brl=NULL, plot.file=NULL, verbose= 1)	
{
	require(phangorn)
	require(RColorBrewer)
	if(is.null(ph.mrca))
		ph.mrca				<- mrca(ph)
	if(any(is.na(ph$node.label)))	stop("Found unexpected NA in ph$node.label")
	#
	# compute bootstrap values between all unlinked dead/seroneg pairs
	#
	if(!is.null(bs.unlinkedpairs))
	{
		if(verbose)	cat(paste("\nCompute bootstrap values between all unlinked dead/seroneg pairs"))
		bs.unlinkedpairs[,dummy:=seq_len(nrow(bs.unlinkedpairs))]
		setkey(bs.unlinkedpairs, mrca)
		tmp					<- unique(bs.unlinkedpairs)	
		tmp					<- tmp[, list( mrca=mrca, ntips= sapply( Descendants(ph, mrca, type="tips"), length) )	]		
		tmp					<- subset(tmp,ntips<=2)
		tmp					<- merge(bs.unlinkedpairs, tmp, by="mrca")
		bs.unlinkedpairs	<- subset(tmp, select=c(tip1, tip2, mrca))
		if(verbose)	cat(paste("\nFound unlinked dead/seroneg pairs, n=",nrow(bs.unlinkedpairs)))
		tmp					<- ph$node.label[ bs.unlinkedpairs[,mrca]-Ntip(ph) ]
		bs.unlinkedpairs[, mrca.bs:=tmp]	
		bs.unlinkedpairs[, bs:=mrca.bs]
	}
	else
		bs.unlinkedpairs	<- NULL
	#
	# compute bootstrap values between all unlinked by space (no restricted to pairs)
	#	
	if(!is.null(bs.unlinked.byspace))
	{
		if(verbose)	cat(paste("\nCompute bootstrap values between all unlinked by space (not restricted to pairs)"))
		setkey(bs.unlinked.byspace, mrca)
		tmp					<- unique(bs.unlinked.byspace)	
		tmp[, dummy:=seq_len(nrow(tmp))]
		tmp[, mrca.bs:=ph$node.label[ tmp[,mrca]-Ntip(ph) ]]	
		tmp					<- tmp[, 	{													
					z			<- Ancestors(ph,mrca)													
					anc.bs		<- ph$node.label[ z-Ntip(ph) ]
					anc.brl		<- dist.brl[ z-Ntip(ph) ]
					z2		<- which( !as.logical(cumsum( as.numeric( anc.brl>=thresh.brl ) )) )
					list(	mrca=mrca, mrca.bs=mrca.bs,
							amrca.bs= ifelse(length(z),max(anc.bs),0), 
							amrca.bs.brl=ifelse(length(z2),max(anc.bs[ z2 ]),0))													
				}, by="dummy"]
		# compute "bs" from "mrca.bs","amrca.bs","amrca.bs.brl"									
		tmp[,bs:=apply( rbind( tmp[, mrca.bs], tmp[,amrca.bs]), 2, max)]		
		tmp[,bs2:=apply( rbind( tmp[, mrca.bs], tmp[,amrca.bs.brl]), 2, max)]
		if(verbose)	cat(paste("\nnumber of bs.unlinked pairs for which brl threshold restricts the bootstrap value, n=",nrow(subset(tmp,bs2>bs)) ))
		bs.unlinked.byspace	<- merge(bs.unlinked.byspace, subset(tmp, select=c(mrca, bs)), by="mrca")
	}
	else
		bs.unlinked.byspace	<- NULL
	#
	# compute bootstrap values between all within patient pairs
	#
	if(verbose)	cat(paste("\nCompute bootstrap values between all TP pairs"))	
	# compute bootstrap values between all within patient sequences by Patient
	tmp					<- ph$node.label[ bs.linked.bypatient[,mrca]-Ntip(ph) ]
	bs.linked.bypatient[,mrca.bs:=tmp]
	bs.linked.bypatient[,dummy:=seq_len(nrow(bs.linked.bypatient))]
	setkey(bs.linked.bypatient,dummy)	
	tmp					<- bs.linked.bypatient[, {													
													tmp			<- Ancestors(ph,mrca)		
													anc.bs		<- ph$node.label[ tmp-Ntip(ph) ]
													anc.brl		<- dist.brl[ tmp-Ntip(ph) ]
													tmp2		<- which( !as.logical(cumsum( as.numeric( anc.brl>=thresh.brl ) )) )
													list(	amrca.bs= ifelse(length(tmp),max(anc.bs),0), 
															amrca.bs.brl=ifelse(length(tmp2),max(anc.bs[ tmp2 ]),0))													
												}, by="dummy"]
	bs.linked.bypatient<- merge(bs.linked.bypatient, tmp, by="dummy")	
	if(verbose)	cat(paste("\nFound bootstrap values between TP pairs, n=",nrow(bs.linked.bypatient)))
	# compute "bs" from "mrca.bs","amrca.bs","amrca.bs.brl"									
	tmp					<- apply( rbind( bs.linked.bypatient[, mrca.bs], bs.linked.bypatient[,amrca.bs]), 2, max)		
	tmp2				<- apply( rbind( bs.linked.bypatient[, mrca.bs], bs.linked.bypatient[,amrca.bs.brl]), 2, max)
	if(verbose)	cat(paste("\nnumber of TP pairs for which brl threshold restricts the bootstrap value, n=",length(which(tmp2>tmp))))
	bs.linked.bypatient[, bs:=tmp]
	# compute BS category
	tmp					<- round(bs.linked.bypatient[, bs]*20,d=0)/2
	bs.linked.bypatient[, bs.cat:=tmp]
	set(bs.linked.bypatient,which(bs.linked.bypatient[, bs.cat==0]),"bs.cat",0.5)
	#
	if(class(bs.linked.bypatient[,tip1])!="character")
		set(bs.linked.bypatient, NULL, "tip1", ph$tip.label[ bs.linked.bypatient[,tip1] ])
	if(class(bs.linked.bypatient[,tip2])!="character")
		set(bs.linked.bypatient, NULL, "tip2", ph$tip.label[ bs.linked.bypatient[,tip2] ])
	#
	# compute time between the two within patient sequences
	#
	if(!is.null(df.seqinfo))
	{
		if(verbose)	cat(paste("\nAdding time difference between seq sampling times"))
		tmp					<- data.table( tip1.PosSeqT= df.seqinfo[bs.linked.bypatient[,tip1],PosSeqT][,PosSeqT], tip2.PosSeqT= df.seqinfo[bs.linked.bypatient[,tip2],PosSeqT][,PosSeqT])
		tmp[, PosSeqT.diff:=tmp[, abs(as.numeric(difftime(tip1.PosSeqT, tip2.PosSeqT, units="days")))/365]]			
		bs.linked.bypatient	<- cbind(bs.linked.bypatient, tmp)
		bs.linked.bypatient	<- subset(bs.linked.bypatient, select=c(Patient, tip1, tip2, bs, bs.cat, PosSeqT.diff))
	}
	else
		bs.linked.bypatient	<- subset(bs.linked.bypatient, select=c(Patient, tip1, tip2, bs, bs.cat))
	#
	# plot boostrap histogram with BS
	#	
	if(!is.null(plot.file) & "PosSeqT.diff"%in%colnames(bs.linked.bypatient))
	{
		breaks.n	<- 20
		if(verbose)	cat(paste("\nPlotting BS distribution with PosSeqT.diff to file",plot.file))
		pdf(width=5,height=6,file=plot.file)
		def.par 	<- par(no.readonly = TRUE)		
		layout( matrix(c(1,1,1,2),ncol=1,nrow=4) )
		par(mar=c(0.5,4,0.5,0.5))		
		cols		<- c(sapply(brewer.pal(4,"Paired"), function(x)  my.fade.col(x, 1))[1:2], "transparent")
		border		<- c("transparent","transparent","black")
		hist(bs.linked.bypatient[,bs], breaks=breaks.n, main='', xlab="", border=border[1], col=cols[1], xaxt='n', freq=1, add=0)		
		hist(subset(bs.linked.bypatient, PosSeqT.diff<1.5)[,bs],breaks=breaks.n, border=border[2], add=1, col=cols[2], freq=1)
		if(!is.null(bs.unlinked.byspace))
			hist(bs.unlinked.byspace[,bs], breaks=breaks.n, main='', xlab="", border=border[3], col=cols[3], xaxt='n', freq=1, add=1)
		legend("topright", bty='n', border=border, legend=c("all pairs of within patient sequences","pairs of within patient sequences\n with difference in time of sampling < 18 mo","all pairs with geographically distant seq"), fill=cols)
		
		tmp			<- bs.linked.bypatient[, list(prop.recent= length(which(PosSeqT.diff<1.5))/length(PosSeqT.diff) ),by=bs.cat]
		setkey(tmp, bs.cat)
		
		par(mar=c(5,4,0.5,0.5))
		plot(1,1,type='n',xlim=range(bs.linked.bypatient[,bs]),ylim=c(0,1),xlab="bootstrap",ylab='prop',bty='n')
		polygon(c(range(bs.linked.bypatient[,bs]),rev(range(bs.linked.bypatient[,bs]))), c(0,0,1,1), col=cols[1], border=NA)
		sapply(seq_len(nrow(tmp)),function(i)
				{										 
					polygon(c(tmp[i,bs.cat]/10-0.05,tmp[i,bs.cat]/10,tmp[i,bs.cat]/10,tmp[i,bs.cat]/10-0.05),c(0,0,rep(tmp[i,prop.recent],2)), border=NA, col=cols[2])
				})
		par(def.par)
		dev.off()
	}
	if(!is.null(plot.file) & !"PosSeqT.diff"%in%colnames(bs.linked.bypatient))
	{
		breaks.n	<- 20
		if(verbose)	cat(paste("\nPlotting BS distribution without PosSeqT.diff to file",plot.file))
		pdf(width=5,height=6,file=plot.file)
		par(mar=c(5,4,0.5,0.5))		
		cols		<- c(sapply(brewer.pal(4,"Paired"), function(x)  my.fade.col(x, 1))[1], "transparent")
		border		<- c("transparent","transparent","black")
		hist(bs.linked.bypatient[,bs], breaks=breaks.n, main='', xlab="bootstrap", border=border[1], col=cols[1], freq=1, add=0)
		if(is.null(bs.unlinked.byspace))
			legend	<- c("all pairs of within patient sequences")
		else
		{
			hist(bs.unlinked.byspace[,bs], breaks=breaks.n, main='', xlab="", border=border[2], col=cols[2], xaxt='n', freq=1, add=1)
			legend	<- c("all pairs of within patient sequences","all pairs with geographically distant seq")
		}
		legend("topright", bty='n', border=border[seq_along(legend)], legend=legend, fill=cols[seq_along(legend)])
		dev.off()
	}
	
	list(bs.linked.bypatient=bs.linked.bypatient, bs.unlinkedpairs=bs.unlinkedpairs, bs.unlinked.byspace=bs.unlinked.byspace)
}
######################################################################################
hivc.phy.get.TP.and.TN<- function(ph, df.all, use.seroneg.as.is= 0, verbose= 1)
{
	#
	# exctract geographically distant seqs that tips in 'ph' and that are assumed to be truly unlinked to NL seqs
	#
	tmp					<- ph$tip.label[ substr(ph$tip.label,1,2)=="TN" ]
	unlinked.byspace	<- data.table(FASTASampleCode=tmp, key="FASTASampleCode")		
	#
	# extract list of truly linked sample codes
	#
	tmp					<- merge(df.all, data.table(FASTASampleCode=ph$tip.label), by="FASTASampleCode")	
	tmp					<- subset(tmp[, length(FASTASampleCode)>1, by=Patient], V1==T, select=Patient)		#patients for which more than one sequence is available
	linked.bypatient	<- merge(tmp, merge(df.all, data.table(FASTASampleCode=ph$tip.label), by="FASTASampleCode"), by="Patient")
	setkey(linked.bypatient, "FASTASampleCode")
	linked.bypatient	<- subset(linked.bypatient, select=c(Patient, FASTASampleCode, PosSeqT))	
	#
	# extract list of truly unlinked sample codes by temporal separation
	#
	df.dead					<- subset(df.all, !is.na(DateDied))
	if(verbose) cat(paste("\nnumber of dead HIV+ individuals",nrow(df.dead)))
	df.dead					<- merge(df.dead, data.table(FASTASampleCode=ph$tip.label), by="FASTASampleCode")
	if(verbose) cat(paste("\nnumber of dead HIV+ individuals in 'ph'",nrow(df.dead)))
	setkey(df.dead, DateDied)	
	#extract seroconverters
	
	tmp							<- df.dead[1,DateDied]
	df.serocon.acc				<- subset(df.all, NegT_Acc=="Yes" & NegT>=tmp)
	#add seroconverters with inaccurate info
	df.serocon.nacc				<- subset(df.all, NegT_Acc=="No" & !is.na(NegT) & NegT>=tmp )
	if(!use.seroneg.as.is)
	{
		df.serocon.nacc.dy			<- subset(df.serocon.nacc, as.POSIXlt(NegT)$mday==15, )						#for inaccurate days, we (conservatively) assume the patient was only seronegative at the start of the month
		if(verbose)	cat(paste("\nnumber of seroconverters with innacurate days, n=", nrow(df.serocon.nacc.dy)))
		if(nrow(df.serocon.nacc.dy))
		{
			tmp							<- as.POSIXlt(df.serocon.nacc.dy[,NegT] )
			tmp$mday					<- 1
			df.serocon.nacc.dy[,NegT:=as.Date(tmp)]		
		}
		df.serocon.nacc.mody		<- subset(df.serocon.nacc, as.POSIXlt(NegT)$mon==6 & as.POSIXlt(NegT)$mday==1, )		#for inaccurate months and days, we (conservatively) assume the patient was only seronegative at the start of the year
		if(verbose)	cat(paste("\nnumber of seroconverters with innacurate months & days, n=", nrow(df.serocon.nacc.mody)))
		if(nrow(df.serocon.nacc.mody))
		{
			tmp							<- as.POSIXlt(df.serocon.nacc.mody[,NegT] )
			tmp$mon						<- 0
			df.serocon.nacc.mody[,NegT:=as.Date(tmp)]
		}
		df.serocon.nacc					<- rbind(df.serocon.nacc.dy, df.serocon.nacc.mody)
	}
	df.serocon					<- rbind(df.serocon.acc, df.serocon.nacc)
	set(df.serocon, NULL, "FASTASampleCode", as.character(df.serocon[,FASTASampleCode]))
	if(verbose) cat(paste("\nnumber of seroconverters with at least 1 preceeding dead HIV+",nrow(df.serocon)))
	df.serocon					<- merge(df.serocon, data.table(FASTASampleCode=ph$tip.label), by="FASTASampleCode")
	if(verbose) cat(paste("\nnumber of seroconverters in 'ph' with at least 1 preceeding dead HIV+",nrow(df.serocon)))
	#for each seq of seroconverter, extract HIV+ seqs that are dead before seroconversion	
	if( any(as.logical(df.serocon[,is.na(NegT)])) )	warning("Found accurate seroconverters with missing NegT")		
	setkey(df.serocon,NegT)
	unlinked.bytime				<- lapply(seq_len(nrow(df.serocon)), function(i)
				{
					tmp				<- df.serocon[i,NegT]
					tmp				<- subset(df.dead, DateDied<=tmp,select=c(Patient,FASTASampleCode,DateDied))
					tmp2			<- rep(df.serocon[i,FASTASampleCode],nrow(tmp))
					tmp[,query.FASTASampleCode:= tmp2]
					setkey(tmp,"FASTASampleCode")
					tmp					
				})
	names(unlinked.bytime)		<- df.serocon[,FASTASampleCode]											
	#
	#plot number of unlinked HIV+ by date (this date is when someone else is still seroneg)
	#
	if(0)
	{
			#
			y		<- sapply(unlinked.bytime, function(x) nrow(x) )
			y2		<- nrow(unlinked.byspace)
			x		<- df.serocon[,NegT]
			xlim	<- range(x)
			ylim	<- c(0, y2+max(y))
			tmp		<- as.POSIXlt(xlim[1])
			tmp$mday<- 1
			tmp$mon	<- 1
			xlim[1]	<- as.Date(tmp)
			plot(1,1,type='n',bty='n',xlab="time of seroconversion",ylab="number unlinked",xaxt='n',xlim=xlim,ylim=ylim)
			axis.Date(1,Year,at=seq(xlim[1], xlim[2],by="12 months"))			
			polygon(c(x[1],x[length(x)],x[length(x)],x[1]), c(y2,y2,0,0), border=NA, col="blue" )
			polygon(c(x,c(x[length(x)],x[1])), c(y+y2,y2,y2), border=NA, col="grey60" )
			legend("topleft",bty='n',fill=c("blue","grey60"),legend=c("by space","by time"),border=NA)						
	}	
	#
	# get TN for seroneg: convert unlinked.bytime from SampleCode to ph node index and add truly unlinked.byspace
	#
	df.tips							<- data.table(Node=seq_len(Ntip(ph)), FASTASampleCode=ph$tip.label, key="FASTASampleCode" )
	unlinked.byspace				<- merge(unlinked.byspace, df.tips, by="FASTASampleCode")
	setkey(unlinked.byspace,"Node")
	ph.unlinked.bytime				<- lapply(unlinked.bytime, function(x)
				{
					tmp	<- merge(x, df.tips, all.x=1)
					tmp2<- tmp[1,query.FASTASampleCode]
					tmp	<- rbind( subset(tmp, select=c(FASTASampleCode,Node)), unlinked.byspace )
					tmp2<- rep(tmp2, nrow(tmp))
					tmp[,query.FASTASampleCode:= tmp2]
					setkey(tmp, Node)
					tmp
				})	
	names(ph.unlinked.bytime)		<- names(unlinked.bytime)
	#
	# get TN for seropos: use only unlinked.byspace
	#
	ph.seroneg						<- merge( data.table(FASTASampleCode=names(ph.unlinked.bytime)), df.tips, by="FASTASampleCode", all.x=1)		
	ph.seropos						<- subset( df.tips[,c(2,1), with=F], !FASTASampleCode%in%names(ph.unlinked.bytime) 
														& 	substr(FASTASampleCode,1,2)!="TN"
														&	substr(FASTASampleCode,1,8)!="PROT+P51"
														)			
	ph.unlinked.byspace				<- lapply(seq_len(nrow(ph.seropos)), function(i)
			{
				tmp	<- unlinked.byspace
				tmp2<- rep(ph.seropos[i,FASTASampleCode], nrow(tmp))
				tmp[,query.FASTASampleCode:= tmp2]
				setkey(tmp, Node)
				tmp
			})
	names(ph.unlinked.byspace)		<- ph.seropos[,FASTASampleCode]
	#
	# put all unlinked ATHENA seqs together and sort by Node
	#
	ph.unlinked						<- c(ph.unlinked.byspace, ph.unlinked.bytime)
	names(ph.unlinked)				<- c(names(ph.unlinked.byspace),names(ph.unlinked.bytime))
	ph.unlinked.info				<- rbind(ph.seropos,ph.seroneg)		
	setkey(ph.unlinked.info, "Node")
	ph.unlinked						<- lapply(ph.unlinked.info[,FASTASampleCode], function(i){		ph.unlinked[[i]]	})
	names(ph.unlinked)				<- ph.unlinked.info[,FASTASampleCode]
	tmp								<- seq_len(nrow(ph.unlinked.info))
	ph.unlinked.info[,NodeIdx:= tmp]
	ph.unlinked.info				<- merge(ph.unlinked.info, df.all, all.x=1, by="FASTASampleCode")
	setkey(ph.unlinked.info, "Node")		
	#
	#	get data table of all linked ATHENA seqs
	#
	ph.linked						<- merge(linked.bypatient, df.tips)
	ph.linked						<- subset(ph.linked, select=c(FASTASampleCode, Patient, Node))
	ph.linked						<- merge(ph.linked, ph.linked[, length(FASTASampleCode), by=Patient], all.x=1, by="Patient")
	setnames(ph.linked, "V1","nTP")
	setkey(ph.linked, "Node")
	ans								<- list(	unlinked.byspace=unlinked.byspace, unlinked.bytime=unlinked.bytime, linked.bypatient=linked.bypatient,
												ph.linked=ph.linked, ph.unlinked.info=ph.unlinked.info, ph.unlinked=ph.unlinked)
	ans									
}
######################################################################################
hivc.get.hpcsys<- function()
{
	tmp<- system('domainname',intern=T)
	if(!nchar(tmp))	tmp<- "debug"
	tmp
}

