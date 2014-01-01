#' @export
HIVC.COUNTRY.TABLE<- data.table(matrix(c(	"AG","ANTIGUA AND BARBUDA",		"AL","ALBANIA",		"AN", "ANTILLES", "AO","ANGOLA",	"AR","ARGENTINA",	"AT","AUSTRIA",		"AU","AUSTRALIA",				
						"BB","BARBADOS",	"BE","BELGIUM",	"BG","BULGARIA",	"BR","BRAZIL",	"CA","CANADA",	"CH","SWITZERLAND",		"CI","COTE D'IVOIRE",	"CL","CHILE",	
						"CN","CHINA",	"CO","COLOMBIA",	"CU","CUBA",	"CW", "CURACAO", "CV","CAPE VERDE",	"CY","CYPRUS",	"CZ","CZECH REPUBLIC",	"DE","GERMANY",	"DK","DENMARK",	
						"DO","DOMINICAN REPUBLIC",	"DZ","ALGERIA",		"EC","ECUADOR",	"ER", "ERITREA", "ES","SPAIN",	"FI","FINLAND",	"FR","FRANCE",	"GA","GABON",	"GB","UNITED KINGDOM",	
						"GD","GRENADA",	"GE","GIORGIA",	"GH","GHANA",	"GR","GREECE",	"GW","GUINEA-BISSAU",	"GY","GUYANA",	"HN","HONDURAS", "HT","HAITI",	
						"HU","HUNGARY",	"ID", "INDONESIA", "IL","ISRAEL",	"IN","INDIA",	"IT","ITALY",	"JM","JAMAICA",	"JP","JAPAN",	"KE","KENYA",	"KR","SOUTH KOREA",	"LB","LEBANON",				
						"LU","LUXEMBOURG",	"LK", "SRILANKA", "LV","LATVIA",	"MA","MOROCCO",	"ME","MONTENEGRO",	"MG","MADAGASCAR",	"ML","MALI",	"MN","MONGOLIA",	"MX","MEXICO",	"MY","MALAYSIA",
						"NG","NIGERIA",	"NL","NETHERLANDS",	"NO","NORWAY",	"PA","PANAMA",	"PE","PERU",	"PH","PHILIPPINES",	"PL","POLAND",	"PR","PUERTO RICO",	"PT","PORTUGAL",			
						"PY","PARAGUAY",	"RO","ROMANIA",	"RS","SERBIA",	"RU","RUSSIAN FEDERATION",	"SC","SEYCHELLES",	"SD","SUDAN",	"SE","SWEDEN",	"SG","SINGAPORE",	
						"SI","SLOVENIA",	"SK","SLOVAKIA",	"SN","SENEGAL",	"SR","SURINAME",	"SV","EL SALVADOR",	"TH","THAILAND",	"TT","TRINIDAD AND TOBAGO",	"TW","TAIWAN",			
						"UA","UKRAINE",	"UG","UGANDA",	"US","UNITED STATES",	"UY","URUGUAY",	"VE","VENEZUELA",	"VN","VIETNAM",	"YE","YEMEN","ZA","SOUTH AFRICA"),ncol=2,byrow=1,dimnames=list(c(),c("key","country"))))

######################################################################################
project.hivc.check.DateRes.after.HIVPosTest<- function(dir.name= DATA, verbose=1)
{
	require(data.table)
	require(ggplot2)
	require(plyr)
	
	NL.HIV.phases<- as.Date( c("1980-01-01","1984-01-01","1996-01-01","2000-01-01","2004-01-01") )
	NL.possibly.Acute<- c(1,2)
	
	file	<- paste(dir.name,"tmp/ATHENA_2013_03_Sequences.R",sep='/')
	load(file)
	df.seq	<- df	
	file	<- paste(dir.name,"tmp/ATHENA_2013_03_Patients.R",sep='/')
	load(file)
	df.pat	<- df
	
	out		<- lapply( seq_along(df.seq),function(i)
			{
				cat(paste("\nprocess", names(df.seq)[i]))
				#no missing DateRes
				nok.idx<- which( is.na(df.seq[[i]][,"DateRes"]) )
				if(verbose)		cat(paste("\nentries in Sequences", nrow(df.seq[[i]])))
				if(verbose)		cat(paste("\nentries with missing DateRes", length(nok.idx)))
				seq.ok.idx<- which( !is.na(df.seq[[i]][,"DateRes"]) )
				df.seq[[i]]<- df.seq[[i]][seq.ok.idx,]								
				
				#extract min DateRes per patient
				df<- data.table( df.seq[[i]][,c("Patient","DateRes")], key="Patient" )
				df.minDateRes<- df[,min(DateRes),by=Patient]
				setnames(df.minDateRes,"V1","DateRes")
				#print(df.minDateRes)				
				if(verbose)		cat(paste("\npatients with DateRes", nrow(df.minDateRes)))
				
				#no missing MyDatePos1
				if(verbose)		cat(paste("\nentries in Patients", nrow(df.pat)))
				nok.idx<- which( is.na(df.pat[,"MyDatePos1"]) )				
				if(verbose)		cat(paste("\nentries with missing MyDatePos1", length(nok.idx)))				
				ok.idx<- which( !is.na(df.pat[,"MyDatePos1"]) )
				df.pos<- data.table( df.pat[ok.idx,c("Patient","MyDatePos1","MyDatePos1_Acc","AcuteInfection")], key="Patient" )				
				if(verbose)		cat(paste("\npatients with MyDatePos1", nrow(df.pos)))
				
				print(df.pos)
				#exclude Acute for checking HIVPos against DateRes
				df.pos<- subset(df.pos, !AcuteInfection%in%NL.possibly.Acute)
				if(verbose)		cat(paste("\npatients with MyDatePos1 and not known to acute or possibly acute", nrow(df.pos)))				
				
				#merge
				df<- merge(df.minDateRes,df.pos)
				if(verbose)		cat(paste("\npatients with HIVPos and DateRes", nrow(df)))
				#print(df)
				
				df.cmpHIVPos<- df[,difftime(MyDatePos1,DateRes,units="days")]
				df.cmpHIVPos<- as.numeric(df.cmpHIVPos)/12
				tmp<- which(df.cmpHIVPos>0)
				if(verbose)		cat(paste("\npatients with HIVPos>DateRes, n=", length(tmp)))	
				patient.laterHIVPos<- df[tmp,Patient]
				
				#in cases where MyDatePos1>DateRes, the guess of MyDatePos1 if MyDatePos1_Acc=0 may not be appropriate, improve this guess.								
				tmp<- which( df[,df.cmpHIVPos>0 & MyDatePos1_Acc==0 & as.POSIXlt(MyDatePos1)$mday==15] )
				if(verbose)		cat(paste("\npatients with unclear MyDatePos1 and day=15 and HIVPos>DateRes, n=", length(tmp)))				
				patient.betterHIVPos<- subset(df, df.cmpHIVPos>0 & MyDatePos1_Acc==0 & as.POSIXlt(MyDatePos1)$mday==15 & as.POSIXlt(MyDatePos1)$mon==as.POSIXlt(DateRes)$mon )
				if(verbose)		cat(paste("\npatients above for which HIVPos is in same month as DateRes (ie HIVPos can be fixed), n=", nrow(patient.betterHIVPos)))
				tmp<- subset( df, df.cmpHIVPos>0 & MyDatePos1_Acc==0 &  as.POSIXlt(MyDatePos1)$mday==1 & as.POSIXlt(MyDatePos1)$mon==6 )
				if(verbose)		cat(paste("\npatients with unclear MyDatePos1 and day=1, month=7 and HIVPos>DateRes, n=", nrow(tmp)))
				tmp<- subset( df, df.cmpHIVPos>0 & MyDatePos1_Acc==0 &  as.POSIXlt(MyDatePos1)$mday==1 & as.POSIXlt(MyDatePos1)$mon==6 & as.POSIXlt(MyDatePos1)$year==as.POSIXlt(DateRes)$year )
				if(verbose)		cat(paste("\npatients above for which HIVPos is in same year as DateRes (ie HIVPos can be fixed), n=", nrow(tmp)))				
				patient.betterHIVPos<- rbind(patient.betterHIVPos, tmp)
				setkey(patient.betterHIVPos, Patient)						#important to sort for next step, which assumes sorted
				tmp<- df.pos[,MyDatePos1]
				tmp[ which( df.pos[,Patient%in%patient.betterHIVPos$Patient]) ]<- patient.betterHIVPos$DateRes
				df.pos[,orDatePos1:=tmp]
				
				#re-merge
				df<- merge(df.minDateRes,df.pos)
				df.cmpHIVPos<- df[,difftime(orDatePos1,DateRes,units="days")]
				df.cmpHIVPos<- as.numeric(df.cmpHIVPos)/12
				tmp<- which(df.cmpHIVPos>0)
				if(verbose)		cat(paste("\npatients with HIVPos>DateRes after FIX, n=", length(tmp)))
				hist(df.cmpHIVPos[tmp], breaks=100)
				print( df[tmp,] )
				
				#tmp<- which(df.cmpHIVPos>1)
				#if(verbose)		cat(paste("\npatients with HIVPos>DateRes after FIX and 1 month grace, n=", length(tmp)))
				#print( df[tmp,] )
				
				patient.laterHIVPos<- df[tmp,Patient]				
				#if(verbose)		cat(paste("\npatients with HIVPos>DateRes, Patient", paste(df[tmp,Patient],collapse=', ') ))
				#print(df[tmp,])
				#hist(df.cmpHIVPos, breaks=50)
				
				if(0)
				{
					df.earlierHIVPos<- -df.cmpHIVPos[ which(df.cmpHIVPos<=0) ]
					df.earlierHIVPos[ df.earlierHIVPos < 1 ]<- 1		#if sequenced within 1mo after diagnosis, ignore
					df.earlierHIVPos<- sort(df.earlierHIVPos)
					col<- "grey70"
					x<- df.earlierHIVPos
					y<- seq_along(df.earlierHIVPos) / length(df.earlierHIVPos)				
					plot(1,1,type='n',bty='n',xlab="time [months]", ylab="%sequenced after diagnosis", log='x', xlim=range(x),ylim=range(y))				
					apply( matrix(seq(1,max(x),by=6),2) ,2,function(x)
							{
								polygon(c(x,rev(x)),c(0,0,1,1),border=NA,col=col)
							})
					lines(x,y,type='s')
				}
				if(1)
				{											
					tmp		<- numeric(nrow(df))
					tmp[ which( df[,orDatePos1]<NL.HIV.phases[2] ) ]<- 1
					tmp[ which( df[,orDatePos1]>=NL.HIV.phases[2] & df[,orDatePos1]<NL.HIV.phases[3]) ]<- 2
					tmp[ which( df[,orDatePos1]>=NL.HIV.phases[3] & df[,orDatePos1]<NL.HIV.phases[4]) ]<- 3
					tmp[ which( df[,orDatePos1]>=NL.HIV.phases[4] & df[,orDatePos1]<NL.HIV.phases[5]) ]<- 4
					tmp[ which( df[,orDatePos1]>NL.HIV.phases[5] ) ]<- 5					
					tmp		<- data.frame( HIVPos=df[,orDatePos1], TimeToSeq= -df.cmpHIVPos, HIVPhase=tmp )
					tmp		<- tmp[ tmp$TimeToSeq>=0, ]
					tmp$TimeToSeq[tmp$TimeToSeq<1]<- 1
					#cumsum(phase) / nrow
					#print( cumsum( as.numeric( tmp$HIVPhase==1 ) ) )
					#quit("no")
					tmp		<- cbind( tmp[ with(tmp, order(TimeToSeq)), ], seq_len(nrow(tmp)) )
					counts	<- t(ddply(tmp, .(HIVPhase), function(x)
							{								
								x<- x[with(x,order(TimeToSeq)),]
								z<- numeric(nrow(tmp))
								for(j in seq_len(nrow(x)))
								{
									z[x[j,4]:length(z)]<- 1+z[x[j,4]:length(z)]
								}
								z
							}))
					counts<- counts[-1,]		
					tmp<- as.data.frame( cbind(tmp$TimeToSeq,counts) )
					colnames(tmp)<- c("TimeToSeq",paste("ph",1:5,sep=''))
					
					cols	<- rainbow(5)
					z		<- numeric(nrow(tmp))
					xlim	<- c(1,max(tmp$TimeToSeq))
					plot(1,1,type='n',xlab="Time from Diag To Seq [months]",ylab="number seq",xlim=xlim,ylim=c(0,nrow(df)), log='x')					
					for(j in 2:ncol(tmp))
							{
								polygon( c(tmp$TimeToSeq,rev(tmp$TimeToSeq)), c(z,rev(z+tmp[,j])), border=NA, col=cols[j-1] )
								z<- z+tmp[,j]
							}
					legend("topleft",fill=cols,legend=colnames(tmp)[2:6],bty='n', border=NA)		
					#need different data frame see  http://stackoverflow.com/questions/5030389/getting-a-stacked-area-plot-in-r 
					#p<- 	ggplot(data=tmp,aes(x=TimeToSeq,y=val, group=HIVPhase, colour= HIVPhase))
					#p<- 	p+geom_line(aes(col=HIVPhase))
					#p<-		p+geom_area(position = "fill")	
					#print(p)					
					quit("no")
				}												
				list(patient.laterHIVPos=patient.laterHIVPos)
			})
	patient.laterHIVPos	<- unique( unlist( lapply(out,function(x) x[["patient.laterHIVPos"]]) ) )
	if(verbose)		cat(paste("\npatients with HIVPos>DateRes RT or PROT, n=", length(patient.laterHIVPos)))
	if(verbose)		cat(paste("\npatients with HIVPos>DateRes RT or PROT, Patient", paste(patient.laterHIVPos,collapse=', ') ))		
}
######################################################################################
project.hivc.check.DateRes.after.T0<- function(dir.name= DATA, verbose=1)
{
	require(data.table)
	
	file<- paste(dir.name,"tmp/ATHENA_2013_03_Sequences.R",sep='/')
	load(file)
	df.seq<- df

	file<- paste(dir.name,"tmp/ATHENA_2013_03_Regimens.R",sep='/')
	load(file)
	df.reg<- df
	
	lapply( seq_along(df.seq),function(i)
			{
				cat(paste("\nprocess", names(df.seq)[i]))
				#no missing DateRes
				nok.idx<- which( is.na(df.seq[[i]][,"DateRes"]) )
				if(verbose)		cat(paste("\nentries in Sequences", nrow(df.seq[[i]])))
				if(verbose)		cat(paste("\nentries with missing DateRes", length(nok.idx)))
				seq.ok.idx<- which( !is.na(df.seq[[i]][,"DateRes"]) )
				df.seq[[i]]<- df.seq[[i]][seq.ok.idx,]
				#print( range(df.seq[[i]][,"DateRes"]) )
				
				
				df<- data.table( df.seq[[i]][,c("Patient","DateRes")], key="Patient" )
				df.minDateRes<- df[,min(DateRes),by=Patient]
				setnames(df.minDateRes,"V1","DateRes")
				print(df.minDateRes)
				#print( range(df.minDateRes[,DateRes]) )
				if(verbose)		cat(paste("\npatients with DateRes", nrow(df.minDateRes)))
				
				#no missing T0				
				nok.idx<- which( is.na(df.reg[,"T0"]) )
				if(verbose)		cat(paste("\nentries in Regimens", nrow(df.reg)))
				if(verbose)		cat(paste("\nentries with missing T0", length(nok.idx)))
				reg.ok.idx<- which( !is.na(df.reg[,"T0"]) )
				df.reg<- df.reg[reg.ok.idx,]
				
				#extract min T0 per patient
				df<- data.table( df.reg[,c("Patient","T0")], key="Patient" )
				df.minT0<- df[,min(T0),by=Patient]
				setnames(df.minT0,"V1","T0")
				if(verbose)		cat(paste("\npatients with T0", nrow(df.minT0)))
				print(df.minT0)
				
				#merge
				df<- merge(df.minDateRes,df.minT0)
				if(verbose)		cat(paste("\npatients with T0 and DateRes", nrow(df)))
				
				#print(df[6258:6262,])
				df.laterT0<- df[,difftime(T0,DateRes,units="days")]
				df.laterT0<- as.numeric(df.laterT0)/365.25
				tmp<- which(df.laterT0>0)
				if(verbose)		cat(paste("\npatients with T0>DateRes, n=", length(tmp)))				
				#if(verbose)		cat(paste("\npatients with T0>DateRes, Patient", paste(df[tmp,Patient],collapse=', ') ))
				#print(df[tmp,])
				hist(df.laterT0, breaks=50)
				#print(range(df.laterT0))				
			})
}
######################################################################################
project.hivc.check<- function()
{
	if(1) project.hivc.check.DateRes.after.HIVPosTest()
	if(0) project.hivc.check.DateRes.after.T0()	
}
######################################################################################
project.hivc.get.geneticdist.from.sdc<- function(dir.name= DATA)
{	
	tmp<- hivc.clu.geneticdist.cutoff(dir.name=dir.name, plot=1, verbose=1, level.retain.unlinked=0.05)
	print(tmp)
}
######################################################################################
project.hivc.seq.dataset.mDR.mRC.mSH.pLANL<- function()	
{
	require(data.table)
	require(phangorn)
	
	file		<- paste(DATA,"/tmp/","ATHENA_2013_03_NoRCDRAll+LANL_Sequences","_",gsub('/',':',"Fri_Nov_01_16/07/23_2013"),".R",sep='')
	if(verbose) cat(paste("\nread",file))
	tmp			<- load(file)
	seq.len		<- hivc.seq.length(seq.PROT.RT)		
	hist(seq.len, breaks=20)
	seq.amb		<- hivc.seq.proportion.ambiguous(seq.PROT.RT)
	hist(seq.amb, breaks=40)
	
	
	#
	# get clusters for No Recombination + No Drug resistance mutations, single linkage criterion		
	#						
	infilecov		<- "ATHENA_2013_03_AllSeqPatientCovariates"
	infile			<- "ATHENA_2013_03_NoRCDRAll+LANL_Sequences_examlbs500"			
	insignat		<- "Fri_Nov_01_16/07/23_2013"		
	argv			<<- hivc.cmd.preclustering(paste(DATA,"/tmp",sep=''), infile, insignat, paste(DATA,"/derived",sep=''), infilecov, resume=1)				 
	argv			<<- unlist(strsplit(argv,' '))
	nrc.clu.pre		<- hivc.prog.get.clustering.precompute()
	#
	ph				<- nrc.clu.pre$ph
	ph.node.bs		<- as.numeric(ph$node.label)
	
	seq.stat		<- data.table(FASTASampleCode=names(seq.len), tip:=match(seq.stat[,FASTASampleCode],ph$tip.label), length=seq.len, pamb=seq.amb)
	seq.stat		<- seq.stat[,	{
				tmp		<- Ancestors(ph, tip, type='all') - Ntip(ph)
				list(bs.mx=max(ph.node.bs[tmp]), length=length, pamb=pamb, tip=tip)			
			},by=FASTASampleCode]		
	seq.short		<- subset(seq.stat, length<400)
	seq.long		<- subset(seq.stat, length>=400)
	seq.short.nfrgn	<- subset(seq.short,substr(seq.short[,FASTASampleCode],1,2)!="TN")		
	seq.closefrgn	<- subset(seq.stat, substr(seq.stat[,FASTASampleCode],1,8)=="PROT+P51") 
	seq.tn			<- subset(seq.stat, substr(seq.stat[,FASTASampleCode],1,2)=="TN")	
	
	hist(seq.short[,bs.mx], breaks=20)			
	hist(seq.long[,bs.mx], breaks=20)
	hist(seq.short.nfrgn[,bs.mx], breaks=20)
	hist(seq.closefrgn[,bs.mx], breaks=20)
	#	-> longer sequences have larger bootstrap
	#	-> short ATHENA sequences have small bootstap
	#	-> almost all short sequences (500 / 548 ) are TN.	so TNs might not cluster simply because they are short. restrict TNs to foreign sequences > 600
	#	-> all PROT+P51 sequences are long, and do not cluster as well as the NL sequences
	seq.athena.exclude	<- subset(seq.short.nfrgn, bs.mx<0.6 )[, FASTASampleCode]		
	#			
	indir		<- paste(DATA,"tmp",sep='/')		
	infile		<- "ATHENA_2013_03_NoRCDRAll+LANL_Sequences"	
	insignat	<- "Fri_Nov_01_16/07/23_2013"
	outfile		<- "ATHENA_2013_03_-DR-RC-SH+LANL_Sequences"	
	outsignat	<- "Wed_Dec_18_11/37/00_2013"		
	
	file		<- paste(indir,'/',infile,'_',gsub('/',':',insignat),".R",sep='')
	load(file)
	#	exclude short ATHENA sequences that have BS<0.6 in pilot run AND all TN sequences because they are <400 nt
	seq.keep	<- setdiff(rownames(seq.PROT.RT), seq.athena.exclude)
	seq.keep	<- seq.keep[ substr(seq.keep,1,2)!="TN" ]				
	seq.PROT.RT	<- seq.PROT.RT[seq.keep,]
	if(verbose)	cat(paste("\nnumber of long sequences, n=",nrow(seq.PROT.RT)))
	file		<- paste(indir,'/',outfile,'_',gsub('/',':',outsignat),".R",sep='')
	if(verbose)	cat(paste("\nsave new file to",file))
	save(seq.PROT.RT, file=file)
}				
######################################################################################
project.hivc.collectpatientdata<- function(dir.name= DATA, verbose=1, resume=0)
{	
	require(data.table)		
	
	#input files generated with "project.hivc.Excel2dataframe"
	resume			<- 0
	verbose			<- 1
	file.seq		<- paste(dir.name,"derived/ATHENA_2013_03_Sequences.R",sep='/')
	file.patient	<- paste(dir.name,"derived/ATHENA_2013_03_Patients.R",sep='/')
	file.viro		<- paste(dir.name,"derived/ATHENA_2013_03_Viro.R",sep='/')
	file.immu		<- paste(dir.name,"derived/ATHENA_2013_03_Immu.R",sep='/')
	file.treatment	<- paste(dir.name,"derived/ATHENA_2013_03_Regimens.R",sep='/')
	
	#compute file AllSeqPatientCovariates
	file.out.name	<- "ATHENA_2013_03_AllSeqPatientCovariates"
	file.out		<- paste(dir.name,"/derived/",file.out.name,".R",sep='')
	if(resume)												#//load if there is R Master data.table
	{
		options(show.error.messages = FALSE)		
		readAttempt<-try(suppressWarnings(load(file.out)))
		if(!inherits(readAttempt, "try-error"))	cat(paste("\nresume file",file.out))
		if(!inherits(readAttempt, "try-error"))	str(df.all)
		options(show.error.messages = TRUE)		
	}
	if(!resume || inherits(readAttempt, "try-error"))		#else generate R Master data.table
	{
		#get data.table of seqs		
		load(file.seq)			
		df		<- lapply( df,function(x)	data.table(x[,c("Patient","DateRes","SampleCode")], key="SampleCode")	)
		df		<- merge(df[[1]],df[[2]],all.x=1,all.y=1)
		setkey(df, Patient.x)
		tmp		<- which(df[,is.na(Patient.x)])
		set(df, tmp, "Patient.x", df[tmp, Patient.y])
		set(df, tmp, "DateRes.x", df[tmp, DateRes.y])
		set(df, NULL, "SampleCode", gsub(' ','', df[, SampleCode]))
		setnames(df, c("SampleCode","Patient.x","DateRes.x"), c("FASTASampleCode","Patient","PosSeqT"))
		set(df, NULL, "FASTASampleCode", factor(df[,FASTASampleCode]))
		set(df, NULL, "Patient", factor(df[,Patient]))		
		df.all	<- subset(df, select=c(FASTASampleCode,Patient,PosSeqT))
		if(verbose)		cat(paste("\nnumber of sequences found, n=", nrow(df.all)))
		#
		# add Patient data
		#
		if(verbose)		cat(paste("\nadding patient data"))
		load(file.patient)		
		df.all	<- merge(df.all, df, all.x=1, by="Patient")
		setnames(df.all, c("MyDateNeg1","MyDatePos1"), c("NegT","PosT"))
		#	reset PosT_Acc=='No' conservatively to end of year / month
		df.all	<- hivc.db.reset.inaccuratePosT(df.all, nacc.dy.dy= 30, nacc.mody.mo= 11, nacc.mody.dy= 31, verbose=1)
		#	reset NegT_Acc=='No' conservatively to start of year / month
		df.all	<- hivc.db.reset.inaccurateNegT(df.all)
		#	check for clashes in NegT and PosT
		tmp		<- subset(df.all, !is.na(PosT) & !is.na(NegT) & PosT<NegT)
		if(verbose)		cat(paste("\nchecking for clashes in NegT and PosT"))
		if(verbose)		print(tmp)
		if(verbose)		cat(paste("\nmanually resolving clashes -- M41654 has wrong PosT since PosRNA much later -- M12982 has wrong NegT as PosRNA before that time"))
		#	resolving M41654
		tmp		<- which(df.all[,FASTASampleCode=="M4165429052012"])
		set(df.all, tmp, "PosT", df.all[tmp, PosSeqT]) 		
		# 	add preliminary AnyPos_T1	-- since PosT conservative we are good to set irrespective of PosT_Acc	
		df.all[, AnyPos_T1:=PosSeqT]
		tmp		<- which( df.all[, !is.na(PosT) & ( is.na(PosSeqT) | PosT<PosSeqT ) ] )
		if(verbose)		cat(paste("\nbuilding prel AnyPos_T1. Number of seq with !is.na(PosT) & ( is.na(PosSeqT) | PosT<PosSeqT ), n=",length(tmp)))
		set(df.all, tmp, "AnyPos_T1", df.all[tmp, PosT])		
		if(nrow(subset(df.all, is.na(AnyPos_T1))))	stop("unexpected NA in AnyPos_T1")
		#	check for invalid NegT and set to NA	-- we would only know that NegT is before AnyPosT and this is not helpful
		df.all	<- hivc.db.resetNegTbyAnyPosT(df.all)				
		#
		#	add first RNA Virology date
		#
		if(verbose)		cat(paste("\nadding virology data"))
		load(file.viro)
		df.cross	<- merge( subset(df.all, select=c(FASTASampleCode,Patient,AnyPos_T1,PosT,PosSeqT,NegT,NegT_Acc)), df, allow.cartesian=T, by="Patient" )
		#	checking manually AnyPos_T1>PosRNA & lRNA<3
		if(verbose)		cat(paste("\ncheck manually AnyPos_T1>PosRNA & lRNA<3 -- THIS IS ASSUMED OK including 2004G180"))
		tmp	<- subset(df.cross, AnyPos_T1>PosRNA & lRNA<3)
		tmp[,diff:=tmp[, difftime(PosRNA,AnyPos_T1, units="weeks")]]
		print( subset(tmp, diff< -10) )		
		#subset(df.cross, FASTASampleCode=="2004G180")
		#		
		#	checking manually NegT>PosRNA
		#
		if(verbose)		cat(paste("\ncheck manually NegT_Acc=='Yes' & NegT>PosRNA -- THIS IS ASSUMED OK"))
		print( subset(df.cross, NegT_Acc=="Yes" & NegT>PosRNA) )
		if(verbose)		cat(paste("\ncheck manually NegT_Acc=='No' & NegT>PosRNA -- THIS IS ASSUMED OK"))
		print( subset(df.cross, NegT_Acc=="No" & NegT>PosRNA) )
		# 	compute lRNA_T1 and lRNA_TS
		tmp			<- hivc.db.getlRNA.T1andTS(df.cross, lRNA.bTS.quantile= 0.75, lRNA.aTS.quantile= 0.25, lRNAi.min= log10(1e4), verbose=1)
		df.all		<- merge(df.all, tmp, all.x=1, by="FASTASampleCode")
		#	reset NegT by lRNA_T1
		df.all<- hivc.db.resetNegTbyPoslRNA_T1(df.all)
		# 	reset preliminary AnyPos_T1
		tmp		<- df.all[, list(AnyPos_T1=min(AnyPos_T1,PoslRNA_T1, na.rm=1)), by="FASTASampleCode"]
		if(verbose)	cat(paste("\nnumber of seq with PoslRNA_T1<AnyPos_T1, n=",length(tmp),"RESETTING"))
		set(df.all,NULL,"AnyPos_T1",tmp[,AnyPos_T1])		
		#
		#	add CD4 count data
		#
		if(verbose)		cat(paste("\nadding CD4 data"))
		load(file.immu)
		if(length(which(df[,is.na(PosCD4_T1)])))	stop("unexpected NA in PosCD4_T1")
		if(length(which(df[, is.na(PosCD4)]))) stop("unexpected NA in PosCD4")
		if(length(which(df[, is.na(CD4)]))) stop("unexpected NA in CD4")
		df.cross	<- merge( subset(df.all, select=c(FASTASampleCode,Patient,AnyPos_T1,PosT,PosSeqT,PoslRNA_T1,NegT,NegT_Acc)), subset(df, select=c(Patient,PosCD4,CD4)), allow.cartesian=T, by="Patient" )
		# delete entries where NegT>PosCD4
		df.cross	<- hivc.db.resetCD4byNegT(df.cross, with.NegT_Acc.No=1, verbose=1)
		# compute CD4_T1 and CD4_TS -- do not use on df[,CD4_T1] because some CD4 measurements might correspond to healthy patients
		tmp		<- hivc.db.getCD4.T1andTS(df.cross, CD4.HIVNeg.min= 500)
		df.all	<- merge(df.all, tmp, by="FASTASampleCode", all.x=1)
		#
		# manually checked remaining PosCD4_T1 < AnyPos_T1 -- overall plausible
		#
		#	tmp		<- subset(df.all,!is.na(PosCD4_T1) & difftime(PosCD4_T1,AnyPos_T1, units="weeks")< 0)
		#	tmp[, diff:=as.numeric(tmp[,difftime(PosCD4_T1,AnyPos_T1, units="weeks")])]
		# 	subset(tmp,diff< -10)
		#
		if(verbose)		cat(paste("\ncheck manually PosCD4_T1 < AnyPos_T1 -- THIS IS ASSUMED OK"))
		tmp	<- subset(df.all, PosCD4_T1 < AnyPos_T1, c(FASTASampleCode, Patient, AnyPos_T1, PosSeqT, PosT, PosT_Acc, PoslRNA_T1, lRNA_T1,  PosCD4_T1, CD4_T1))
		print( tmp , nrow=400)
		tmp		<- which(df.all[,!is.na(PosCD4_T1) & PosCD4_T1<AnyPos_T1])
		if(verbose)		cat(paste("\nnumber of seq with !is.na(PosCD4_T1) & PosCD4_T1<AnyPos_T1, n=",length(tmp)))
		if(verbose)		cat(paste("\nnumber of patients with !is.na(PosCD4_T1) & PosCD4_T1<AnyPos_T1, n=",length(unique(df.all[tmp,Patient]))))
		set(df.all, tmp, "AnyPos_T1", df.all[tmp,PosCD4_T1])
		#
		#	add Treatment dates
		#
		if(verbose)		cat(paste("\nadding treatment data"))
		load(file.treatment)
		df			<- subset(df, select=c(Patient, StartTime, StopTime, TrI, TrCh.failure, TrCh.adherence, TrCh.patrel, TrI.n, TrI.mo, TrI.p, AnyT_T1, AnyT_T1_Acc, HAART_T1))
		tmp			<- subset(df, select=c(Patient, TrI.n, AnyT_T1, AnyT_T1_Acc ))
		setkey(tmp, Patient)
		df.all		<- merge(df.all, unique(tmp), all.x=1, by="Patient")		
		#	compare treatment history relative to PosSeqT		
		df.cross	<- merge( subset(df.all, select=c(FASTASampleCode,Patient,PosSeqT)), df, allow.cartesian=T, by="Patient" )
		tmp			<- hivc.db.getTrIMo(df.cross)
		df.all		<- merge(df.all, tmp, all.x=1, by="FASTASampleCode")
		#		
		setkey(df.all, PosSeqT)
		setkey(df.all, Patient)
		#		
		#	manually checked remaining AnyT_T1<AnyPos_T1 -- overall plausible
		#
		if(verbose)		cat(paste("\ncheck manually AnyT_T1<AnyPos_T1 -- THIS IS ASSUMED OK"))
		tmp	<- subset(df.all, !is.na(AnyT_T1) & AnyT_T1<AnyPos_T1, c(Patient, FASTASampleCode, PosCD4_T1, CD4_T1, TrI.n, PosSeqT, AnyPos_T1,  AnyT_T1, AnyT_T1_Acc))
		tmp[, diff:=as.numeric(tmp[,difftime(AnyT_T1,AnyPos_T1, units="weeks")])]
		subset(tmp,diff< -10)
		tmp		<- which( df.all[,!is.na(AnyT_T1) & AnyT_T1<AnyPos_T1] )
		if(verbose)		cat(paste("\nnumber of seq with !is.na(AnyT_T1) & AnyT_T1<AnyPos_T1, n=",length(tmp),"SET AnyPos_T1 to lower value"))
		set(df.all, tmp, "AnyPos_T1", df.all[tmp,AnyT_T1])
		df.all[, AnyT_T1_Acc:=NULL]
		#	round numbers
		#set(df.all, NULL, "TrI.p", round(df.all[,TrI.p],d=2))
		#set(df.all, NULL, "TrI.mo", round(df.all[,TrI.mo],d=1))
		if(verbose)	cat(paste("\nsave to file",file.out))
		#setkey(df.all, FASTASampleCode)
		save(df.all,file=file.out)
		str(df.all)		
	}	
	#
	#	new diagnoses by CD4
	#
	df.newdiag		<- copy(subset(df.all, select=c(Patient,AnyPos_T1, CD4_T1)))
	setkey(df.newdiag,Patient)
	df.newdiag		<- unique(df.newdiag)
	df.newdiagCD4 	<- hivc.db.getplot.newdiagnosesbyCD4(df.newdiag, plot.file= paste(dir.name,"/derived/",file.out.name,"_NewDiagByCD4.pdf",sep=''), plot.ylab= "New diagnoses HIV-1 subtype B,\n seq available")
	#
	#	seem in care by risk group
	#
	df.living		<- copy(subset(df.all, select=c(Patient, AnyPos_T1, Trm, DateDied, DateLastContact)))
	setkey(df.living,Patient)
	df.living		<- unique(df.living)
	tmp				<- hivc.db.getplot.livingbyexposure(df.living, plot.file=paste(dir.name,"/derived/",file.out.name,"_Seen4CareByExpGroup.pdf",sep=''), plot.ylab="Seen for care with HIV-1 subtype B,\n seq available", db.endtime=2013.3, db.diff.lastcontact2died=0.5, db.diff.lastcontact2now= 2.3, verbose=1)
	
	quit("no")
	#
	#compute coverage of all data covariates
	#
	#coverage by seq
	df.covbyseq	<- hivc.db.getcoverage(df)	
	#coverage by patient
	df			<- subset(df.all, select= c(	Patient, DateBorn, Sex, CountryBorn, RegionOrigin, DateDied, Subtype, isAcute, 
												NegT, NegT_Acc, PosT, PosT_Acc, CountryInfection, Trm,  DateLastContact, RegionHospital, DateFirstEverCDCC,
												isDead, PoslRNA_T1, lRNA_T1, PosCD4_T1, CD4_T1, AnyT_T1))
	setkey(df, Patient)
	df			<- unique(df)
	df.covbypat	<- hivc.db.getcoverage(df)							
	
	quit("no")
	#
	#compute file All1stPatientCovariates
	#
	file.out		<- paste(dir.name,"derived/ATHENA_2013_03_All1stPatientCovariates.R",sep='/')
	if(resume)												#//load if there is R Master data.table
	{
		options(show.error.messages = FALSE)		
		readAttempt<-try(suppressWarnings(load(file.out)))
		if(!inherits(readAttempt, "try-error"))	cat(paste("\nresume file",file.out))
		if(!inherits(readAttempt, "try-error"))	str(df.all)
		options(show.error.messages = TRUE)		
	}
	if(!resume || inherits(readAttempt, "try-error"))		#else generate R Master data.table
	{
		#get PosSeqPROT   PosSeqRT		
		load(file.seq)
		df.seq<- df	
		df.minDateRes	<- lapply( seq_along(df.seq),function(i)
				{
					cat(paste("\nprocess", names(df.seq)[i]))
					#extract entries without missing DateRes
					df.miss<- data.table(subset(df.seq[[i]], is.na(DateRes), c(Patient,DateRes,SampleCode)), key="Patient")
					if( nrow(df.miss)!=length(unique(df.miss[,Patient])) ) stop("handling of missing DateRes not appropriate")
					if(verbose)		cat(paste("\npatients without DateRes", nrow(df.miss)))
					
					df<- subset(df.seq[[i]], !is.na(DateRes), c(Patient,DateRes,SampleCode))
					df<- data.table(df, key="Patient")								
					if(verbose)		cat(paste("\nentries without missing DateRes", nrow(df)))							
					#extract min DateRes per patient	
					df<- df[,.SD[which.min(DateRes)],by=Patient]
					if(verbose)		cat(paste("\npatients with DateRes", nrow(df)))					
					df<- rbindlist(list(df,subset( df.miss, !Patient%in%df$Patient )))
					setkey(df, Patient)
					df[, "SampleCode"]	<- factor(df[,SampleCode])
					setnames(df,"SampleCode",paste("Seq",names(df.seq)[i],sep=''))
					setnames(df,"DateRes",paste("Pos",names(df.seq)[i],sep=''))							
					df
				})
		df.minDateRes<- merge(df.minDateRes[[1]],df.minDateRes[[2]],all=1)
		if(verbose)		cat(paste("\npatients RT or PROT", nrow(df.minDateRes)))
		if(nrow(df.minDateRes)!=length(unique( df.minDateRes[,Patient] )))	stop("non-unique patients at this point")		
		#str(df.minDateRes)
	
		#add Patient data		
		load(file.patient)		
		df.all<- df[df.minDateRes]
		if(verbose)		cat(paste("\npatients in combined data table", nrow(df.all)))
		
		#add Treatment dates
		load(file.treatment)
		df		<- data.table(df, key="Patient")
		setnames(df, "T0","HAART_T1")
		df		<- subset(df, select= c(Patient, HAART_T1, StartTime))		
		if(verbose)		cat(paste("\n\nadding treatment data\nentries in regimens data table", nrow(df)))
		df		<- subset(df, !is.na(StartTime) )				
		df		<- df[, { tmp<- which.min(StartTime); list(HAART_T1=HAART_T1[tmp], AnyT_T1=StartTime[tmp])}, by=Patient]		
		if( nrow(subset(df,!is.na(HAART_T1) & HAART_T1<AnyT_T1)) )	stop("found AnyT_T1 that is older than HAART_T1")
		if(verbose)		cat(paste("\npatients with at least one non-missing treatment date", nrow(df)))
		df		<- subset(df,select=c(Patient,AnyT_T1))
		setnames(df, "AnyT_T1","AnyTherT1")		
		df.all	<- df[df.all]
		#setnames(df, "T0","PREHAART_T1")
		
		#add first RNA Virology date		
		load(file.viro)
		df		<- data.table(df, key="Patient")		
		#str(df)		
		setnames(df, "DateRNA","PosRNA")
		if(verbose)		cat(paste("\n\nadding virology\nentries in viro data table", nrow(df)))
		df		<- subset(df, !is.na(PosRNA) & Undetectable!="Yes" )
		if(verbose)		cat(paste("\nentries in viro data table without missing or non-detectable", nrow(df)))
		df		<- df[,{tmp<- which.min(PosRNA); list(PosRNA_T1= PosRNA[tmp], RNA_T1= RNA[tmp]) }, by=Patient]		
		if(verbose)		cat(paste("\npatients with at least one non-missing and detectable RNA date", nrow(df)))
		if(0)
		{
			print(range(df[,PosRNA], na.rm=1))		
		}
		df.all<- df[df.all]
		
		#add first CD4 count date		
		load(file.immu)
		df		<- data.table(df, key="Patient")
		setnames(df, "DateImm","PosCD4")
		if(verbose)		cat(paste("\nadding immunology\nentries in immu data table", nrow(df)))
		df<- subset(df, !is.na(PosCD4) )
		if(verbose)		cat(paste("\nentries in immu data table, non-missing", nrow(df)))
		df		<- df[,{tmp<- which.min(PosCD4); list(PosCD4_T1= PosCD4[tmp], CD4_T1= CD4A[tmp]) }, by=Patient]
		if(verbose)		cat(paste("\npatients with at least one non-missing CD4 date", nrow(df)))
		if(0)
		{
			print(range(df[,PosCD4_T1], na.rm=1))		
		}
		df.all<- df[df.all]
		df.all<- subset(df.all,select=c(Patient,Trm,SeqPROT,SeqRT,PosPROT,PosRT,PosT,PosT_Acc,PosCD4_T1,PosRNA_T1, NegT,NegT_Acc, AnyTherT1, isAcute,isDead,Died, RNA_T1, CD4_T1))
		
		
		#if 	PosPROT!=PosRT		consider only earliest sequence and set other to missing
		tmp<- subset(df.all, is.na(PosPROT) & is.na(PosRT))
		if(verbose)		cat(paste("\n\nkeep only earliest sequence PROT or RT if PosPROT!=PosRT\npatients with is.na(PosPROT) & is.na(PosRT):", nrow(tmp)))				
		tmp<- subset(df.all, !is.na(PosPROT) & !is.na(PosRT) & PosPROT!=PosRT)
		if(verbose)		cat(paste("\npatients with PosPROT!=PosRT & !is.na(PosPROT) & !is.na(PosRT):", nrow(tmp)))
		tmp<- subset(df.all, !is.na(PosPROT) & !is.na(PosRT) & PosPROT<PosRT, Patient)
		df.all[tmp,c("PosRT","SeqRT")]<- as.factor(NA)
		tmp<- subset(df.all, !is.na(PosPROT) & !is.na(PosRT) & PosPROT>PosRT, Patient)
		df.all[tmp,c("PosPROT","SeqPROT")]<- as.factor(NA)	
		tmp<- subset(df.all, PosPROT!=PosRT & !is.na(PosPROT) & !is.na(PosRT))
		if(verbose)		cat(paste("\npatients with PosPROT!=PosRT & !is.na(PosPROT) & !is.na(PosRT)	AFTER FIX:", nrow(tmp)))		
		#replace PosPROT, PosRT with PosSeq
		tmp					<- df.all[,PosPROT]
		tmp[ is.na(tmp) ]	<- df.all[is.na(tmp),PosRT]	
		df.all[,PosSeq:=tmp]		
		df.all<- subset(df.all, select=c(Patient,Trm,SeqPROT,SeqRT,PosSeq,PosT,PosT_Acc,PosCD4_T1,PosRNA_T1,NegT,NegT_Acc, AnyTherT1, isAcute,isDead,Died,RNA_T1, CD4_T1))
		
		
		if(verbose)		cat(paste("\nsave df.all to", file.out))				
		save(df.all,file=file.out)
		str(df.all)
	}
}
######################################################################################
project.hivc.Excel2dataframe.Regimen<- function(dir.name= DATA, verbose=1)
{
	NA.time			<- c("01/01/1911","01/11/1911","11/11/1911")	
	MAX.time		<- c("")
	TR.notyet		<- "30/03/2013"
	TR.failure 		<- c(21,31,32)
	TR.adherence	<- c(47)
	TR.patrel		<- c(23, 24, 33, 34, 36)	#either patient s decision, toxicity or interaction with other medication
	#read REGIMEN csv data file and preprocess
	file			<- paste(dir.name,"derived/ATHENA_2013_03_Regimens.csv",sep='/')
	df				<- read.csv(file, stringsAsFactors=FALSE)							
	
	date.var		<- c("T0","StartTime","StopTime")		
	for(x in date.var)
	{
		cat(paste("\nprocess Time", x))
		nok.idx			<- which( df[,x]==NA.time[1] )
		if(verbose)	cat(paste("\nentries with format ",NA.time[1],", n=", length(nok.idx), "set to NA"))
		if(length(nok.idx))	
			df[nok.idx,x]	<- NA
		nok.idx			<- which( df[,x]==NA.time[2] )
		if(verbose)	cat(paste("\nentries with format ",NA.time[2],", n=", length(nok.idx), "set to NA"))
		if(length(nok.idx))
			df[nok.idx,x]	<- NA
		nok.idx			<- which( df[,x]==NA.time[3] )
		if(verbose)	cat(paste("\nentries with format ",NA.time[3],", n=", length(nok.idx), "set to NA"))
		if(length(nok.idx))
			df[nok.idx,x]	<- NA			
		nok.idx			<- which( df[,x]==MAX.time[1] )
		if(verbose)	cat(paste("\nentries with format ",MAX.time[1],", n=", length(nok.idx), "set to max time 01/01/2999"))
		if(length(nok.idx))
			df[nok.idx,x]	<- TR.notyet
		df[,x]			<- as.Date(df[,x], format="%d/%m/%Y")	
	}
	TR.notyet				<- as.Date(TR.notyet, format="%d/%m/%Y")
	
	df						<- data.table(df, key="Patient")
	setnames(df, "T0","HAART_T1")
	set(df,NULL,"Patient",factor(df[,Patient]))
	if(verbose)	cat(paste("\nnumber of entries, n=",nrow(df)))
	#	fix entry 	M17493  1996-01-01 manually
	tmp						<- which(df[, Patient=="M17493" & StartTime=="1996-01-01"])
	if(verbose)	cat(paste("\nfix entry 	M17493  1996-01-01 manually to somewhere in 1996"))
	set(df, tmp, "StartTime", as.Date("1996-07-01"))		
	#	fix entry 	M41688  2011-01-01 manually
	tmp						<- which(df[, Patient=="M41688" & StartTime=="2011-01-01"])
	if(verbose)	cat(paste("\nfix entry 	M41688  2011-01-01 manually to 2011-11-01 (Ard)"))
	set(df, tmp, "StartTime", as.Date("2011-11-01"))	
	set(df, which(df[, Patient=="M41688"]), "HAART_T1", as.Date("2011-11-01"))		
	#	fix entry 	M41688  2011-08-23 manually
	tmp						<- which(df[, Patient=="M42092" & StartTime=="2011-08-23"])
	if(verbose)	cat(paste("\nfix entry 	M42092  2011-08-23 to 2012-08-23"))
	set(df, tmp, "StartTime", as.Date("2012-08-23"))
	set(df, tmp, "HAART_T1", as.Date("2012-08-23"))
	#	fix entry 	M41688  2010-10-23 manually
	tmp						<- which(df[, Patient=="M42186" & StartTime=="2010-10-23"])
	if(verbose)	cat(paste("\nfix entry 	M42186  2010-10-23 to 2012-10-23"))
	set(df, tmp, "StartTime", as.Date("2012-10-23"))
	set(df, tmp, "HAART_T1", as.Date("2012-10-23"))		
	#	fix entry 	M37531  2006-08-15 manually
	tmp						<- which(df[, Patient=="M37531" & StartTime=="2006-08-15"])
	if(verbose)	cat(paste("\nfix entry 	M37531  2006-08-15 to 2008-08-15"))
	set(df, tmp, "StartTime", as.Date("2008-08-15"))
	tmp						<- which(df[, Patient=="M37531"])
	set(df, tmp, "HAART_T1", as.Date("2008-08-15"))				
	#	set potentially inaccurate StartTime	
	tmp																										<- rep(0, nrow(df))
	tmp[is.na(df[,StartTime])]																				<- NA
	tmp[which( df[, !is.na(StartTime) & as.POSIXlt(StartTime)$mday==15] )]									<- 1
	tmp[which( df[, !is.na(StartTime) & as.POSIXlt(StartTime)$mon==6 & as.POSIXlt(StartTime)$mday==1 ] )]	<- 2
	if(verbose) cat(paste("\nnumber of inaccurate StartTime entries, n=",length(tmp[na.omit(tmp>0)])))
	df[,StartTime_Acc:= factor(tmp,levels=c(0,1,2),labels=c("Acc","NAccD","NAccMD"))]
	#	set potentially inaccurate StopTime
	tmp																										<- rep(0, nrow(df))
	tmp[is.na(df[,StopTime])]																				<- NA
	tmp[which( df[, !is.na(StopTime) & as.POSIXlt(StopTime)$mday==15] )]									<- 1
	tmp[which( df[, !is.na(StopTime) & as.POSIXlt(StopTime)$mon==6 & as.POSIXlt(StopTime)$mday==1 ] )]		<- 2
	if(verbose) cat(paste("\nnumber of inaccurate StopTime entries, n=",length(tmp[na.omit(tmp>0)])))
	df[,StopTime_Acc:= factor(tmp,levels=c(0,1,2),labels=c("Acc","NAccD","NAccMD"))]
	#
	#	removing Patients not on treatment
	#
	tmp						<- which(df[, is.na(StartTime) & StopTime==TR.notyet & NoDrug==0])
	if(verbose)	cat(paste("\nnumber of entries with is.na(StartTime) & StopTime==TR.notyet & NoDrug==0",length(tmp),"SETTING TO NA"))
	set(df, tmp, "StopTime", NA)
	tmp						<- which(df[, StartTime==TR.notyet & StopTime==TR.notyet])
	if(verbose)	cat(paste("\nnumber of entries with StartTime==TR.notyet & StopTime==TR.notyet",length(tmp),"SETTING TO NA"))
	set(df, tmp, "StartTime", NA)
	set(df, tmp, "StopTime", NA)
	tmp						<- which(df[, StartTime==TR.notyet])
	if(verbose)	cat(paste("\nnumber of entries with StartTime==TR.notyet & StopTime!=TR.notyet",length(tmp),"MISCLASSIFIED StartTime - setting to NA"))
	set(df, tmp, "StartTime", NA)		
	df						<- subset(df, !is.na(StopTime))
	if(verbose)	cat(paste("\nnumber of entries with !is.na(StartTime) & !is.na(StopTime), n=",nrow(df)))		
	#
	#	TR.interrupted
	#
	if(nrow(df[which(is.na(df[,NoDrug]) & StartTime!=TR.notyet),])) stop("unexpected NA in NoDrug when on treatment")
	tmp								<- rep(0, nrow(df))
	tmp[ which( df[,NoDrug==0] ) ]	<- 1
	df[, TrI:= factor(tmp,levels=c(0,1),labels=c("No","Yes"))]
	#
	#	process treatment change reasons
	#
	if(any(!is.na(df[, Reason7])))	stop("unexpected !NA after Reason7")
	df.TrCh.noreason		<- which( df[, is.na(Reason1)&is.na(Reason2)&is.na(Reason3)&is.na(Reason4)&is.na(Reason5)&is.na(Reason6)] )		
	#	TR.failure
	tmp						<-	rep(0, nrow(df))
	tmp[ df.TrCh.noreason ]	<- NA
	tmp[ which(df[, Reason1%in%TR.failure | Reason2%in%TR.failure | Reason3%in%TR.failure | Reason4%in%TR.failure | Reason5%in%TR.failure | Reason6%in%TR.failure ]) ]<- 1
	df[, TrCh.failure:= factor(tmp,levels=c(0,1),labels=c("No","Yes"))]
	#	TR.adherence
	tmp						<-	rep(0, nrow(df))
	tmp[ df.TrCh.noreason ]	<- NA
	tmp[ which(df[, Reason1%in%TR.adherence | Reason2%in%TR.adherence | Reason3%in%TR.adherence | Reason4%in%TR.adherence | Reason5%in%TR.adherence | Reason6%in%TR.adherence ]) ]<- 1
	df[, TrCh.adherence:= factor(tmp,levels=c(0,1),labels=c("No","Yes"))]
	#	TR.patient related
	tmp						<-	rep(0, nrow(df))
	tmp[ df.TrCh.noreason ]	<- NA
	tmp[ which(df[, Reason1%in%TR.patrel | Reason2%in%TR.patrel | Reason3%in%TR.patrel | Reason4%in%TR.patrel | Reason5%in%TR.patrel | Reason6%in%TR.patrel ]) ]<- 1
	df[, TrCh.patrel:= factor(tmp,levels=c(0,1),labels=c("No","Yes"))]
	#
	#	removing Patient entries before treatment started
	#
	tmp						<- which( df[,is.na(StartTime) & TrI=="Yes"]  )
	tmp						<- df[tmp,][, df.idx:=tmp] 
	tmp						<- merge(df, subset(tmp,select=c(Patient,df.idx)), by="Patient")		
	tmp						<- tmp[, 	{
				x<- matrix( as.numeric(c(StartTime, StopTime)), ncol=2 )
				x<- x[order(x[,2]),]
				list(NAStartTime.canberemoved=which(is.na(x[,1]))==1, df.idx=df.idx[1])					
			}, by=Patient]
	tmp						<- subset(tmp, NAStartTime.canberemoved)
	if(verbose)	cat(paste("\nnumber of entries with is.na(StartTime) & TrI=='Yes' before start of treatment , n=",nrow(tmp),"REMOVE"))
	set(df, tmp[,df.idx], "StopTime", NA)		
	df						<- subset(df, !is.na(StopTime))
	if(verbose)	cat(paste("\nnumber of entries with !is.na(StartTime) & !is.na(StopTime), n=",nrow(df)))
	#
	#	fix inconsistent timings
	#
	tmp		<- which(df[,StartTime>StopTime])		
	tmp		<- cbind( tmp, sapply(tmp,function(x)		which(df[, Patient==df[x,Patient] & StopTime==df[x,StartTime]])	) )		#second col contains StopTime
	for(i in seq_len(nrow(tmp)))
	{
		if(verbose)	cat(paste("\nprocess",df[tmp[i,1],Patient]))
		if(df[tmp[i,1],StartTime_Acc=="Acc"])			#set StartTime_Acc=="NAccD"
		{
			set(df,tmp[i,1],"StartTime_Acc","NAccD")			
			set(df,tmp[i,2],"StopTime_Acc","NAccD")			
		}
		if(df[tmp[i,1],StartTime_Acc=="NAccD"])			#reset StopTime		StartTime_Acc=="NAccD"
		{
			z		<- as.POSIXlt(df[tmp[i,2],StopTime])
			z$mday	<- 1
			z		<- as.Date(z)
			z		<- max( df[tmp[i,2],StartTime],z )
			set(df,tmp[i,1],"StartTime",z)
			set(df,tmp[i,2],"StopTime",z)
		}
		if(df[tmp[i,1],StartTime]>df[tmp[i,1],StopTime])
		{
			set(df,tmp[i,1],"StartTime_Acc","NAccMD")			
			set(df,tmp[i,2],"StopTime_Acc","NAccMD")
		}
		if(df[tmp[i,1],StartTime_Acc=="NAccMD"])		#set to midpoint
		{
			z<- df[tmp[i,2],StartTime] + difftime(df[tmp[i,1],StopTime],df[tmp[i,2],StartTime],units="days") / 2
			set(df,tmp[i,1],"StartTime",z)
			set(df,tmp[i,2],"StopTime",z)
		}		
	}
	#
	#	get 1st sort by Patient, 2nd sort by StopTime -> treatment history by Patient is now in order
	#
	setkey(df, StopTime)
	setkey(df, Patient)	
	#subset( df,StartTime_Acc=="NAccMD" & StopTime_Acc!="NAccMD")
	#subset( df, Patiet=="M10027")		
	#
	#	simple statistics of Patient history: number of treatment interruptions, total length of interruption in months and proportion of treatment interruption in treatment history 		
	#
	tmp		<- df[, 	{ 
				x		<- data.table( StartTime, StopTime, TrI, StartTime_Acc )
				x		<- subset(x, !is.na(StartTime))		#discard first time period if unknown StartTime
				z		<- subset(x, TrI=="Yes")
				list( 	TrI.n		= nrow(z), 
						TrI.mo		= as.numeric(sum(z[,difftime(StopTime,StartTime,units="days")/30])), 
						TrI.p		= as.numeric(sum(z[,difftime(StopTime,StartTime,units="days")/30])) / as.numeric(sum(x[,difftime(StopTime,StartTime,units="days")/30])), 
						HAART_T1	= HAART_T1[1], 
						AnyT_T1		= subset(x, TrI=="No")[,StartTime][1], 
						AnyT_T1_Acc	= subset(x, TrI=="No")[,StartTime_Acc][1]		)																					
			}, by=Patient]
	if( nrow(subset(tmp,!is.na(HAART_T1) & HAART_T1<AnyT_T1)) )	stop("found AnyT_T1 that is older than HAART_T1")		
	df		<- merge(subset(tmp,select=c(Patient,AnyT_T1,AnyT_T1_Acc,TrI.n, TrI.mo, TrI.p)), df, all.y=1,by="Patient")
	#
	#	reset AnyT_T1 conservatively
	#
	nacc				<- which(df[,AnyT_T1_Acc=="NAccD"])		
	tmp					<- as.POSIXlt(df[nacc,AnyT_T1] )
	tmp$mday			<- 30
	set(df, nacc, "AnyT_T1", as.Date(tmp))
	nacc				<- which(df[,AnyT_T1_Acc=="NAccMD"])
	tmp					<- as.POSIXlt(df[nacc,AnyT_T1] )
	tmp$mday			<- 31
	tmp$mon				<- 11
	set(df, nacc, "AnyT_T1", as.Date(tmp))		
	#
	file		<- paste(dir.name,"derived/ATHENA_2013_03_Regimens.R",sep='/')
	if(verbose) cat(paste("\nsave to", file))
	save(df, file=file)
}
######################################################################################
project.hivc.Excel2dataframe.CD4<- function(dir.name= DATA, verbose=1)
{
	NA.time			<- c("","01/01/1911","11/11/1911","24/06/1923")		
	verbose			<- 1
	#read CD4 csv data file and preprocess
	file			<- paste(dir.name,"derived/ATHENA_2013_03_Immu.csv",sep='/')
	df				<- read.csv(file, stringsAsFactors=FALSE)											
	date.var		<- c("DateImm")		
	for(x in date.var)
	{
		cat(paste("\nprocess Time", x))
		nok.idx			<- which( df[,x]==NA.time[1] )
		if(verbose)	cat(paste("\nentries with format ",NA.time[1],", n=", length(nok.idx)))
		#if(verbose && length(nok.idx))	cat(paste("\nentries with format ",NA.time[1],", Patient", paste(df[nok.idx,"Patient"],collapse=', ')))			
		if(length(nok.idx))	
			df[nok.idx,x]	<- NA
		nok.idx			<- which( df[,x]==NA.time[2] )
		if(verbose)	cat(paste("\nentries with format ",NA.time[2],", n=", length(nok.idx)))
		if(length(nok.idx))
			df[nok.idx,x]	<- NA
		nok.idx			<- which( df[,x]==NA.time[3] )
		if(verbose)	cat(paste("\nentries with format ",NA.time[3],", n=", length(nok.idx)))
		if(length(nok.idx))
			df[nok.idx,x]	<- NA
		nok.idx			<- which( df[,x]==NA.time[4] )
		if(verbose)	cat(paste("\nentries with format ",NA.time[4],", n=", length(nok.idx)))
		if(length(nok.idx))
			df[nok.idx,x]	<- NA
		df[,x]			<- as.Date(df[,x], format="%d/%m/%Y")	
	}		
	df		<- data.table(df, key="Patient")
	setnames(df, "DateImm","PosCD4")
	set(df, NULL, "Patient", factor(df[,Patient]))
	if(verbose) cat(paste("\nnumber of entries read, n=",nrow(df)))
	tmp		<- which(df[,is.na(CD4A)])
	if(verbose) cat(paste("\nnumber of entries with is.na(CD4A), n=",length(tmp),"SETTING PosCD4 to NA"))
	set(df, tmp, "PosCD4", NA)
	df		<- subset(df,!is.na(PosCD4))
	if(verbose) cat(paste("\nnumber of entries with !is.na(PosCD4), n=",nrow(df)))
	#
	#	data corrections from Ard
	#
	tmp		<- which(df[, Patient=="M11392" & PosCD4=="2000-10-30" & CD4A==3001])
	if(verbose) cat(paste("\nnumber of entries with corrected CD4 units Patient=='M11392' & CD4A>3001, n=",length(tmp),"double entry, remove"))
	set(df, tmp, "CD4A", NA)	
	tmp		<- which(df[, Patient=="M26334" & PosCD4=="2008-09-11" & CD4A==2230])
	if(verbose) cat(paste("\nnumber of entries with corrected CD4 units Patient==M26334 & PosCD4==2008-09-11 & CD4A==2230, n=",length(tmp),"double entry, remove"))
	set(df, tmp, "CD4A", NA)		
	tmp		<- which(df[, 	Patient=='M26293' & PosCD4=='2003-09-12' & CD4A==800])
	if(verbose) cat(paste("\nnumber of entries with corrected PosCD4  Patient=='M26293' & PosCD4=='2003-09-12' & CD4A==800, n=",length(tmp),"set to 16/9/2003"))
	set(df, tmp, "PosCD4", as.Date('2003-09-16'))
	tmp		<- which(df[, 	Patient=='M11233' & PosCD4=='2001-01-10' & CD4A==24])
	if(verbose) cat(paste("\nnumber of entries with corrected PosCD4  Patient=='M11233' & PosCD4=='2001-01-10' & CD4A==24, n=",length(tmp),"set to 25/4/1996"))
	set(df, tmp, "PosCD4", as.Date('1996-04-25'))
	tmp		<- which(df[, 	Patient=='M13124' & PosCD4=='1998-08-07'])
	if(verbose) cat(paste("\nnumber of entries with corrected PosCD4  Patient=='M13124' & PosCD4=='1998-08-07', n=",length(tmp),"set to 7/9/1998"))
	set(df, tmp, "PosCD4", as.Date('1998-09-07'))	
	tmp		<- which(df[, 	Patient=='M13124' & PosCD4=='2001-09-14'])
	if(verbose) cat(paste("\nnumber of entries with corrected PosCD4  Patient=='M13124' & PosCD4=='2001-09-14', n=",length(tmp),"set to 14/9/2000"))
	set(df, tmp, "PosCD4", as.Date('2000-09-14'))
	tmp		<- which(df[, 	Patient=="M10544" & PosCD4=="2003-02-17" & CD4A==850	|
							Patient=="M11099" & PosCD4=="1997-12-30" & CD4A==1240	|
							Patient=="M11133" & PosCD4=="2003-06-16" & CD4A==170	|
							Patient=="M11137" & PosCD4=="2003-06-25" & CD4A==460	|
							Patient=="M11167" & PosCD4=="2006-09-04" & CD4A==400	|	
							Patient=="M11351" & PosCD4=="1996-10-08" & CD4A==150	|
							Patient=="M11713" & PosCD4=="1996-07-03" & CD4A==0.37	|
							Patient=="M12577" & PosCD4=="2000-09-25" & CD4A==210	|
							Patient=="M12884" & PosCD4=="1997-11-04" & CD4A==350	|
							Patient=="M13051" & PosCD4=="1998-06-08" & CD4A==460	|
							Patient=="M13124" & PosCD4=="2001-10-09" & CD4A==1.17	|
							Patient=="M13124" & PosCD4=="2003-02-12" & CD4A==0.74	|
							Patient=="M13124" & PosCD4=="2003-03-21" & CD4A==0.59	|
							Patient=="M13124" & PosCD4=="2003-06-17" & CD4A==0.61	|
							Patient=="M13126" & PosCD4=="2001-01-08" & CD4A==0.5	|
							Patient=="M13126" & PosCD4=="2001-03-05" & CD4A==0.43	|
							Patient=="M13126" & PosCD4=="2003-01-17" & CD4A==0.48	|
							Patient=="M13126" & PosCD4=="2003-04-28" & CD4A==0.48	|
							Patient=="M13126" & PosCD4=="2003-07-24" & CD4A==0.46	|	#no CD4 anymore at 24/7/2003, delete both
							Patient=="M13126" & PosCD4=="2003-07-24" & CD4A==460	|
							Patient=="M13298" & PosCD4=="1997-12-11" & CD4A==760	|
							Patient=="M14834" & PosCD4=="1998-12-10" & CD4A==990	|
							Patient=="M14945" & PosCD4=="2000-09-20" & CD4A==1620	|
							Patient=="M14986" & PosCD4=="2001-03-29" & CD4A==640	|
							Patient=="M14995" & PosCD4=="1999-01-12" & CD4A==370	|
							Patient=="M15071" & PosCD4=="1999-10-12" & CD4A==100	|
							Patient=="M15234" & PosCD4=="1998-09-02" & CD4A==25		|
							Patient=="M15234" & PosCD4=="1998-11-04" & CD4A==32		|
							Patient=="M15234" & PosCD4=="1998-12-16" & CD4A==35		|
							Patient=="M16018" & PosCD4=="1998-09-09" & CD4A==1400	|
							Patient=="M16570" & PosCD4=="2003-06-17" & CD4A==280	|
							Patient=="M16622" & PosCD4=="2000-09-27" & CD4A==100	|
							Patient=="M17154" & PosCD4=="2011-01-05" & CD4A==495	|
							Patient=="M17819" & PosCD4=="2000-03-15" & CD4A==1010	|
							Patient=="M17951" & PosCD4=="1999-07-22" & CD4A==530	|
							Patient=="M18712" & PosCD4=="2000-07-31" & CD4A==600	|
							Patient=="M25530" & PosCD4=="2002-04-09" & CD4A==59		|
							Patient=="M25530" & PosCD4=="2002-04-19" & CD4A==66		|	
							Patient=="M28189" & PosCD4=="2007-09-04" & CD4A==31		|
							Patient=="M30605" & PosCD4=="2011-08-22" & CD4A==83		|
							Patient=="M31573" & PosCD4=="2011-03-24" & CD4A==0		|	#import from hospital db, remove unlikely entry
							Patient=="M33353" & PosCD4=="2006-03-22" & CD4A==101	|
							Patient=="M33924" & PosCD4=="2007-11-01" & CD4A==0		|
							Patient=="M37294" & PosCD4=="2011-07-18" & CD4A==820	|
							Patient=="M39055" & PosCD4=="2012-02-06" & CD4A==6850							
							])
	if(verbose) cat(paste("\nnumber of entries with incorrect CD4  should be 45, n=",length(tmp),"set to NA"))	
	set(df, tmp, "CD4A", NA)	
	#	adjust likely wrong units > 20000
	tmp		<- which(df[, CD4A>20000])
	if(verbose) cat(paste("\nnumber of entries with likely wrong CD4 units > 20000, n=",length(tmp),"DIVIDE BY 1e3"))
	if(verbose) print(df[tmp,])
	set(df, tmp, "CD4A", round(df[tmp, CD4A]*1e-3,d=0))
	#	adjust likely wrong units > 10000
	if(verbose) cat(paste("\npatient M11368"))
	tmp		<- which(df[, CD4A>10000 & Patient=="M11368"])
	set(df, tmp, "CD4A", round(df[tmp, CD4A]*1e-3,d=0))
	if(verbose) cat(paste("\npatient M32958"))
	tmp		<- which(df[, Patient=="M32958" & PosCD4<="2010-05-03"])
	set(df, tmp, "CD4A", round(df[tmp, CD4A]*1e-1,d=0))
	if(verbose) cat(paste("\npatient M20633"))
	tmp		<- which(df[, Patient=="M20633" & CD4A>5000])
	set(df, tmp, "CD4A", round(df[tmp, CD4A]*1e-3,d=0))
	tmp		<- which(df[, CD4A>3000])
	if(verbose) cat(paste("\nnumber of entries with likely wrong CD4 units > 3000, n=",length(tmp),"DIVIDE BY 10"))
	if(verbose) print(df[tmp,])
	set(df, tmp, "CD4A", round(df[tmp, CD4A]*1e-1,d=0))
	#
	#	check above 1700 manually
	#
	tmp<- merge(df, subset(df, CD4A>1700, Patient), by="Patient")
	tmp<- tmp[,	list(CD4.med= median(CD4A), CD4.max=max(CD4A)),by="Patient"]
	print( subset(tmp, CD4.med*1.5<CD4.max) )
	#	divide by 10
	tmp		<- which(df[, Patient%in%c("M10212","M14927","M15431","M15519","M20720","M26334","M27643","M27571") & CD4A>1700])
	if(verbose) cat(paste("\nnumber of entries with likely wrong CD4 units > 1700  M10212 M27571 M14927  M15431  M15519  M20720  M26334  M27643, n=",length(tmp),"DIVIDE BY 10"))
	set(df, tmp, "CD4A", round(df[tmp, CD4A]*1e-1,d=0))
	if(verbose) cat(paste("\npatient M17554"))
	tmp		<- which(df[, Patient=="M17554" & CD4A%in%c(2760,2170)])
	set(df, tmp, "CD4A", round(df[tmp, CD4A]*1e-1,d=0))
	if(verbose) cat(paste("\npatient M17554"))
	tmp		<- which(df[, Patient=="1500" & CD4A%in%c(1500)])
	set(df, tmp, "CD4A", round(df[tmp, CD4A]*1e-1,d=0))
	#	set to NA
	tmp		<- which(df[, Patient%in%c("M12953","M13340","M26537","M26669","M35668") & CD4A>1900])
	if(verbose) cat(paste("\nnumber of entries with likely wrong CD4 units > 1900  M12953   M13340  M26537  M26669  M35668, n=",length(tmp),"SET to NA"))
	set(df, tmp, "CD4A", NA)
	tmp		<- which(df[, Patient%in%c("M15743") & CD4A>2500])
	if(verbose) cat(paste("\nnumber of entries with likely wrong CD4 units > 2500  M15743, n=",length(tmp),"SET to NA"))
	set(df, tmp, "CD4A", NA)
	tmp		<- which(df[, Patient%in%c("M30155") & CD4A>1000])
	if(verbose) cat(paste("\nnumber of entries with likely wrong CD4 units > 1000  M30155, n=",length(tmp),"SET to NA"))
	set(df, tmp, "CD4A", NA)
	#
	df		<- subset(df, !is.na(CD4A))
	if(verbose) cat(paste("\nnumber of entries with !is.na(PosCD4), n=",nrow(df)))
	#
	#	check for too small entries
	#
	tmp		<- which(df[,CD4A<1 & CD4A>0])
	if(verbose) cat(paste("\nnumber of entries with too small CD4 CHECK MANUALLY, n=",length(tmp),"SET TO *1e3"))
	if(verbose)	print(df[tmp,])
	set(df, tmp, "CD4A", df[tmp,CD4A]*1e3)
	#
	#	check for double entries
	#
	tmp		<- df[,	{
						z	<- which(as.numeric(difftime(PosCD4[-1],PosCD4[-length(PosCD4)],units="days"))==0)
						z	<- rbind(z,z+1)
						z2	<- rep( apply( z, 2, function(z2)	abs(CD4A[z2[1]]-CD4A[z2[2]])	), each=2 )
						z3	<- rep( apply( z, 2, function(z2)	mean(c(CD4A[z2[1]],CD4A[z2[2]]))	), each=2 )
						z	<- as.numeric(z)
						list( PosCD4=PosCD4[z], CD4A=CD4A[z], CD4d=z2, CD4mean=z3 ) 					
					}, by="Patient"]
	if(verbose) cat(paste("\nfound double entries, n=",nrow(tmp),"SET TO MEAN -- NOT ALWAYS OK"))
	if(verbose)	print(tmp)	
	for( i in seq.int(1,nrow(tmp),2))
	{
		z<- which(df[, Patient==tmp[i,Patient] & PosCD4==tmp[i,PosCD4]])
		set(df,z[1],"CD4A", mean(df[z,CD4A]))
		set(df,z[-1],"CD4A", NA)		
	}
	df		<- subset(df, !is.na(CD4A))
	if(verbose) cat(paste("\nnumber of entries with !is.na(PosCD4), n=",nrow(df)))
	#	
	setkey(df, PosCD4)
	setkey(df, Patient)
	#
	df		<- df[, 	{
				z<- which.min(PosCD4)
				list(PosCD4=PosCD4, CD4=CD4A, PosCD4_T1=PosCD4[z], CD4_T1=CD4A[z] ) 	
			},by=Patient]
	
	
	
	file		<- paste(dir.name,"derived/ATHENA_2013_03_Immu.R",sep='/')
	if(verbose) cat(paste("\nsave to", file))
	save(df, file=file)		
}
######################################################################################
project.hivc.Excel2dataframe.Viro<- function()		
{
	verbose				<- 1
	dir.name			<- DATA
	DB.locktime			<- HIVC.db.locktime
	
	#need for checking of VL data
	file.treatment		<- paste(dir.name,"derived/ATHENA_2013_03_Regimens.R",sep='/')
	load(file.treatment)
	df.treat			<- subset(df, select=c(Patient, StartTime, StopTime, AnyT_T1, TrI, TrCh.failure, TrCh.adherence, TrCh.patrel))
	
	NA.time				<- c("","01/01/1911","11/11/1911","24/06/1923")		
	RNA.min				<- 400	#seems to be standard value
	RNA.max				<- 5e6
	lRNA.min.infectious	<- log10(1e4)
	lRNA.min.early		<- log10(1e5)
	lRNA.max.b4early	<- log10(2e4)
	
	#read VIROLOGY csv data file and preprocess
	file			<- paste(dir.name,"derived/ATHENA_2013_03_Viro.csv",sep='/')
	df				<- read.csv(file, stringsAsFactors=FALSE)									
	df$Undetectable	<- factor(df$Undetectable, levels=c(0,1,2),labels=c("No","Yes","LargerThan"))
	date.var		<- c("DateRNA")		
	for(x in date.var)
	{
		cat(paste("\nprocess Time", x))
		nok.idx			<- which( df[,x]==NA.time[1] )
		if(verbose)	cat(paste("\nentries with format ",NA.time[1],", n=", length(nok.idx)))
		#if(verbose && length(nok.idx))	cat(paste("\nentries with format ",NA.time[1],", Patient", paste(df[nok.idx,"Patient"],collapse=', ')))			
		if(length(nok.idx))	
			df[nok.idx,x]	<- NA
		nok.idx			<- which( df[,x]==NA.time[2] )
		if(verbose)	cat(paste("\nentries with format ",NA.time[2],", n=", length(nok.idx)))
		if(length(nok.idx))
			df[nok.idx,x]	<- NA
		nok.idx			<- which( df[,x]==NA.time[3] )
		if(verbose)	cat(paste("\nentries with format ",NA.time[3],", n=", length(nok.idx)))
		if(length(nok.idx))
			df[nok.idx,x]	<- NA
		nok.idx			<- which( df[,x]==NA.time[4] )
		if(verbose)	cat(paste("\nentries with format ",NA.time[4],", n=", length(nok.idx)))
		if(length(nok.idx))
			df[nok.idx,x]	<- NA
		df[,x]			<- as.Date(df[,x], format="%d/%m/%Y")	
	}		
	
	df		<- data.table(df, key="Patient")
	if(verbose) cat(paste("\nnumber of entries read, n=",nrow(df)))
	setnames(df, "DateRNA","PosRNA")
	set(df, NULL, "Patient", factor(df[,Patient]))
	#
	#	checking manually NegT>PosRNA
	#
	if(verbose)		cat(paste("\nset entry Patient=='M38400' & as.character(PosRNA)=='2005-11-24' to NA -- seems unlikely"))
	if(verbose)		print( subset(df, Patient=="M38400") )
	tmp		<- which( df[, Patient=="M38400" & as.character(PosRNA)=="2005-11-24"] )
	if(verbose)		cat(paste("\nsetting number of entries to NA, n=",length(tmp)))
	set(df,tmp,"PosRNA",NA)		
	if(verbose)		cat(paste("\nset entry Patient=='M36146' & as.character(PosRNA)=='2005-11-19' to NA -- seems unlikely"))
	if(verbose)		print( subset(df, Patient=="M36146") )
	tmp		<- which( df[, Patient=="M36146" & as.character(PosRNA)=="2005-11-19"] )
	if(verbose)		cat(paste("\nsetting number of entries to NA, n=",length(tmp)))
	set(df,tmp,"PosRNA",NA)
	# remove is.na(PosRNA) and !is.na(RNA)
	df		<- subset(df, !is.na(PosRNA) & !is.na(RNA))
	if(verbose)		cat(paste("\nnumber of entries with !is.na(PosRNA) & !is.na(RNA), n=",nrow(df)))		
	#
	#	#checking for duplicates -- do not matter that much - leave if any
	#
	#subset(df[,	{
	#				z<- diff(RNA); tmp<- (RNA>1e4 & !is.na(Undetectable))[-1]; length(which(z[tmp]==0))
	#			}				
	#			,by="Patient"], V1>0 )
	#
	#	combine Undetectable=="LargerThan" with Undetectable=="No"
	#
	tmp		<- which(df[,Undetectable=="LargerThan"])
	if(verbose)		cat(paste("\nsetting Undetectable=='LargerThan' to Undetectable=='No', n=",length(tmp)))
	set(df,tmp,"Undetectable","No")
	if(any(df[,is.na(Undetectable)]))	stop("unexpected is.na(Undetectable)")		
	set(df, NULL, "Undetectable",factor(as.numeric(df[,Undetectable]), levels=c(1,2),labels=c("No","Yes")))		
	#
	#	set Undetectable=="Yes" and RNA<RNA.min to Undetectable=="No" and RNA.min
	#
	tmp<- which( df[, Undetectable=="Yes" & RNA<1e4] )
	if(verbose)		cat(paste("\nsetting Undetectable=='Yes' and RNA<RNA.min to Undetectable=='No' and RNA.min, n=",length(tmp)))
	set(df,tmp,"Undetectable","No")
	set(df,tmp,"RNA",RNA.min)
	#
	#	set RNA<RNA.min to RNA.min
	#		
	tmp<- which( df[, RNA<RNA.min] )
	if(verbose)		cat(paste("\nsetting RNA<RNA.min to RNA.min, n=",length(tmp)))		
	set(df,tmp,"RNA",RNA.min)
	#
	#	wrong units ? adjust manually
	#				
	tmp		<- which(df[, Patient%in%c("M11314","M11331","M40782","M14759") & RNA>5e6])
	if(verbose) cat(paste("\nnumber of entries with likely wrong RNA units > 5e6  -- M11314  M11331  M40782  M14759, n=",length(tmp),"SET to 5e6"))
	set(df, tmp, "RNA", 5e6)
	#
	tmp		<- which(df[, Patient=="M27377" & PosRNA>"2007-08-11"])
	if(verbose) cat(paste("\nnumber of entries with likely wrong RNA units -- patient M27377 after 2007-08-11, n=",length(tmp),"DIV by 1e3"))
	set(df, tmp, "RNA", df[tmp,RNA]*1e-3)
	#
	tmp		<- which(df[, Patient%in%c("M13134","M18385","M18767","M20308","M35814","M35852","M36515","M41877") & RNA>1e6])
	if(verbose) cat(paste("\nnumber of entries with likely wrong RNA units > 1e6  -- M13134  M18385  M18767  M20308  M35814  M35852  M36515  M41877, n=",length(tmp),"DIV by 10"))
	set(df, tmp, "RNA", df[tmp,RNA]*1e-1)
	#
	tmp		<- which(df[, Patient%in%c("M38031") & RNA>1e5])
	if(verbose) cat(paste("\nnumber of entries with likely wrong RNA units > 1e5  -- M38031, n=",length(tmp),"DIV by 10"))
	set(df, tmp, "RNA", df[tmp,RNA]*1e-1)
	#
	tmp		<- which(df[, Patient%in%c("M38036","M33131","M33668","M34200","M34302","M20350") & RNA>1e6])
	if(verbose) cat(paste("\nnumber of entries with likely wrong RNA units > 1e6  -- M38036  M33131  M33668  M34200  M34302 M20350, n=",length(tmp),"DIV by 100"))
	set(df, tmp, "RNA", df[tmp,RNA]*1e-2)
	#
	tmp		<- which(df[, Patient%in%c("M11995","M13266","M14486","M15621","M16588") & RNA>5e6])
	if(verbose) cat(paste("\nnumber of entries with likely wrong RNA units > 5e6  -- M11995 M13266 M14486 M15621 M16588, n=",length(tmp),"DIV by 10"))
	set(df, tmp, "RNA", df[tmp,RNA]*1e-1)		
	#
	tmp		<- which(df[, Patient%in%c("M16570") & RNA>2e6])
	if(verbose) cat(paste("\nnumber of entries with likely wrong RNA units > 2e6  -- M16570, n=",length(tmp),"DIV by 10"))
	set(df, tmp, "RNA", df[tmp,RNA]*1e-1)
	#
	tmp		<- which(df[, Patient%in%c("M17655","M37746","M30733","M27377") & RNA>2e6])
	if(verbose) cat(paste("\nnumber of entries with likely wrong RNA units > 2e6  -- M17655 M37746 M30733 M27377, n=",length(tmp),"NA"))
	set(df, tmp, "RNA", NA)
	#
	tmp		<- which(df[, Patient%in%c("M17554") & RNA>1e6])
	if(verbose) cat(paste("\nnumber of entries with likely wrong RNA units > 1e6  -- M17554, n=",length(tmp),"DIV by 100"))
	set(df, tmp, "RNA", df[tmp,RNA]*1e-2)
	#
	tmp		<- which(df[, Patient%in%c("M10607","M10969","M11428","M28707","M31388","M31455","M32401","M32877","M33406","M33839","M33918","M34062","M35280","M30788") & RNA>=1e6])
	if(verbose) cat(paste("\nnumber of entries with likely wrong RNA units > 1e6  -- M10607 M10969 M11428 M28707 M31388 M31455 M32401 M32877 M33406 M33839 M33918 M34062 M35280  M30788, n=",length(tmp),"DIV by 10"))
	set(df, tmp, "RNA", df[tmp,RNA]*1e-1)
	#
	# double check entries manually in range >5e6
	#		
	tmp		<- merge(subset(df, RNA>5e6, c(Patient,PosRNA)), df.treat, by="Patient")
	tmp		<- tmp[, 	{
				z<- which( difftime(StopTime,PosRNA,units="days")>0 )
				z<- z[1]
				list(PosRNA=PosRNA[z], StartTime=StartTime[z], StopTime=StopTime[z], TrI=TrI[z], TrCh.failure=TrCh.failure[z], TrCh.adherence=TrCh.adherence[z], TrCh.patrel=TrCh.patrel[z])				
			}, by="Patient"]
	setnames(tmp, "PosRNA","PosRNAh")
	tmp		<- merge(tmp, df,by="Patient")		
	tmp[, print(data.table(Patient,RNA,PosRNA,PosRNAh,StartTime,StopTime,TrI,TrCh.failure,TrCh.adherence,TrCh.patrel)), by="Patient"]
	#
	# double check entries manually in range >RNA<5e6 & RNA>2e6
	#		
	tmp		<- merge(subset(df, RNA<5e6 & RNA>2e6, c(Patient,PosRNA)), df.treat, by="Patient")
	tmp		<- tmp[, 	{
				z<- which( difftime(StopTime,PosRNA,units="days")>0 )
				z<- z[1]
				list(PosRNA=PosRNA[z], StartTime=StartTime[z], StopTime=StopTime[z], TrI=TrI[z], TrCh.failure=TrCh.failure[z], TrCh.adherence=TrCh.adherence[z], TrCh.patrel=TrCh.patrel[z])				
			}, by="Patient"]
	setnames(tmp, "PosRNA","PosRNAh")
	tmp		<- merge(tmp, df,by="Patient")		
	tmp[, print(data.table(Patient,RNA,PosRNA,PosRNAh,StartTime,StopTime,TrI,TrCh.failure,TrCh.adherence,TrCh.patrel)), by="Patient"]
	#
	tmp		<- which(df[, Patient=="M12612" & PosRNA=="1996-06-13"])
	if(verbose) cat(paste("\nset  M12612 1996-06-13 to NA, n=",length(tmp)))
	set(df, tmp, "RNA", NA)
	tmp		<- which(df[, Patient=="M17044" & PosRNA=="2004-07-19"])
	if(verbose) cat(paste("\nset  M17044 2004-07-19 to NA, n=",length(tmp)))
	set(df, tmp, "RNA", NA)
	#
	# check entries manually with Undetectable=="Yes" & RNA>1e4
	#		
	tmp		<- merge(subset(df, Undetectable=="Yes" & RNA>1e4, c(Patient,PosRNA)), df.treat, by="Patient")
	tmp		<- tmp[, 	{
				z<- which( difftime(StopTime,PosRNA,units="days")>0 )
				z<- z[1]
				list(PosRNA=PosRNA[z], StartTime=StartTime[z], StopTime=StopTime[z], TrI=TrI[z], TrCh.failure=TrCh.failure[z], TrCh.adherence=TrCh.adherence[z], TrCh.patrel=TrCh.patrel[z])				
			}, by="Patient"]
	setnames(tmp, "PosRNA","PosRNAh")
	tmp		<- merge(tmp, df,by="Patient")		
	tmp[, print(data.table(Patient,RNA,PosRNA,Undetectable,StartTime,StopTime,TrI,TrCh.failure,TrCh.adherence,TrCh.patrel)), by="Patient"]
	#
	tmp		<- which(df[, Patient=="M30584" & PosRNA=="2004-10-14"])
	if(verbose) cat(paste("\nset  M30584 2004-10-14 to NA, n=",length(tmp)))
	set(df, tmp, "RNA", NA)
	tmp		<- which(df[, Patient=="M26369" & PosRNA=="2003-10-13"])[1]
	if(verbose) cat(paste("\nset  M26369 2003-10-13 to NA, n=",length(tmp)))
	set(df, tmp, "RNA", NA)
	tmp		<- which(df[, Patient=="M16566" & PosRNA=="2000-11-16"])
	if(verbose) cat(paste("\nset  M16566 2000-11-16 to NA, n=",length(tmp)))
	set(df, tmp, "RNA", NA)
	tmp		<- which(df[, Patient=="M12818" & PosRNA=="1996-09-18" & RNA>=1e4])
	if(verbose) cat(paste("\nset  M12818 1996-09-18 to NA, n=",length(tmp)))
	set(df, tmp, "RNA", NA)		
	#
	#	set Undetectable=='Yes' & RNA>1e4
	#
	tmp		<- which(df[, Undetectable=="Yes" & RNA>=1e4])
	if(verbose) cat(paste("\nsetting Undetectable=='Yes' & RNA>1e4 to Undetectable=='No', n=",length(tmp)))
	set(df, tmp, "Undetectable", "No")		
	#
	#	set RNA>RNA.max to RNA.max
	#		
	tmp		<- which(df[, RNA>RNA.max])
	if(verbose) cat(paste("\nsetting RNA>RNA.max to RNA.max, n=",length(tmp)))
	set(df, tmp, "RNA", RNA.max)		
	#
	df		<- subset(df,!is.na(RNA), select=c(Patient, PosRNA, RNA))
	if(verbose) cat(paste("\nentries with !is.na(RNA), n=",nrow(df)))
	#
	#	average duplicate entries out
	#						
	if(verbose) cat(paste("\nremoving duplicate entries"))
	df		<- df[,	{
				z		<- as.numeric(difftime(PosRNA[-1], PosRNA[-length(PosRNA)],units="days"))					
				if(length(which(z==0)))
				{
					RNA.d				<- sapply(which(z==0), function(i)  mean( RNA[seq.int(i,i+1)] ) )					
					PosRNA.d			<- PosRNA[z!=0]		#keep the last of the duplicates
					z2					<- RNA				
					z2[which(z==0)+1]	<- RNA.d			#overwrite the last of the duplicates
					RNA.d				<- z2[z!=0]			#keep the last of the duplicates
					#print(z); print(which(z==0)); print(RNA.d); print(list( PosRNA= PosRNA.d, RNA= RNA.d))
				}
				else
					ans<- list( PosRNA= PosRNA, RNA= RNA)
				ans					
			},by="Patient"]
	if(verbose) cat(paste("\nnumber of entries, n=",nrow(df)))
	setkey(df, PosRNA)
	setkey(df, Patient)
	#
	df[,"lRNA":=round(log10( df[,RNA] ), d=3)]		
	#	
	#	compute several statistics on lRNA life history
	#	PoslRNA_T1		time of first lRNA
	#	lRNA_T1			first lRNA 
	#	lRNA.i 			proportion of time spent above 'lRNA.min.infectious'
	#	lRNA.hb4tr_LT 	last time of lRNA above 'lRNA.min.early'
	#	lRNA.early		is there increasing lRNA before treatment 
	#	
	df					<- df[, {		
									z<- which.min(PosRNA)
									list(PosRNA=PosRNA, RNA=RNA, lRNA=lRNA, PoslRNA_T1=PosRNA[z], lRNA_T1=lRNA[z])
								},by=Patient]
	tmp					<- subset(df.treat,select=c(Patient,AnyT_T1))
	setkey(tmp,Patient)		 
	df					<- merge(df, unique(tmp), all.x=1, by="Patient")
	set(df,which(df[,is.na(AnyT_T1)]),"AnyT_T1",DB.locktime)
	df[, lRNA.infectious:=lRNA>=lRNA.min.infectious]
	df[, lRNA.high.b4tr	:=lRNA>=lRNA.min.early & PosRNA<AnyT_T1]
	
	tmp	<- df[,		{
						z				<- data.table(PosRNA, StopRNA=c(PosRNA[-1],DB.locktime), lRNA, AnyT_T1, lRNA.infectious, lRNA.high.b4tr)
						p.infectious	<- z[,as.numeric(difftime(StopRNA, PosRNA,units="days")/30)]							#difftime between subsequent PosRNA
						p.infectious	<- sum(p.infectious[ which(z[,lRNA.infectious]) ]) / sum(p.infectious)	#prop of time in lRNA.infectious
						lt.highb4tr		<- subset(z,lRNA.high.b4tr)
						if(nrow(lt.highb4tr))
							lt.highb4tr	<- lt.highb4tr[,StopRNA][nrow(lt.highb4tr)]	#last time before tr that VL was high
						else
							lt.highb4tr	<- as.Date(NA)
						if(!is.na(lt.highb4tr) && lt.highb4tr>AnyT_T1[1])
							lt.highb4tr	<- AnyT_T1[1]
						early			<- ifelse(	any(lRNA.high.b4tr)  && any( lRNA[seq_len( which(lRNA.high.b4tr)[1] )]<lRNA.max.b4early )	,TRUE,FALSE)
						list(lRNA.i= p.infectious, lRNA.hb4tr_LT=lt.highb4tr, lRNA.early= early )										
					}, by="Patient"]
	df<- merge(subset(df,select=c(Patient, PosRNA, RNA, lRNA, PoslRNA_T1, lRNA_T1)), subset(tmp,select=c(Patient,lRNA.i,lRNA.hb4tr_LT,lRNA.early)), all.x=1, by="Patient")
	
	file		<- paste(dir.name,"derived/ATHENA_2013_03_Viro.R",sep='/')
	if(verbose) cat(paste("\nsave to", file))
	save(df, file=file)
}	
######################################################################################
project.hivc.Excel2dataframe<- function(dir.name= DATA, min.seq.len=21, verbose=1)
{
	if(0)
	{
		#read SEQUENCE csv data file and preprocess				
		names.GeneCode	<- c("PROT","RT")
		NA.DateRes		<- as.Date("1911-11-11") 
		
		file			<- paste(dir.name,"derived/ATHENA_2013_03_Sequences.csv",sep='/')
		df				<- read.csv(file, stringsAsFactors=FALSE)				
		proc.GeneCode	<- c(1,2)
		df				<- lapply(proc.GeneCode,function(gene)
				{
					cat(paste("\nprocess GeneCode", gene))
					tmp				<- df[ df[,"GeneCode"]==gene, c("Patient","SampleCode","DateRes","Sequence"), drop=0 ]
					tmp[,"DateRes"]	<- as.Date(tmp[,"DateRes"], format="%d/%m/%Y")
					
					nok.idx<- which( tmp[,"DateRes"]==NA.DateRes )					
					if(verbose) cat(paste("\nentries with missing DateRes, n=", length(nok.idx)))					
					if(verbose) cat(paste("\nentries with missing DateRes, SampleCode", paste(tmp[nok.idx,"SampleCode"],collapse=', ')))
					tmp[nok.idx,"DateRes"]<- NA
					if(verbose) cat(paste("\nrange of DateRes is",paste(range(tmp[,"DateRes"], na.rm=1),collapse=', ')))
					cat(paste("\nfound n=", nrow(tmp)))
					seq.len			<- nchar(tmp[,"Sequence"])
					nok.idx		<- which(seq.len<min.seq.len)
					seq.ok.idx		<- which(seq.len>=min.seq.len)
					if(verbose)	cat(paste("\ndiscarding sequences with insufficient length, n=",length(nok.idx),'\n'))
					if(verbose)	cat(paste("\ndiscarding sequence with insufficient length, SampleCode",paste(tmp[nok.idx,"SampleCode"],collapse=', ')))
					tmp				<- tmp[seq.ok.idx,]
					if(verbose) cat(paste("\nfinal n=", nrow(tmp)))
					tmp
				})		
		names(df)	<- names.GeneCode
		file		<- paste(dir.name,"derived/ATHENA_2013_03_Sequences.R",sep='/')
		cat(paste("\nsave to", file))
		save(df, file=file)
		quit("no")
	}
	if(0)
	{
		project.hivc.Excel2dataframe.Regimen()						
	}
	if(0)
	{
		NA.Acute			<- c(NA,9)
		NA.CountryInfection	<- c(NA,"")
		NA.time				<- c("","01/01/1911","11/11/1911")		
		NA.transmission		<- 900
		#read PATIENT csv data file and preprocess
		file				<- paste(dir.name,"derived/ATHENA_2013_03_Patient.csv",sep='/')
		df					<- read.csv(file, stringsAsFactors=FALSE)									
		df$isDead			<- as.numeric( df[,"DateDied"]!="")
		df[which(df[,"Transmission"]==NA.transmission),"Transmission"]<- NA

		date.var			<- c("DateBorn","MyDateNeg1","MyDatePos1","DateDied","DateLastContact","DateFirstEverCDCC")		
		for(x in date.var)
		{
			cat(paste("\nprocess Time", x))
			nok.idx			<- which( df[,x]==NA.time[1] )
			if(verbose)	cat(paste("\nentries with format ",NA.time[1],", n=", length(nok.idx)))
			#if(verbose && length(nok.idx))	cat(paste("\nentries with format ",NA.time[1],", Patient", paste(df[nok.idx,"Patient"],collapse=', ')))
			
			if(length(nok.idx))	
				df[nok.idx,x]	<- NA
			nok.idx			<- which( df[,x]==NA.time[2] )
			if(verbose)	cat(paste("\nentries with format ",NA.time[2],", n=", length(nok.idx)))
			if(length(nok.idx))
				df[nok.idx,x]	<- NA
			nok.idx			<- which( df[,x]==NA.time[3] )
			if(verbose)	cat(paste("\nentries with format ",NA.time[3],", n=", length(nok.idx)))
			if(length(nok.idx))
				df[nok.idx,x]	<- NA
			df[,x]			<- as.Date(df[,x], format="%d/%m/%Y")	
		}
		file		<- paste(dir.name,"derived/ATHENA_2013_03_Patients.R",sep='/')
		if(verbose) cat(paste("\nsave to", file))
		
		df<- data.table(df)
		setnames(df, c("MyDateNeg1_Acc","MyDatePos1_Acc","AcuteInfection","Transmission","HospitalRegion"), c("NegT_Acc","PosT_Acc","isAcute","Trm","RegionHospital"))
		set(df, NULL, "Patient", factor(df[,Patient]))
		set(df, which( df[,isAcute%in%NA.Acute] ), "isAcute", NA )		
		set(df, NULL, "isAcute", factor(df[,isAcute], levels=c(0,1,2), labels=c("No","Yes","Maybe")) )
		set(df, NULL, "Sex", factor(df[,Sex], levels=c(1,2), labels=c("M","F")))
		set(df, NULL, "Subtype", factor(df[,Subtype]))
		set(df, NULL, "CountryBorn", factor(df[,CountryBorn]))
		set(df, which( df[,CountryInfection%in%NA.CountryInfection] ), "CountryInfection", NA_character_ )		
		set(df, NULL, "CountryInfection", factor(df[,CountryInfection]))
		set(df, NULL, "RegionOrigin", factor(df[,RegionOrigin]))
		set(df, NULL, "NegT_Acc", factor(df[,NegT_Acc], levels=c(0,1), labels=c("No","Yes")))
		set(df, NULL, "PosT_Acc", factor(df[,PosT_Acc], levels=c(0,1), labels=c("No","Yes")))		
		set(df, NULL, "isDead", factor(df[,isDead], levels=c(0,1), labels=c("No","Yes")))
		set(df, NULL, "Trm", factor(df[, Trm], levels=c(100, 101,  102,  202, 103,  104,  105,  106,  107, 108,  110), labels= c("MSM","BI","HET","HETfa","IDU","BLOOD","NEEACC", "PREG", "BREAST", "OTH", "SXCH")) )
		set(df, NULL, "RegionHospital", factor(df[,RegionHospital], levels=c(1,2,3,4,5,6), labels=c("Amst","N","E","S","W","Curu")))
		setkey(df,Patient)
		str(df)
		save(df, file=file)		
	}
	if(0)
	{
		project.hivc.Excel2dataframe.Viro()	
	}
	if(0)
	{
		project.hivc.Excel2dataframe.CD4()		
	}
}
######################################################################################
#create PROT+RT data set of first sequences from all patients
hivc.prog.get.bootstrapseq<- function(check.any.bs.identical=0)
{	
	library(ape)
	library(data.table)
	library(hivclust)
	
	indir				<- outdir		<- paste(DATA,"tmp",sep='/')
	infile				<- "ATHENA_2013_03_FirstCurSequences_PROTRT"
	signat.out			<- signat.in	<- "Sat_May_11_14/23/46_2013"
	verbose				<- resume		<- 1
	opt.bootstrap.by	<- "codon"
	bs					<- 0
	
	if(exists("argv"))
	{
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,6),
									indir= return(substr(arg,8,nchar(arg))),NA)	}))
		if(length(tmp)>0) indir<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,7),
									outdir= return(substr(arg,9,nchar(arg))),NA)	}))
		if(length(tmp)>0) outdir<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,7),
									infile= return(substr(arg,9,nchar(arg))),NA)	}))
		if(length(tmp)>0) infile<- tmp[1]				
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,9),
									insignat= return(substr(arg,11,nchar(arg))),NA)	}))
		if(length(tmp)>0) signat.in<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,10),
									outsignat= return(substr(arg,12,nchar(arg))),NA)	}))
		if(length(tmp)>0) signat.out<- tmp[1]		
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,7),
									resume= return(as.numeric(substr(arg,9,nchar(arg)))),NA)	}))
		if(length(tmp)>0) resume<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,2),
									v= return(as.numeric(substr(arg,4,nchar(arg)))),NA)	}))
		if(length(tmp)>0) verbose<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,10),
									bootstrap= return(as.numeric(substr(arg,12,nchar(arg)))),NA)	}))
		if(length(tmp)>0) bs<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,3),
									by= return(substr(arg,5,nchar(arg))),NA)	}))
		if(length(tmp)>0) opt.bootstrap.by<- tmp[1]
	}
	if(1)
	{
		print( indir ) 
		print(outdir)
		print(infile)
		print(signat.in)
		print(signat.out)
		print(verbose)
		print(resume)
		print(bs)
		print(opt.bootstrap.by)
		print(signat.in)
		print(signat.out)
	}
	if(!opt.bootstrap.by%in%c("nucleotide","codon"))	stop("Unexpected opt.bootstrap.by")		
	pattern 	<- paste(infile,"_",gsub('/',':',signat.out),".phylip.",sprintf("%03d",bs),sep='')
	file		<- list.files(path=outdir, pattern=pattern, full.names=1)
	if(!resume || !length(file))	
	{					
		file		<- paste(outdir,"/",infile,"_",gsub('/',':',signat.out),".R",sep='')
		if(verbose) cat(paste("\nread",file))
		tmp			<- load(file)
		if(length(tmp)!=1)	stop("Unexpected lenght of loaded objects")
		eval(parse(text=paste("seq.PROT.RT<- ",tmp,sep='')))		
		print(seq.PROT.RT)
		#print(bs)
		if(bs)		#keep bs0 intact
		{
			dummy			<- 0
			any.eq			<- 1
			j				<- 0
			while(any.eq)
			{
				j			<- j+1
				if(opt.bootstrap.by=="codon")
				{
					bs.blocks.n	<- floor( ncol(seq.PROT.RT )/3)
					bs.blocks.s	<- sample(seq_len(bs.blocks.n),bs.blocks.n,replace=T)-1
					bs.seq.s	<- as.numeric( sapply(bs.blocks.s,function(x)		3*x+c(1,2,3)		) )
				}
				else if(opt.bootstrap.by=="nucleotide")
				{
					bs.seq.s	<- sample(seq_len(ncol(seq.PROT.RT )),ncol(seq.PROT.RT ),replace=T)
				}
				seq.BS		<- seq.PROT.RT[,bs.seq.s]				
				if(check.any.bs.identical)
				{
					if(verbose) cat(paste("\ncheck for identity proposed boostrap seq alignment no",j))
					#check no seqs identical								
					for(i1 in seq_len(nrow(seq.BS)-1))
					{
						seq1		<- seq.BS[i1,]
						tmp			<- 1-sapply(seq.int(i1+1,nrow(seq.BS)),function(i2)
													{		
														.C("hivc_dist_ambiguous_dna", seq1, seq.BS[i2,], ncol(seq1), dummy )[[4]]			
													})
						#print(tmp)
						if(any(tmp==0))
						{
							print(tmp)
							break
						}									
						if(i1==nrow(seq.BS)-1)
							any.eq	<- 0
					}
					if(verbose) cat(paste("\nchecked for identity proposed boostrap seq alignment no",j,"is any identical",any.eq))
				}
				else
					any.eq	<- 0
			}					
		}
		else
		{
			cat(paste("\nkeep boostrap seq alignment no",bs,"as original"))
			seq.BS	<- seq.PROT.RT
		}
		file		<- paste(outdir,"/",infile,"_",gsub('/',':',signat.out),".phylip.",sprintf("%03d",bs),sep='')
		cat(paste("\nsave boostrap seq alignment to",file))
		hivc.seq.write.dna.phylip(seq.BS, file=file)
	}
	else
		cat("\nfound boostrap sequence alignment")
}
######################################################################################
#create PROT+RT data set of first sequences from all patients
hivc.prog.get.allseq<- function()
{	
	library(ape)
	library(data.table)
	library(RFLPtools)
	library(hivclust)
	
	indir		<- paste(DATA,"derived",sep='/')
	outdir		<- paste(DATA,"tmp",sep='/')
	infile		<- "ATHENA_2013_03_All1stPatientCovariates.R"
	signat.out	<- "Sat_Jun_16_17/23/46_2013"
	signat.in	<- "Wed_May__1_17/08/15_2013"
	verbose		<- 1 
	resume		<- 0
	
	if(exists("argv"))
	{
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,6),
									indir= return(substr(arg,8,nchar(arg))),NA)	}))
		if(length(tmp)>0) indir<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,7),
									outdir= return(substr(arg,9,nchar(arg))),NA)	}))
		if(length(tmp)>0) outdir<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,7),
									infile= return(substr(arg,9,nchar(arg))),NA)	}))
		if(length(tmp)>0) infile<- tmp[1]				
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,9),
									insignat= return(substr(arg,11,nchar(arg))),NA)	}))
		if(length(tmp)>0) signat.in<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,10),
									outsignat= return(substr(arg,12,nchar(arg))),NA)	}))
		if(length(tmp)>0) signat.out<- tmp[1]		
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,7),
									resume= return(as.numeric(substr(arg,9,nchar(arg)))),NA)	}))
		if(length(tmp)>0) resume<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,2),
									v= return(as.numeric(substr(arg,4,nchar(arg)))),NA)	}))
		if(length(tmp)>0) verbose<- tmp[1]
	}
	if(1)
	{
		print(indir)
		print(outdir)
		print(infile)
		print(signat.in)
		print(signat.out)
		print(verbose)
		print(resume)
	}
	pattern 	<- gsub('/',':',paste("ATHENA_2013_03_All+LANL_Sequences_",signat.out,".R$",sep=''))
	file		<- list.files(path=outdir, pattern=pattern, full.names=1)
	if(!resume || !length(file))	
	{						
		#create file of sequences that are in preliminary cluster. use these as seed to enrich data set to get more reliable boostrap
		dir.name		<- DATA
		infile.tree		<- "ATHENA_2013_03_FirstCurSequences_PROTRT_examlbs100"
		infile.seq		<- "ATHENA_2013_03_FirstCurSequences_PROTRT"
		signat.in		<- "Sat_May_11_14/23/46_2013"
		
		
		infile.enrich	<- "LosAlamos_HIV1B_Prot_P51_n4243_gapsyes"
		outfile.enrich	<- "LosAlamos_HIV1B_Prot_P51_n4243_gapsno"
		infile.enrich	<- "LosAlamos_HIV1B_Prot_n85928_gapsyes"
		outfile.enrich	<- "LosAlamos_HIV1B_Prot_n85928_gapsno"

		if(0) 
		{
			insignat				<- "Sat_Jun_15_18/23/46_2013"
			#build BLAST database
			#remove all gaps from FASTA file
			file<- paste(paste(dir.name,"original",sep='/'),'/',paste(infile.enrich,'_',gsub('/',':',insignat),".fasta", sep=''), sep='')						
			seq.enrich				<- read.dna( file, format="fa" )
			seq.enrich				<- hivc.seq.rmgaps(seq.enrich, rm.only.col.gaps=0)
			names(seq.enrich)		<- gsub('-','NA',names(seq.enrich))					
			file<- paste(paste(dir.name,"derived",sep='/'),'/',paste(outfile.enrich,'_',gsub('/',':',insignat),".fasta", sep=''), sep='')			
			write.dna(seq.enrich, file=file, format="fasta", append=0, colsep='', colw=80, blocksep=0)
			#build BLAST database
			cmd						<- hivc.cmd.blast.makedb(paste(dir.name,"derived",sep='/'), outfile.enrich, signat=gsub('/',':',insignat), with.mask=0, verbose=0)
			cat ( cmd )
			#cat( system(cmd, intern=TRUE) )
		}
		if(0)
		{	
			plot							<- 0
			#identify query sequences for BLAST search from large trial clusters
			signat.in						<- "Fri_May_24_12/59/06_2013"						
			#read tree to construct large trial clusters
			file							<- paste(dir.name,"tmp",paste(infile.tree,'_',gsub('/',':',signat.in),".newick",sep=''),sep='/')
			cat(paste("read file",file))
			ph								<- ladderize( read.tree(file) )
			#read bootstrap support values		
			ph.node.bs						<- as.numeric( ph$node.label )
			ph.node.bs[is.na(ph.node.bs)]	<- 0
			ph.node.bs						<- ph.node.bs/100
			ph$node.label					<- ph.node.bs
			dist.brl						<- hivc.clu.brdist.stats(ph, eval.dist.btw="leaf")			
			#produce trial clustering
			thresh.bs						<- 0.9
			thresh.brl						<- 0.105		#subst rate 3.4 10^-3  so for 15 years either way have 10.5 * 10^-2
			clustering						<- hivc.clu.clusterbythresh(ph, thresh.nodesupport=thresh.bs, thresh.brl=thresh.brl, dist.brl=dist.brl, nodesupport=ph.node.bs,retval="all")
			#get tip names in trial cluster
			ph.tips.in.cluster				<- which( !is.na(clustering[["clu.mem"]][seq_len(Ntip(ph))]) )
			ph.tips.in.cluster.names		<- ph$tip.label[ph.tips.in.cluster]
			#get sequences in trial cluster			
			file							<- paste(dir.name,"tmp",paste(infile.seq,'_',gsub('/',':',signat.in),".R",sep=''),sep='/')
			load(file)
			if(verbose) cat(paste("\nload trial sequences from",file,sep=''))
			seq.Trial.PROT.RT				<- seq.PROT.RT[ph.tips.in.cluster.names,]
			seq.Trial.PROT.RT				<- hivc.seq.replace(seq.Trial.PROT.RT, code.from='?', code.to='n')			
			if(verbose) cat(paste("\nfound trial sequences, n=",nrow(seq.Trial.PROT.RT),sep=''))			
			outsignat						<- "Sat_Jun_15_18/23/46_2013"
			file							<- paste(paste(dir.name,"tmp",sep='/'),'/',paste(infile.tree,"_intrialcluster_",gsub('/',':',outsignat),".fasta", sep=''), sep='')
			if(verbose) cat(paste("\nwrite trial sequences to",file,sep=''))
			write.dna(seq.Trial.PROT.RT, file=file, format="fasta", append=0, colsep='', colw=ncol(seq.Trial.PROT.RT), blocksep=0)
			
			#BLAST search against PROT+P51 database
			indir							<- paste(dir.name,"tmp",sep='/')
			infile							<- paste(infile.tree,"_intrialcluster",sep='')
			insignat						<- "Sat_Jun_15_18/23/46_2013"
			dbdir							<- paste(dir.name,"derived",sep='/')
			dbfile							<- "LosAlamos_HIV1B_Prot_P51_n4243_gapsno"			
			outsignat						<- "Sat_Jun_15_18/23/46_2013"
			cmd<- hivc.cmd.blast(indir, infile, gsub('/',':',insignat), dbdir, dbfile, gsub('/',':',insignat), outdir=indir, outsignat=gsub('/',':',outsignat), blast.max_target_seqs=20)			
			cat(cmd)
			system(cmd, intern=TRUE)	#wait until finished	
			#read BLAST file against PROT+P51 database and take unique best 10 hits
			file							<- paste(paste(dir.name,"tmp",sep='/'),'/',paste(infile.tree,"_intrialcluster_",gsub('/',':',signat.out),".blast", sep=''), sep='')
			seq.enrich						<- hivc.seq.blast.read(file=file)
			if(plot)
			{
				tmp							<- 1-seq.enrich[,identity]/100
				summary(tmp)
				hist(tmp)
			}
			seq.enrich.unique				<- data.table(subject.id=unique(seq.enrich[,subject.id]), key="subject.id")			
			tmp								<- seq.enrich.unique[, substr(subject.id,3,4)]
			tmp[tmp=="NA"]					<- NA
			seq.enrich.unique[,subject.location:=tmp]
			tmp								<- sapply(seq.enrich.unique[, strsplit(subject.id,'.',fixed=1) ], function(x) x[3])
			tmp[tmp=="NA"]					<- NA
			seq.enrich.unique[,subject.PosSeq:=as.numeric(tmp)]
			if(plot)
			{
				hist(seq.enrich.unique[,subject.PosSeq])
				seq.enrich.unique.loc		<- sort(table(seq.enrich.unique[,subject.location]),decreasing=TRUE)
				barplot(seq.enrich.unique.loc, cex.names=.75)
			}			
			#best PROT+p51 hits to enrich NL data set
			indir				<- paste(dir.name,"derived",sep='/')
			infile				<- "LosAlamos_HIV1B_Prot_P51_n4243_gapsno"
			insignat			<- "Sat_Jun_15_18/23/46_2013"
			file				<- paste(indir,'/',paste(infile,'_',gsub('/',':',insignat),".fasta", sep=''), sep='')
			seq.enrich.PROT.P51	<- read.dna( file, format="fa" )
			tmp					<- which( names(seq.enrich.PROT.P51)%in%seq.enrich.unique[,subject.id] )
			tmp2				<- lapply( tmp, function(i)	seq.enrich.PROT.P51[[i]] )
			names(tmp2)			<- paste("PROT+P51",names(seq.enrich.PROT.P51)[tmp],sep="_")
			class(tmp2)			<- "DNAbin"
			seq.enrich.PROT.P51	<- tmp2
			
			
			#BLAST search against PROT database
			indir							<- paste(dir.name,"tmp",sep='/')
			infile							<- paste(infile.tree,"_intrialcluster",sep='')			
			dbdir							<- paste(dir.name,"derived",sep='/')			
			dbfile							<- "LosAlamos_HIV1B_Prot_n85928_gapsno"
			outsignat						<- "Sat_Jun_16_17/23/46_2013"
			cmd<- hivc.cmd.blast(indir, infile, gsub('/',':',insignat), dbdir, dbfile, gsub('/',':',insignat), outdir=indir, outsignat=gsub('/',':',outsignat), blast.max_target_seqs=20)			
			cat(cmd)
			#cat( system(cmd, intern=TRUE) )
			
			#read BLAST file against PROT database and take unique best 10 hits that are geographically distant as TRUE NEGATIVES
			file							<- paste(paste(dir.name,"tmp",sep='/'),'/',paste(infile.tree,"_intrialcluster_",gsub('/',':',outsignat),".blast", sep=''), sep='')
			seq.enrich						<- hivc.seq.blast.read(file=file)
			if(plot)
			{
				tmp							<- 1-seq.enrich[,identity]/100
				summary(tmp)
				hist(tmp)
			}
			seq.enrich.unique				<- data.table(subject.id=unique(seq.enrich[,subject.id]), key="subject.id")			
			tmp								<- seq.enrich.unique[, substr(subject.id,3,4)]
			tmp[tmp=="NA"]					<- NA
			seq.enrich.unique[,subject.location:=tmp]
			tmp								<- sapply(seq.enrich.unique[, strsplit(subject.id,'.',fixed=1) ], function(x) x[3])
			tmp[tmp=="NA"]					<- NA
			seq.enrich.unique[,subject.PosSeq:=as.numeric(tmp)]
			if(plot)
			{
				hist(seq.enrich.unique[,subject.PosSeq])
				seq.enrich.unique.loc			<- sort(table(seq.enrich.unique[,subject.location]),decreasing=TRUE)
				barplot(seq.enrich.unique.loc, cex.names=.75)
			}
			
			#to build HIVC.COUNTRY.TABLE, use this to search LANL
			#tmp								<- seq.enrich.unique[,min(subject.id),by=subject.location]
			#tmp								<- sapply(tmp[, strsplit(V1,'.',fixed=1) ], function(x) x[5])			
			#extract seq.enrich.unique that are very unlikely to transmit to NL		seq.enrich.unique.loc.excl taken from Bezemer2013
			seq.enrich.unique.loc.level		<- "0.01"			
			seq.enrich.unique.loc.excl		<- structure(list(	`0.01` = c("GREENLAND", "CHILE", "ALBANIA", "Suriname", "SINGAPORE", "PANAMA", "SENEGAL", "FINLAND", "TAIWAN", "JAPAN", 
																			"SLOVAKIA", "LUXEMBOURG", "AUSTRALIA", "THAILAND", "SERBIA", "ROMANIA", "NORWAY", "VENEZUELA", "SLOVENIA", "RUSSIAN FEDERATION", 
																			"AUSTRIA", "SWEDEN", "FRANCE", "SOUTH KOREA", "CYPRUS", "MONTENEGRO", "CUBA", "HONDURAS", "DENMARK", "POLAND", "PORTUGAL", "CHINA", 
																			"SWITZERLAND", "BRAZIL", "ARGENTINA", "GERMANY", "BELGIUM", "CANADA","CZECH REPUBLIC", "UNITED STATES", "SPAIN", "ITALY", "UNITED KINGDOM","NETHERLANDS"), 
																`0.05` = c("SERBIA", "ROMANIA", "NORWAY", "VENEZUELA", "SLOVENIA","RUSSIAN FEDERATION", "AUSTRIA", "SWEDEN", "FRANCE", "SOUTH KOREA", 
																			"CYPRUS", "MONTENEGRO", "CUBA", "HONDURAS", "DENMARK", "POLAND","PORTUGAL", "CHINA", "SWITZERLAND", "BRAZIL", "ARGENTINA", "GERMANY",
																			"BELGIUM", "CANADA", "CZECH REPUBLIC", "UNITED STATES", "SPAIN","ITALY", "UNITED KINGDOM","NETHERLANDS")), .Names = c("0.01", "0.05"))																										
			seq.enrich.unique.loc.incl		<- subset(HIVC.COUNTRY.TABLE, !country%in%seq.enrich.unique.loc.excl[[seq.enrich.unique.loc.level]] )
			seq.enrich.unique.loc.incl		<- subset(seq.enrich.unique.loc.incl, !country%in%c("BULGARIA","GRENADA","GIORGIA","GREECE","HUNGARY","UKRAINE","LATVIA") )						
			seq.enrich.unlikelytransmission	<- subset(seq.enrich.unique,subject.location%in%as.character(seq.enrich.unique.loc.incl[,key]) )
			if(plot)
			{				
				hist(seq.enrich.unlikelytransmission[,subject.PosSeq])
				seq.enrich.unique.loc			<- sort(table(seq.enrich.unlikelytransmission[,subject.location]),decreasing=TRUE)
				barplot(seq.enrich.unique.loc, cex.names=.75)
			}
			#best geographically distant PROT hits to enrich NL data set
			indir				<- paste(dir.name,"derived",sep='/')
			infile				<- "LosAlamos_HIV1B_Prot_n85928_gapsno"
			insignat			<- "Sat_Jun_15_18/23/46_2013"
			file				<- paste(indir,'/',paste(infile,'_',gsub('/',':',insignat),".fasta", sep=''), sep='')
			seq.enrich.TN		<- read.dna( file, format="fa" )
			tmp					<- which( names(seq.enrich.TN)%in%seq.enrich.unlikelytransmission[,subject.id] )
			tmp2				<- lapply( tmp, function(i)	seq.enrich.TN[[i]] )
			names(tmp2)			<- paste("TN",names(seq.enrich.TN)[tmp],sep="_")
			class(tmp2)			<- "DNAbin"
			seq.enrich.TN		<- tmp2

			#combine all sequence data sets and save
			seq.enrich			<- c(seq.enrich.PROT.P51,seq.enrich.TN)
			outdir				<- paste(dir.name,"derived",sep='/')
			outsignat			<- "Sat_Jun_16_17/23/46_2013"			
			file				<- paste(outdir,"/LosAlamos_EnrichSequences_For_ATHENA201303_",gsub('/',':',outsignat),".R",sep='')
			if(verbose) cat(paste("\nwrite to",file))
			save(seq.enrich, file=file)						
		}
		if(0)
		{
			#get all ATHENA sequences into PROT+RT format, and add foreign sequences
			indir		<- paste(dir.name,"tmp",sep='/')
			outdir		<- paste(dir.name,"tmp",sep='/')
			insignat	<- "Wed_May__1_17/08/15_2013"
			outsignat	<- "Sat_Jun_16_17/23/46_2013"
			outfile		<- "ATHENA_2013_03_All+LANL_Sequences"
			
			if(verbose)	cat(paste("\ncreate LosAlamos_EnrichSequences_For_ATHENA201303 file"))						
			#create matrix of all PROT+RT sequences				
			pattern 	<- gsub('/',':',paste(insignat,".clustalo$",sep=''))
			files		<- list.files(path=indir, pattern=pattern, full.names=1)
			#read all sequences and add a missing one with name "NA" if not in union of all possible sequences take at a sampling date		
			seq.PROT	<- read.dna( files[ grep("PROT",files) ], format="fa", as.matrix=1 )
			if(verbose) cat(paste("\nfound PROT sequences, n=",nrow(seq.PROT),"\n"))
			seq.RT		<- read.dna( files[ grep("RT",files) ], format="fa", as.matrix=1 )
			if(verbose) cat(paste("\nfound RT sequences, n=",nrow(seq.RT),"\n"))
			seq.nam.all	<- union( rownames(seq.PROT), rownames(seq.RT) )
			if(verbose) cat(paste("\nfound unique sampling IDs, n=",length(seq.nam.all),"\n"))
			seq.nam.PROT<- seq.nam.all
			seq.nam.RT	<- seq.nam.all
			seq.nam.PROT[ !seq.nam.all%in%rownames(seq.PROT) ]	<- "NA"		#those sampling dates missing among PROT get name NA
			seq.nam.RT[ !seq.nam.all%in%rownames(seq.RT) ]		<- "NA"
			#prepare PROT and RT sequences: add a missing one with name "NA" 		
			seq.PROT	<- read.dna( files[ grep("PROT",files) ], format="fa", as.matrix=1 )
			tmp			<- as.DNAbin( matrix(rep('-',ncol(seq.PROT)),1,ncol(seq.PROT), dimnames=list(c("NA"),c())) )
			seq.PROT	<- rbind(seq.PROT,tmp)
			seq.RT		<- read.dna( files[ grep("RT",files) ], format="fa", as.matrix=1 )						 				
			tmp			<- as.DNAbin( matrix(rep('-',ncol(seq.RT)),1,ncol(seq.RT), dimnames=list(c("NA"),c())) )
			seq.RT		<- rbind(seq.RT,tmp)
			#add missing sequence where needed to combine PROT and RT 
			seq.PROT	<- seq.PROT[seq.nam.PROT,]
			seq.RT		<- seq.RT[seq.nam.RT,]
			rownames(seq.PROT)	<- seq.nam.all
			rownames(seq.RT)	<- seq.nam.all
			#combine PROT and RT
			seq.PROT.RT	<- cbind(seq.PROT,seq.RT)
			if(verbose) cat(paste("\ncombined PROT and RT sequences, n=",nrow(seq.PROT.RT),"\n"))
			print(seq.PROT.RT)
			
			#load EnrichSequences
			indir				<- paste(dir.name,"derived",sep='/')
			insignat			<- "Sat_Jun_16_17/23/46_2013"
			file				<- paste(indir,"/LosAlamos_EnrichSequences_For_ATHENA201303_",gsub('/',':',insignat),".R",sep='')
			if(verbose) cat(paste("\nloading file",file))
			load(file)	
			print(seq.enrich)
			#load HXB2 reference sequence
			data( refseq_hiv1_hxb2 )
			hxb2				<- as.character( data.table( hxb2 )[, HXB2.K03455 ] )
			hxb2				<- hxb2[seq_len(length(hxb2)-2)]
			
			#write all sequences to fasta file
			file		<- paste(outdir,'/',outfile,'_',gsub('/',':',outsignat),".fasta",sep='')			
			if(verbose) cat(paste("\nwrite reference sequence to",file,"\n"))
			cat( paste(">HXB2\n",paste( hxb2, collapse='',sep='' ),"\n",sep=''), file=file, append=0)			
			if(verbose) cat(paste("\nappend ATHENA combined PROT and RT sequences to",file,"\n"))
			write.dna(seq.PROT.RT, file=file, format="fasta", append=1, colsep='', colw=length(hxb2), blocksep=0)
			if(verbose) cat(paste("\nappend LosAlamos_EnrichSequences to",file,"\n"))
			write.dna(seq.enrich, file=file, format="fasta", append=1, colsep='', colw=length(hxb2), blocksep=0)				 
		}		
		if(1)	#curate alignment with reference 
		{
			indir								<- paste(dir.name,"tmp",sep='/')
			infile								<- "ATHENA_2013_03_All+LANL_Sequences"
			outdir								<- paste(dir.name,"tmp",sep='/')
			outfile								<- "ATHENA_2013_03_CurAll+LANL_Sequences"
			insignat							<- "Sat_Jun_16_17/23/46_2013"
			outsignat							<- "Sat_Jun_16_17/23/46_2013"
			file								<- paste(indir,'/',infile,'_',gsub('/',':',insignat),".fasta.clustalo",sep='')			
			if(verbose) cat(paste("\nread ",file))
			seq.PROT.RT							<- read.dna(file, format="fasta", as.character=1)
			query.yes							<- hivc.seq.find(seq.PROT.RT, 2273, c("g","-","-"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,2273:2275]<- matrix( c("-","-","g"), nrow=length(query.yes), ncol=3, byrow=1 )
			query.yes							<- hivc.seq.find(seq.PROT.RT, 2348, c("?","g"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,2348:2349]<- matrix( c("g","?"), nrow=length(query.yes), ncol=2, byrow=1 )			
			query.yes							<- hivc.seq.find(seq.PROT.RT, 2350, c("-","-","-","t"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,2350:2353]<- matrix( c("t","-","-","-"), nrow=length(query.yes), ncol=4, byrow=1 )
			#always delete drug resistance from pos 2356 because this is really cryptic; leave g at 2355 always OK
			seq.PROT.RT[,2356:2363]				<- matrix( c("-","-","-","-","-","-","-","-"), nrow=nrow(seq.PROT.RT), ncol=8, byrow=1 )			
			write.dna(seq.PROT.RT, file=paste(indir,'/',infile,'_',gsub('/',':',insignat),".fasta.clustalo1",sep=''), format="fasta", colsep='', colw=ncol(seq.PROT.RT), blocksep=0)			
			#manual edits seq 3151 3243 9049
			query.yes							<- hivc.seq.find(seq.PROT.RT, 2624, c("t","-"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,2624:2625]<- matrix( c("-","t"), nrow=length(query.yes), ncol=2, byrow=1 )			
			query.yes							<- hivc.seq.find(seq.PROT.RT, 2634, c("-","g"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,2634:2635]<- matrix( c("g","-"), nrow=length(query.yes), ncol=2, byrow=1 )			
			seq.PROT.RT.sort.by					<- apply(seq.PROT.RT,1,function(x)		which(rev(x)!="-" )[1]  )
			seq.PROT.RT							<- seq.PROT.RT[sort(seq.PROT.RT.sort.by, index.return=1)$ix,]
			write.dna(seq.PROT.RT, file=paste(indir,'/',infile,'_',gsub('/',':',insignat),".fasta.clustalo2",sep=''), format="fasta", colsep='', colw=ncol(seq.PROT.RT), blocksep=0)
			
			query.yes							<- hivc.seq.find(seq.PROT.RT, 2759, c("-","g","a"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,2759:2761]<- matrix( c("g","a","-"), nrow=length(query.yes), ncol=3, byrow=1 )
			query.yes							<- hivc.seq.find(seq.PROT.RT, 2759, c("-","a","a"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,2759:2761]<- matrix( c("a","a","-"), nrow=length(query.yes), ncol=3, byrow=1 )
			query.yes							<- hivc.seq.find(seq.PROT.RT, 2760, c("-","a"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,2760:2761]<- matrix( c("a","-"), nrow=length(query.yes), ncol=2, byrow=1 )			
			query.yes							<- hivc.seq.find(seq.PROT.RT, 2760, c("-","g"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,2760:2761]<- matrix( c("g","-"), nrow=length(query.yes), ncol=2, byrow=1 )			
			query.yes							<- hivc.seq.find(seq.PROT.RT, 2760, c("c","-"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,2760:2761]<- matrix( c("-","c"), nrow=length(query.yes), ncol=2, byrow=1 )			
			query.yes							<- hivc.seq.find(seq.PROT.RT, 2761, c("-","c"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,2761:2762]<- matrix( c("c","-"), nrow=length(query.yes), ncol=2, byrow=1 )			
			#always delete drug resistance from pos 2752 because this is really cryptic; leave a at 2751 always OK
			seq.PROT.RT[,2752:2759]				<- matrix( c("-","-","-","-","-","-","-","-"), nrow=nrow(seq.PROT.RT), ncol=8, byrow=1 )
			write.dna(seq.PROT.RT, file=paste(indir,'/',infile,'_',gsub('/',':',insignat),".fasta.clustalo3",sep=''), format="fasta", colsep='', colw=ncol(seq.PROT.RT), blocksep=0)
			#double fixup needed
			query.yes							<- hivc.seq.find(seq.PROT.RT, 2759, c("-","g","a"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,2759:2761]<- matrix( c("g","a","-"), nrow=length(query.yes), ncol=3, byrow=1 )
			query.yes							<- hivc.seq.find(seq.PROT.RT, 2759, c("-","g","g"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,2759:2761]<- matrix( c("g","g","-"), nrow=length(query.yes), ncol=3, byrow=1 )
			query.yes							<- hivc.seq.find(seq.PROT.RT, 2759, c("-","a","a"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,2759:2761]<- matrix( c("a","a","-"), nrow=length(query.yes), ncol=3, byrow=1 )
			query.yes							<- hivc.seq.find(seq.PROT.RT, 2759, c("-","r","a"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,2759:2761]<- matrix( c("r","a","-"), nrow=length(query.yes), ncol=3, byrow=1 )
			seq.PROT.RT[,2752:2760]				<- matrix( c("-","-","-","-","-","-","-","-","-"), nrow=nrow(seq.PROT.RT), ncol=9, byrow=1 )
			query.yes							<- hivc.seq.find(seq.PROT.RT, 2760, c("-","g"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,2760:2761]<- matrix( c("g","-"), nrow=length(query.yes), ncol=2, byrow=1 )			
			query.yes							<- hivc.seq.find(seq.PROT.RT, 2760, c("-","a"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,2760:2761]<- matrix( c("a","-"), nrow=length(query.yes), ncol=2, byrow=1 )			
			seq.PROT.RT[,2752:2760]				<- matrix( c("-","-","-","-","-","-","-","-","-"), nrow=nrow(seq.PROT.RT), ncol=9, byrow=1 )
			write.dna(seq.PROT.RT, file=paste(indir,'/',infile,'_',gsub('/',':',insignat),".fasta.clustalo4",sep=''), format="fasta", colsep='', colw=ncol(seq.PROT.RT), blocksep=0)
			query.yes							<- hivc.seq.find(seq.PROT.RT, 3149, c("-"))
			if(length(query.yes))
			{	
				seq.PROT.RT[query.yes,3149:3156]<- seq.PROT.RT[query.yes,3150:3157]
				seq.PROT.RT[query.yes,3157]		<- "-"
			}	
			write.dna(seq.PROT.RT, file=paste(indir,'/',infile,'_',gsub('/',':',insignat),".fasta.clustalo5",sep=''), format="fasta", colsep='', colw=ncol(seq.PROT.RT), blocksep=0)				
			query.yes							<- hivc.seq.find(seq.PROT.RT, 3156, c("-","a"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,3156:3157]<- matrix( c("a","-"), nrow=length(query.yes), ncol=2, byrow=1 )			
			query.yes							<- hivc.seq.find(seq.PROT.RT, 3157, c("a","-"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,3157:3158]<- matrix( c("-","a"), nrow=length(query.yes), ncol=2, byrow=1 )
			#rm col 3213 - not in HB2 and we keep the standard numbering
			seq.PROT.RT							<- seq.PROT.RT[,-3213]
			write.dna(seq.PROT.RT, file=paste(indir,'/',infile,'_',gsub('/',':',insignat),".fasta.clustalo6",sep=''), format="fasta", colsep='', colw=ncol(seq.PROT.RT), blocksep=0)
			query.yes							<- hivc.seq.find(seq.PROT.RT, 3243, c("-"))
			if(length(query.yes))
			{	
				seq.PROT.RT[query.yes,3241:3243]<- seq.PROT.RT[query.yes,3240:3242]
				seq.PROT.RT[query.yes,3240]		<- "-"
			}	
			#fixup: rm col 2671 to 2676 - not in HB2 and we return to the standard numbering
			seq.PROT.RT							<- seq.PROT.RT[,c(1:2670,2677:ncol(seq.PROT.RT))]
			write.dna(seq.PROT.RT, file=paste(indir,'/',infile,'_',gsub('/',':',insignat),".fasta.clustalo7",sep=''), format="fasta", colsep='', colw=ncol(seq.PROT.RT), blocksep=0)						
			query.yes							<- hivc.seq.find(seq.PROT.RT, 3433, c("-","t"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,3433:3434]<- matrix( c("t","-"), nrow=length(query.yes), ncol=2, byrow=1 )			
			query.yes							<- hivc.seq.find(seq.PROT.RT, 3433, c("-","-","-","-","t"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,3433:3437]<- matrix( c("t","-","-","-","-"), nrow=length(query.yes), ncol=5, byrow=1 )
			query.yes							<- hivc.seq.find(seq.PROT.RT, 3434, c("-","-","-","-","a"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,3434:3438]<- matrix( c("a","-","-","-","-"), nrow=length(query.yes), ncol=5, byrow=1 )
			query.yes							<- hivc.seq.find(seq.PROT.RT, 3434, c("g","-","-","-","-"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,3434:3438]<- matrix( c("-","-","-","-","g"), nrow=length(query.yes), ncol=5, byrow=1 )
			write.dna(seq.PROT.RT, file=paste(indir,'/',infile,'_',gsub('/',':',insignat),".fasta.clustalo8",sep=''), format="fasta", colsep='', colw=ncol(seq.PROT.RT), blocksep=0)
			query.yes							<- hivc.seq.find(seq.PROT.RT, 3434, c("-","-","-","-","g","g"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,3434:3439]<- matrix( c("g","-","-","-","-","g"), nrow=length(query.yes), ncol=6, byrow=1 )						
			query.yes							<- hivc.seq.find(seq.PROT.RT, 3435, c("-","-","-","-","a"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,3435:3439]<- matrix( c("a","-","-","-","-"), nrow=length(query.yes), ncol=5, byrow=1 )
			query.yes							<- hivc.seq.find(seq.PROT.RT, 3435, c("-","-","-","-","r"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,3435:3439]<- matrix( c("r","-","-","-","-"), nrow=length(query.yes), ncol=5, byrow=1 )			
			query.yes							<- hivc.seq.find(seq.PROT.RT, 3436, c("-","-","-","c"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,3436:3439]<- matrix( c("c","-","-","-"), nrow=length(query.yes), ncol=4, byrow=1 )
			query.yes							<- hivc.seq.find(seq.PROT.RT, 3436, c("-","-","-","t"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,3436:3439]<- matrix( c("t","-","-","-"), nrow=length(query.yes), ncol=4, byrow=1 )
			query.yes							<- hivc.seq.find(seq.PROT.RT, 3436, c("-","-","c"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,3436:3438]<- matrix( c("c","-","-"), nrow=length(query.yes), ncol=3, byrow=1 )
			query.yes							<- hivc.seq.find(seq.PROT.RT, 3436, c("-","-","t"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,3436:3438]<- matrix( c("t","-","-"), nrow=length(query.yes), ncol=3, byrow=1 )
			query.yes							<- hivc.seq.find(seq.PROT.RT, 3436, c("-","-","y"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,3436:3438]<- matrix( c("y","-","-"), nrow=length(query.yes), ncol=3, byrow=1 )			
			query.yes							<- hivc.seq.find(seq.PROT.RT, 3437, c("-","a"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,3437:3438]<- matrix( c("a","-"), nrow=length(query.yes), ncol=2, byrow=1 )
			query.yes							<- hivc.seq.find(seq.PROT.RT, 3438, c("-","g"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,3438:3439]<- matrix( c("g","-"), nrow=length(query.yes), ncol=2, byrow=1 )
			query.yes							<- hivc.seq.find(seq.PROT.RT, 3438, c("a","-"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,3437:3438]<- matrix( c("-","a"), nrow=length(query.yes), ncol=2, byrow=1 )
			
			write.dna(seq.PROT.RT, file=paste(indir,'/',infile,'_',gsub('/',':',insignat),".fasta.clustalo9",sep=''), format="fasta", colsep='', colw=ncol(seq.PROT.RT), blocksep=0)
			query.yes							<- hivc.seq.find(seq.PROT.RT, 3491, c("c","a","-"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,3491:3493]<- matrix( c("-","c","a"), nrow=length(query.yes), ncol=3, byrow=1 )															
			query.yes							<- hivc.seq.find(seq.PROT.RT, 3501, c("-","t","a"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,3501:3503]<- matrix( c("t","a","-"), nrow=length(query.yes), ncol=3, byrow=1 )
			query.yes							<- hivc.seq.find(seq.PROT.RT, 3501, c("a","-"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,3501:3502]<- matrix( c("-","a"), nrow=length(query.yes), ncol=2, byrow=1 )			
			query.yes							<- hivc.seq.find(seq.PROT.RT, 3503, c("-"))
			for(i in query.yes)
				seq.PROT.RT[i,3501:3503]		<- c("-",seq.PROT.RT[i,3501:3502])			
			query.yes							<- hivc.seq.find(seq.PROT.RT, 3504, c("-","-","-","t"))	
			for(i in query.yes)
				seq.PROT.RT[i,3504:3524]		<- c(seq.PROT.RT[i,3507:3524],"-","-","-")
			query.yes							<- hivc.seq.find(seq.PROT.RT, 3506, c("-","-","-","t","g"))
			for(i in query.yes)
				seq.PROT.RT[i,3506:3526]		<- c(seq.PROT.RT[i,3509:3526],"-","-","-")
			query.yes							<- hivc.seq.find(seq.PROT.RT, 3507, c("-","-","-","g","a"))
			for(i in query.yes)
				seq.PROT.RT[i,3507:3537]		<-	c(seq.PROT.RT[i,3510:3537],"-","-","-")			
			query.yes							<- hivc.seq.find(seq.PROT.RT, 3543, c("-","c"))
			if(length(query.yes))
			{	
				seq.PROT.RT[query.yes,3543:3550]<- seq.PROT.RT[query.yes,3544:3551]
				seq.PROT.RT[query.yes,3551]		<- "-"
			}
			query.yes							<- hivc.seq.find(seq.PROT.RT, 3549, c("-","c","a"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,3549:3551]<- matrix( c("c","a","-"), nrow=length(query.yes), ncol=3, byrow=1 )
			query.yes							<- hivc.seq.find(seq.PROT.RT, 3549, c("-","g","a"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,3549:3551]<- matrix( c("g","a","-"), nrow=length(query.yes), ncol=3, byrow=1 )
			query.yes							<- hivc.seq.find(seq.PROT.RT, 3549, c("-","g","t"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,3549:3551]<- matrix( c("g","t","-"), nrow=length(query.yes), ncol=3, byrow=1 )
			query.yes							<- hivc.seq.find(seq.PROT.RT, 3551, c("g","g","c","-"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,3551:3554]<- matrix( c("-","g","g","c"), nrow=length(query.yes), ncol=4, byrow=1 )
			query.yes							<- hivc.seq.find(seq.PROT.RT, 3495, c("g","-","-","-","g"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,3495:3499]<- matrix( c("g","g","-","-","-"), nrow=length(query.yes), ncol=5, byrow=1 )
			write.dna(seq.PROT.RT, file=paste(indir,'/',infile,'_',gsub('/',':',insignat),".fasta.clustalo10",sep=''), format="fasta", colsep='', colw=ncol(seq.PROT.RT), blocksep=0)
			
			query.yes							<- hivc.seq.find(seq.PROT.RT, 3497, c("-","-","-","a"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,3497:3500]<- matrix( c("a","-","-","-"), nrow=length(query.yes), ncol=4, byrow=1 )						
			query.yes							<- hivc.seq.find(seq.PROT.RT, 3497, c("-","-","-","g"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,3497:3500]<- matrix( c("g","-","-","-"), nrow=length(query.yes), ncol=4, byrow=1 )									
			query.yes							<- hivc.seq.find(seq.PROT.RT, 3499, c("t","-","-"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,3499:3501]<- matrix( c("-","-","t"), nrow=length(query.yes), ncol=3, byrow=1 )			
			query.yes							<- hivc.seq.find(seq.PROT.RT, 3578, c("g","a","g","-"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,3578:3581]<- matrix( c("-","g","a","g"), nrow=length(query.yes), ncol=4, byrow=1 )
			query.yes							<- hivc.seq.find(seq.PROT.RT, 3578, c("g","-"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,3578:3579]<- matrix( c("-","g"), nrow=length(query.yes), ncol=2, byrow=1 )
			query.yes							<- hivc.seq.find(seq.PROT.RT, 3551, c("g","g"))
			for(i in query.yes)
				seq.PROT.RT[i,3551:3600]		<-	c("-",seq.PROT.RT[i,3551:3599])
			query.yes							<- hivc.seq.find(seq.PROT.RT, 3552, c("g","c"))
			for(i in query.yes)
				seq.PROT.RT[i,3552:3600]		<-	c("-",seq.PROT.RT[i,3552:3599])
			query.yes							<- hivc.seq.find(seq.PROT.RT, 3498, c("c","g","t"))
			for(i in query.yes)
				seq.PROT.RT[i,3498:3520]		<-	c("-",seq.PROT.RT[i,3498:3519])
			query.yes							<- hivc.seq.find(seq.PROT.RT, 3549, c("a","a","-"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,3549:3551]<- matrix( c("-","-","-"), nrow=length(query.yes), ncol=3, byrow=1 )
			query.yes							<- hivc.seq.find(seq.PROT.RT, 3551, c("g","a","c"))
			for(i in query.yes)
				seq.PROT.RT[i,3551:3560]		<-	c("-",seq.PROT.RT[i,3551:3559])
			query.yes							<- hivc.seq.find(seq.PROT.RT, 3551, c("g","a","y","-"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,3551:3554]<- matrix( c("-","g","a","y"), nrow=length(query.yes), ncol=4, byrow=1 )
			query.yes							<- hivc.seq.find(seq.PROT.RT, 3549, c("-","c","c","g","g"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,3549:3551]<- matrix( c("c","-","c"), nrow=length(query.yes), ncol=3, byrow=1 )
			write.dna(seq.PROT.RT, file=paste(indir,'/',infile,'_',gsub('/',':',insignat),".fasta.clustalo11",sep=''), format="fasta", colsep='', colw=ncol(seq.PROT.RT), blocksep=0)
			
			#manual curation on the end
			query.yes							<- hivc.seq.find(seq.PROT.RT, 3422, c("t","-","-"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,3422:3424]<- matrix( c("-","-","t"), nrow=length(query.yes), ncol=3, byrow=1 )
			query.yes							<- hivc.seq.find(seq.PROT.RT, 3494, c("g","g","-"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,3494:3496]<- matrix( c("-","g","g"), nrow=length(query.yes), ncol=3, byrow=1 )
			query.yes							<- hivc.seq.find(seq.PROT.RT, 3498, c("-","-","g"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,3498:3500]<- matrix( c("g","-","-"), nrow=length(query.yes), ncol=3, byrow=1 )
			query.yes							<- hivc.seq.find(seq.PROT.RT, 3578, c("g","a","-"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,3578:3580]<- matrix( c("-","g","a"), nrow=length(query.yes), ncol=3, byrow=1 )
			query.yes							<- hivc.seq.find(seq.PROT.RT, 3618, c("g","-","-","-","-"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,3618:3622]<- matrix( c("-","-","-","-","g"), nrow=length(query.yes), ncol=5, byrow=1 )
			query.yes							<- hivc.seq.find(seq.PROT.RT, 3499, c("-","c"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,3499:3500]<- matrix( c("c","-"), nrow=length(query.yes), ncol=2, byrow=1 )
			query.yes							<- hivc.seq.find(seq.PROT.RT, 3504, c("a","-"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,3504:3505]<- matrix( c("-","a"), nrow=length(query.yes), ncol=2, byrow=1 )
			query.yes							<- hivc.seq.find(seq.PROT.RT, 3437, c("-","a"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,3437:3438]<- matrix( c("a","-"), nrow=length(query.yes), ncol=2, byrow=1 )
			query.yes							<- hivc.seq.find(seq.PROT.RT, 3438, c("c","-"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,3438:3439]<- matrix( c("-","c"), nrow=length(query.yes), ncol=2, byrow=1 )
			query.yes							<- hivc.seq.find(seq.PROT.RT, 3444, c("c","-"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,3444:3445]<- matrix( c("-","c"), nrow=length(query.yes), ncol=2, byrow=1 )
			#cut before PROT pol start and at max len in database
			seq.PROT.RT							<- seq.PROT.RT[,2253:ncol(seq.PROT.RT)] 
			seq.PROT.RT							<- seq.PROT.RT[,1:1624]
			seq.PROT.RT.sort.by					<- apply(seq.PROT.RT,1,function(x)		which(rev(x)!="-" )[1]  )
			seq.PROT.RT							<- seq.PROT.RT[sort(seq.PROT.RT.sort.by, index.return=1)$ix,]
			#always delete drug resistance at G333D/E because this is at the end of the alignment and alignment unreliable - could have picked up other ends
			seq.PROT.RT[,1294:1299]				<- matrix( c("-","-","-","-","-","-"), nrow=nrow(seq.PROT.RT), ncol=6, byrow=1 )			
			write.dna(seq.PROT.RT, file=paste(indir,'/',infile,'_',gsub('/',':',insignat),".fasta.clustalo12",sep=''), format="fasta", colsep='', colw=ncol(seq.PROT.RT), blocksep=0)
			query.yes							<- hivc.seq.find(seq.PROT.RT, 1300, c("g","g","g"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,1300:1302]<- matrix( c("g","g","-"), nrow=length(query.yes), ncol=3, byrow=1 )
			query.yes							<- hivc.seq.find(seq.PROT.RT, 1300, c("g","g","t"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,1300:1302]<- matrix( c("g","g","-"), nrow=length(query.yes), ncol=3, byrow=1 )
			query.yes							<- hivc.seq.find(seq.PROT.RT, 1300, c("g","g","a"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,1300:1302]<- matrix( c("g","g","-"), nrow=length(query.yes), ncol=3, byrow=1 )
			query.yes							<- hivc.seq.find(seq.PROT.RT, 1300, c("g","a","c"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,1300:1302]<- matrix( c("-","-","-"), nrow=length(query.yes), ncol=3, byrow=1 )
			query.yes							<- hivc.seq.find(seq.PROT.RT, 1300, c("g","g","-","c"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,1300:1303]<- matrix( c("g","g","g","c"), nrow=length(query.yes), ncol=4, byrow=1 )			
			write.dna(seq.PROT.RT, file=paste(indir,'/',infile,'_',gsub('/',':',insignat),".fasta.clustalo14",sep=''), format="fasta", colsep='', colw=ncol(seq.PROT.RT), blocksep=0)						
			query.yes							<- hivc.seq.find(seq.PROT.RT, 1300, c("a","a","c","-"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,1300:1303]<- matrix( c("-","-","-","-"), nrow=length(query.yes), ncol=4, byrow=1 )
			query.yes							<- hivc.seq.find(seq.PROT.RT, 1300, c("g","a","t","-"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,1300:1303]<- matrix( c("g","-","-","-"), nrow=length(query.yes), ncol=4, byrow=1 )
			query.yes							<- hivc.seq.find(seq.PROT.RT, 1300, c("t","g","t","-"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,1300:1303]<- matrix( c("t","g","-","-"), nrow=length(query.yes), ncol=4, byrow=1 )
			query.yes							<- hivc.seq.find(seq.PROT.RT, 1300, c("a","t","c","-"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,1300:1303]<- matrix( c("-","-","-","-"), nrow=length(query.yes), ncol=4, byrow=1 )
			query.yes							<- hivc.seq.find(seq.PROT.RT, 1300, c("-","g","t","-"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,1300:1303]<- matrix( c("-","g","-","-"), nrow=length(query.yes), ncol=4, byrow=1 )
			query.yes							<- hivc.seq.find(seq.PROT.RT, 1300, c("g","a","y","-"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,1300:1303]<- matrix( c("g","-","-","-"), nrow=length(query.yes), ncol=4, byrow=1 )
			
			seq.PROT.RT							<- as.DNAbin( seq.PROT.RT )
			seq.PROT.RT							<- hivc.seq.replace(seq.PROT.RT, code.from='?', code.to='n')
			seq.PROT.RT							<- seq.PROT.RT[-1,]												#remove HXB2
			rownames(seq.PROT.RT)				<- gsub(' ','',rownames(seq.PROT.RT))
			
			file								<- paste(outdir,'/',outfile,'_',gsub('/',':',outsignat),".phylip",sep='')
			if(verbose)	cat(paste("\nwrite phylip file to",file))
			hivc.seq.write.dna.phylip(seq.PROT.RT, file=file)			
			file								<- paste(outdir,'/',outfile,'_',gsub('/',':',outsignat),".R",sep='')
			if(verbose)	cat(paste("\nwrite R file to",file))
			save(seq.PROT.RT, file=file)
			file								<- paste(outdir,'/',outfile,'_',gsub('/',':',outsignat),".fasta",sep='')			
			if(verbose)	cat(paste("\nwrite fasta file to",file))
			write.dna(seq.PROT.RT, file=file, format="fasta", colsep='', colw=ncol(seq.PROT.RT), blocksep=0)						
		}
		quit("no")				
	}
	else
	{
		if(verbose)	cat(paste("\nload",file))
		load(file)
	}	
}
######################################################################################
#create PROT+RT data set of first sequences from all patients
hivc.prog.get.firstseq<- function()
{	
	library(ape)
	library(data.table)
	library(hivclust)
	
	indir		<- outdir		<- paste(DATA,"tmp",sep='/')
	infile		<- "ATHENA_2013_03_SeqMaster.R"
	signat.out	<- signat.in	<- "Wed_May__1_17/08/15_2013"
	verbose		<- resume		<- 1
	
	if(exists("argv"))
	{
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,6),
									indir= return(substr(arg,8,nchar(arg))),NA)	}))
		if(length(tmp)>0) indir<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,7),
									outdir= return(substr(arg,9,nchar(arg))),NA)	}))
		if(length(tmp)>0) outdir<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,7),
									infile= return(substr(arg,9,nchar(arg))),NA)	}))
		if(length(tmp)>0) infile<- tmp[1]				
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,9),
									insignat= return(substr(arg,11,nchar(arg))),NA)	}))
		if(length(tmp)>0) signat.in<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,10),
									outsignat= return(substr(arg,12,nchar(arg))),NA)	}))
		if(length(tmp)>0) signat.out<- tmp[1]		
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,7),
									resume= return(as.numeric(substr(arg,9,nchar(arg)))),NA)	}))
		if(length(tmp)>0) resume<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,2),
									v= return(as.numeric(substr(arg,4,nchar(arg)))),NA)	}))
		if(length(tmp)>0) verbose<- tmp[1]
	}
	if(1)
	{
		print(indir)
		print(outdir)
		print(infile)
		print(signat.in)
		print(signat.out)
		print(verbose)
		print(resume)
	}
	pattern 	<- gsub('/',':',paste("FirstAliSequences_PROTRT_",signat.in,".R$",sep=''))
	file		<- list.files(path=outdir, pattern=pattern, full.names=1)
	if(!resume || !length(file))	
	{		
		if(verbose)	cat(paste("\ncreate FirstAliSequences_PROTRT file"))			
		file		<- paste(indir,infile,sep='/')
		if(verbose)	cat(paste("\nload",file,"\n"))
		load(file)
		#str(df.all)
		
		#get correct order of sequence SampleCodes corresponding to first seq of Patient
		seq.PROT.nam<- as.character( df.all[,SeqPROT] )
		seq.PROT.nam[ which(is.na(seq.PROT.nam)) ]	<- "NA"
		seq.RT.nam	<- as.character( df.all[,SeqRT] )
		seq.RT.nam[ which(is.na(seq.RT.nam)) ]		<- "NA"
					
		pattern 	<- gsub('/',':',paste(signat.in,".clustalo$",sep=''))
		files		<- list.files(path=indir, pattern=pattern, full.names=1)
		#read all sequences and add a missing one with name "NA" 		
		seq.PROT	<- read.dna( files[ grep("PROT",files) ], format="fa", as.matrix=1 )
		tmp			<- as.DNAbin( matrix(rep('-',ncol(seq.PROT)),1,ncol(seq.PROT), dimnames=list(c("NA"),c())) )
		seq.PROT	<- rbind(seq.PROT,tmp)
		seq.RT		<- read.dna( files[ grep("RT",files) ], format="fa", as.matrix=1 )						 				
		tmp			<- as.DNAbin( matrix(rep('-',ncol(seq.RT)),1,ncol(seq.RT), dimnames=list(c("NA"),c())) )
		seq.RT		<- rbind(seq.RT,tmp)
		
		seq.PROT	<- seq.PROT[seq.PROT.nam,]
		seq.RT		<- seq.RT[seq.RT.nam,]
		rownames(seq.PROT)<- as.character( df.all[,Patient] )
		rownames(seq.RT)<- as.character( df.all[,Patient] )
		
		seq.PROT.RT	<- cbind(seq.PROT,seq.RT)
		if(verbose) print(seq.PROT.RT)
		file		<- paste(outdir,"/ATHENA_2013_03_FirstAliSequences_PROTRT_",gsub('/',':',signat.out),".R",sep='')
		if(verbose) cat(paste("\nwrite to",file))
		save(seq.PROT.RT, file=file)	
	}
	else
	{
		if(verbose)	cat(paste("\nload",file))
		load(file)
	}		
	if(0){	print("DEBUG FirstAliSequences"); seq.PROT.RT<- seq.PROT.RT[1:10,]	}
	if(0)	#create full fasta file with reference sequence and align once more
	{
		file		<- paste(outdir,"/ATHENA_2013_03_FirstAliSequences_HXB2PROTRT_",gsub('/',':',signat.out),".fasta",sep='')
		if(verbose) cat(paste("\nwrite to ",file))		
		data( refseq_hiv1_hxb2 )
		hxb2		<- data.table( hxb2 )
		hxb2		<- as.character( hxb2[, HXB2.K03455 ] )
		cat( paste(">HXB2\n",paste( hxb2[ seq.int(1,length(hxb2)-2) ], collapse='',sep='' ),"\n",sep=''), file=file, append=1)
		write.dna(seq.PROT.RT, file=file, format="fasta", append=1, colsep='', colw=length(hxb2)-2, blocksep=0)
				
		file<- paste("ATHENA_2013_03_FirstAliSequences_HXB2PROTRT_",gsub('/',':',signat.out),".fasta",sep='')
		cmd<- hivc.cmd.clustalo(outdir, file, signat='')
		system(cmd)		
	}
	if(0)	#curate alignment with reference 
	{
		file								<- paste(outdir,"/ATHENA_2013_03_FirstAliSequences_HXB2PROTRT_",gsub('/',':',signat.out),".fasta.clustalo",sep='')
		if(verbose) cat(paste("\nread ",file))
		seq.PROT.RT							<- read.dna(file, format="fasta", as.character=1)
		query.yes							<- hivc.seq.find(seq.PROT.RT, 2550, c("-","-","-"))
		if(length(query.yes))
			seq.PROT.RT[query.yes,2550:2552]<- matrix( c("c","c","y"), nrow=length(query.yes), ncol=3, byrow=1 )
		query.yes							<- hivc.seq.find(seq.PROT.RT, 2550, c("c","c","y","-"))
		if(length(query.yes))
			seq.PROT.RT[query.yes,2550:2553]<- matrix( c("-","-","-","-"), nrow=length(query.yes), ncol=4, byrow=1 )
		query.yes							<- hivc.seq.find(seq.PROT.RT, 2751, c("-","-","-","a"))
		if(length(query.yes))
			seq.PROT.RT[query.yes,2751:2754]<- matrix( c("a","-","-","-"), nrow=length(query.yes), ncol=4, byrow=1 )
		query.yes							<- hivc.seq.find(seq.PROT.RT, 3143, c("-"))
		for(i in query.yes)	
			seq.PROT.RT[i,3143:3151]		<-	c(seq.PROT.RT[i,3144:3151],"-") 	#align pos 3143 and move gap into 3rd codon
		query.yes							<- hivc.seq.find(seq.PROT.RT, 3151, c("a","-"))
		if(length(query.yes))
			seq.PROT.RT[query.yes,3151:3152]<- matrix( c("-","a"), nrow=length(query.yes), ncol=2, byrow=1 )
		#always delete drug resistance insertion at pos 69 because this is really cryptic
		seq.PROT.RT[,2752:2754]				<- matrix( c("-","-","-"), nrow=length(query.yes), ncol=3, byrow=1 )
		query.yes							<- hivc.seq.find(seq.PROT.RT, 3237, c("-"))
		for(i in query.yes)	
			seq.PROT.RT[i,3233:3237]<- c("-",seq.PROT.RT[i,3233:3236]) 		#align pos 3237 and move gap into 3rd codon
		query.yes							<- hivc.seq.find(seq.PROT.RT, 3433, c("-","t"))
		if(length(query.yes))
			seq.PROT.RT[query.yes,3433:3434]<- matrix( c("t","-"), nrow=length(query.yes), ncol=2, byrow=1 )
		query.yes							<- hivc.seq.find(seq.PROT.RT, 3433, c("-","-","-","-","t"))
		if(length(query.yes))
			seq.PROT.RT[query.yes,3433:3437]<- matrix( c("t","-","-","-","-"), nrow=length(query.yes), ncol=5, byrow=1 )
		query.yes							<- hivc.seq.find(seq.PROT.RT, 3434, c("-","-","-","-","a"))
		if(length(query.yes))
			seq.PROT.RT[query.yes,3434:3438]<- matrix( c("a","-","-","-","-"), nrow=length(query.yes), ncol=5, byrow=1 )
		query.yes							<- hivc.seq.find(seq.PROT.RT, 3503, c("-"))
		for(i in query.yes)
			seq.PROT.RT[i,3501:3503]		<- c("-",seq.PROT.RT[i,3501:3502])
		query.yes							<- hivc.seq.find(seq.PROT.RT, 3504, c("-","-","-","t"))	
		for(i in query.yes)
			seq.PROT.RT[i,3504:3524]		<- c(seq.PROT.RT[i,3507:3524],"-","-","-")
		query.yes							<- hivc.seq.find(seq.PROT.RT, 3506, c("-","-","-","t","g"))
		for(i in query.yes)
			seq.PROT.RT[i,3506:3526]		<- c(seq.PROT.RT[i,3509:3526],"-","-","-")
		query.yes							<- hivc.seq.find(seq.PROT.RT, 3507, c("-","-","-","g","a"))
		for(i in query.yes)
			seq.PROT.RT[i,3507:3537]		<-	c(seq.PROT.RT[i,3510:3537],"-","-","-")
		query.yes							<- hivc.seq.find(seq.PROT.RT, 3504, c("g","a","c"))
		for(i in query.yes)
			seq.PROT.RT[i,3504:3530]		<-	c("-","-","-",seq.PROT.RT[i,3504:3527])
		query.yes							<- hivc.seq.find(seq.PROT.RT, 3495, c("g","-","-","-","g"))
		if(length(query.yes))
			seq.PROT.RT[query.yes,3495:3499]<- matrix( c("g","g","-","-","-"), nrow=length(query.yes), ncol=5, byrow=1 )
		query.yes							<- hivc.seq.find(seq.PROT.RT, 3501, c("-","t","a"))
		if(length(query.yes))
			seq.PROT.RT[query.yes,3501:3503]<- matrix( c("t","a","-"), nrow=length(query.yes), ncol=3, byrow=1 )
		query.yes							<- hivc.seq.find(seq.PROT.RT, 3501, c("a","-"))
		if(length(query.yes))
			seq.PROT.RT[query.yes,3501:3502]<- matrix( c("-","a"), nrow=length(query.yes), ncol=2, byrow=1 )
		query.yes							<- hivc.seq.find(seq.PROT.RT, 3551, c("g","g"))
		for(i in query.yes)
			seq.PROT.RT[i,3551:3600]		<-	c("-",seq.PROT.RT[i,3551:3599])
		query.yes							<- hivc.seq.find(seq.PROT.RT, 3552, c("g","c"))
		for(i in query.yes)
			seq.PROT.RT[i,3552:3600]		<-	c("-",seq.PROT.RT[i,3552:3599])
		query.yes							<- hivc.seq.find(seq.PROT.RT, 3501, c("-","t","t"))
		for(i in query.yes)
			seq.PROT.RT[i,3501:3600]		<-	c("-",seq.PROT.RT[i,3501:3599])
		query.yes							<- hivc.seq.find(seq.PROT.RT, 3498, c("c","g","t"))
		for(i in query.yes)
			seq.PROT.RT[i,3498:3520]		<-	c("-",seq.PROT.RT[i,3498:3519])
		query.yes							<- hivc.seq.find(seq.PROT.RT, 3549, c("-","c","a"))
		if(length(query.yes))
			seq.PROT.RT[query.yes,3549:3551]<- matrix( c("c","a","-"), nrow=length(query.yes), ncol=3, byrow=1 )
		query.yes							<- hivc.seq.find(seq.PROT.RT, 3543, c("-","c","a","g","g","g","g","c","a"))
		if(length(query.yes))
			seq.PROT.RT[query.yes,3543:3551]<- matrix( c("c","a","g","g","g","g","c","a","-"), nrow=length(query.yes), ncol=9, byrow=1 )
		query.yes							<- hivc.seq.find(seq.PROT.RT, 3543, c("-","c","a","g","g","a","g","c","a"))
		if(length(query.yes))
			seq.PROT.RT[query.yes,3543:3551]<- matrix( c("c","a","g","g","a","g","c","a","-"), nrow=length(query.yes), ncol=9, byrow=1 )
		query.yes							<- hivc.seq.find(seq.PROT.RT, 3543, c("-","c","a","k","g","g","a","c","a"))
		if(length(query.yes))
			seq.PROT.RT[query.yes,3543:3551]<- matrix( c("c","a","k","g","g","a","c","a","-"), nrow=length(query.yes), ncol=9, byrow=1 )
		query.yes							<- hivc.seq.find(seq.PROT.RT, 3543, c("-","c","a","g","g","g","g","g","a"))
		if(length(query.yes))
			seq.PROT.RT[query.yes,3543:3551]<- matrix( c("c","a","g","g","g","g","g","a","-"), nrow=length(query.yes), ncol=9, byrow=1 )
		query.yes							<- hivc.seq.find(seq.PROT.RT, 3549, c("a","a","-"))
		if(length(query.yes))
			seq.PROT.RT[query.yes,3549:3551]<- matrix( c("-","-","-"), nrow=length(query.yes), ncol=3, byrow=1 )
		query.yes							<- hivc.seq.find(seq.PROT.RT, 3551, c("g","a","c"))
		for(i in query.yes)
			seq.PROT.RT[i,3551:3560]		<-	c("-",seq.PROT.RT[i,3551:3559])
		query.yes							<- hivc.seq.find(seq.PROT.RT, 3551, c("g","a","y","-"))
		if(length(query.yes))
			seq.PROT.RT[query.yes,3551:3554]<- matrix( c("-","g","a","y"), nrow=length(query.yes), ncol=4, byrow=1 )
		query.yes							<- hivc.seq.find(seq.PROT.RT, 3549, c("-","c","c","g","g"))
		if(length(query.yes))
			seq.PROT.RT[query.yes,3549:3551]<- matrix( c("c","-","c"), nrow=length(query.yes), ncol=3, byrow=1 )
		#some more manual editing at the end
		#cut before PROT pol start and at max len in database
		seq.PROT.RT							<- seq.PROT.RT[,2253:ncol(seq.PROT.RT)] 
		seq.PROT.RT							<- seq.PROT.RT[,1:1624]
		seq.PROT.RT.sort.by					<- apply(seq.PROT.RT,1,function(x)		which(rev(x)!="-" )[1]  )
		seq.PROT.RT							<- seq.PROT.RT[sort(seq.PROT.RT.sort.by, index.return=1)$ix,]
		#some more manual editing at all ends after sorting
		#always delete drug resistance at G333D/E because this is at the end of the alignment and alignment unreliable - could have picked up other ends
		seq.PROT.RT[,1294:1299]				<- matrix( c("-","-","-","-","-","-"), nrow=nrow(seq.PROT.RT), ncol=6, byrow=1 )

		seq.PROT.RT							<- as.DNAbin( seq.PROT.RT )
		file								<- paste(outdir,"/ATHENA_2013_03_FirstAliSequences_HXB2PROTRT_",gsub('/',':',signat.out),".fasta.clustalo2",sep='')		
		write.dna(seq.PROT.RT, file=file, format="fasta", colsep='', colw=ncol(seq.PROT.RT), blocksep=0)
	}
	if(0)	#create final HXB2PROTRT R and phylip files; need phylip for ExaML
	{
		signat.in	<- "Wed_May__1_17/08/15_2013"
		signat.out	<- "Sat_May_11_14:23:46_2013"
		file<- paste(outdir,"/ATHENA_2013_03_FirstAliSequences_HXB2PROTRT_",gsub('/',':',signat.in),".fasta.clustalo4",sep='')
		seq.PROT.RT	<- read.dna(file, format="fasta")
		print(seq.PROT.RT)
						
		seq.PROT.RT	<- seq.PROT.RT[ rownames(seq.PROT.RT)!="HXB2", ]
		
		file		<- paste(outdir,"/ATHENA_2013_03_FirstCurSequences_PROTRT_",gsub('/',':',signat.out),".R",sep='')
		if(verbose) cat(paste("\nwrite to",file))
		save(seq.PROT.RT, file=file)
		
		file		<- paste(outdir,"/ATHENA_2013_03_FirstCurSequences_PROTRT_",gsub('/',':',signat.out),".phylip",sep='')
		if(verbose) cat(paste("\nwrite to ",file))
		hivc.seq.write.dna.phylip(seq.PROT.RT, file=file)				
	}
	if(1)	#retain only third codon positions
	{
		signat.in	<- "Sat_May_11_14:23:46_2013"
		signat.out	<- "Sat_May_11_14:23:46_2013"

		file		<- paste(outdir,"/ATHENA_2013_03_FirstCurSequences_PROTRT_",gsub('/',':',signat.in),".R",sep='')
		if(verbose) cat(paste("\nload ",file))			
		load(file)
		
		codon3.idx	<- seq.int(3,ncol(seq.PROT.RT),3)
		seq.PROT.RT3<- seq.PROT.RT[, codon3.idx]
		file		<- paste(outdir,"/ATHENA_2013_03_FirstCurSequences_PROTRTCD3_",gsub('/',':',signat.out),".R",sep='')
		if(verbose) cat(paste("\nwrite to ",file))
		save(seq.PROT.RT3, file=file)
		
		file		<- paste(outdir,"/ATHENA_2013_03_FirstCurSequences_PROTRTCD3_",gsub('/',':',signat.out),".phylip",sep='')
		if(verbose) cat(paste("\nwrite to ",file))
		hivc.seq.write.dna.phylip(seq.PROT.RT3, file=file)
	}
}
######################################################################################
project.hivc.clustering.get.linked.and.unlinked<- function(dir.name= DATA)
{
	require(data.table)
	require(ape)
	if(1)	#extract unlinked and linked pairs -- this is version Sat_Jun_16_17/23/46_2013
	{
		verbose				<- 1
		
		#load all+enriched sequences
		indir				<- paste(dir.name,"tmp",sep='/')
		infile				<- "ATHENA_2013_03_CurAll+LANL_Sequences"
		insignat			<- "Sat_Jun_16_17/23/46_2013"		
		file				<- paste(indir,'/',infile,'_',gsub('/',':',insignat),".R",sep='')
		if(verbose)	cat(paste("\nread file",file))
		load(file)
		print(seq.PROT.RT)
		#exctract geographically distant seqs that are assumed to be truly unlinked to NL seqs
		seq.PROT.RT.TN		<- seq.PROT.RT[substr(rownames(seq.PROT.RT),1,2)=="TN", ]
		#extract ATHENA seqs
		seq.PROT.RT.NL		<- seq.PROT.RT[substr(rownames(seq.PROT.RT),1,2)!="TN", ]
		seq.PROT.RT.NL		<- seq.PROT.RT.NL[substr(rownames(seq.PROT.RT.NL),1,8)!="PROT+P51", ]
		
		#need patient id and PosSeq for each sample code		
		indir				<- paste(dir.name,"derived",sep='/')
		infile.seqinfo		<- "ATHENA_2013_03_Sequences"
		file				<- paste(indir,'/',paste(infile.seqinfo,".R", sep=''), sep='')
		if(verbose) cat(paste("\nloading file",file))
		load(file)	
		df.PROT				<- as.data.table(df[["PROT"]][,c("Patient","SampleCode","DateRes")])
		setkey(df.PROT,"SampleCode")
		df.RT				<- as.data.table(df[["RT"]][,c("Patient","SampleCode","DateRes")])
		setkey(df.RT,"SampleCode")
		df					<- merge(df.PROT,df.RT,all=TRUE)
		if( !all( df[,Patient.x==Patient.y], na.rm=1) ) stop("Patient names per sampling code in PROT and RT don t equal")
		if( !all( df[,DateRes.x==DateRes.y], na.rm=1) ) stop("Sampling dates per sampling code in PROT and RT don t equal")
		tmp					<- df[,DateRes.x]
		tmp2				<- df[,is.na(DateRes.x)]
		tmp[tmp2]			<- df[tmp2,DateRes.y]	
		df[,PosSeqT:= tmp]
		tmp					<- df[,Patient.x]
		tmp2				<- df[,is.na(Patient.x)]
		tmp[tmp2]			<- df[tmp2,Patient.y]	
		df[,Patient:= tmp]
		if(verbose) cat(paste("\nfound ATHENA seqs info, n=",nrow(df)))
		df.seqinfo			<- subset(df, !is.na(PosSeqT), select=c(SampleCode, Patient, PosSeqT))
		tmp					<- df.seqinfo[,gsub(' ','',SampleCode,fixed=1)]
		df.seqinfo[,FASTASampleCode:=tmp]
		
		if(verbose) cat(paste("\nfound ATHENA seqs with known sampling date, n=",nrow(df.seqinfo)))
		if(verbose) cat(paste("\nfound ATHENA patients whose sequences have known sampling date, n=",length(unique(df.seqinfo[,Patient]))))
		#extract list of truly linked sample codes
		setkey(df.seqinfo, "Patient")
		tmp					<- subset(df.seqinfo[, length(SampleCode)>1, by=Patient], V1==T, select=Patient)
		linked.bypatient	<- merge(tmp, df.seqinfo, all.x=1)
		setkey(linked.bypatient, "FASTASampleCode")
		linked.bypatient	<- subset(linked.bypatient, select=c(Patient, FASTASampleCode, PosSeqT))
		
		#extract list of truly unlinked seqs -- temporal separation
		indir				<- paste(dir.name,"derived",sep='/')
		infile.patient		<- "ATHENA_2013_03_All1stPatientCovariates"
		file				<- paste(indir,'/',paste(infile.patient,".R", sep=''), sep='')
		if(verbose) cat(paste("\nloading file",file))
		load(file)	
		#extract seqs of dead invidividuals
		df.dead						<- subset(df.all, !is.na(Died), c(Patient,Died))
		df.dead						<- merge(df.dead, df.seqinfo, by="Patient")
		setkey(df.dead,Died)
		#extract seroconverters
		df.serocon.acc				<- subset(df.all, NegT_Acc=="Yes" & NegT>=df.dead[1,Died],)
		#add seroconverters with inaccurate info
		df.serocon.nacc				<- subset(df.all, NegT_Acc=="No" & !is.na(NegT) & NegT>=df.dead[1,Died], )
		df.serocon.nacc.dy			<- subset(df.serocon.nacc, as.POSIXlt(NegT)$mday==15, )						#for inaccurate days, we (conservatively) assume the patient was only seronegative at the start of the month
		tmp							<- as.POSIXlt(df.serocon.nacc.dy[,NegT] )
		tmp$mday					<- 1
		df.serocon.nacc.dy[,NegT:=as.Date(tmp)]
		df.serocon.nacc.mody		<- subset(df.serocon.nacc, as.POSIXlt(NegT)$mon==6 & as.POSIXlt(NegT)$mday==1, )		#for inaccurate months and days, we (conservatively) assume the patient was only seronegative at the start of the year
		tmp							<- as.POSIXlt(df.serocon.nacc.mody[,NegT] )
		tmp$mon						<- 0
		df.serocon.nacc.mody[,NegT:=as.Date(tmp)]
		df.serocon					<- rbind(df.serocon.acc, df.serocon.nacc.dy, df.serocon.nacc.mody)		
		if(verbose) cat(paste("\nnumber of seroconverters with at least 1 preceeding dead HIV+",nrow(df.serocon)))
		#for each accurate seroconverter, extract HIV+ seqs that are dead before seroconversion
		if( any(as.logical(df.serocon[,is.na(NegT)])) )	warning("Found accurate seroconverters with missing NegT")
		df.serocon					<- merge(df.serocon,df.seqinfo,all.x=1,by="Patient")		
		setkey(df.serocon,NegT)
		unlinked.bytime				<- lapply(seq_len(nrow(df.serocon)), function(i)
				{
					tmp				<- subset(df.dead, Died<=df.serocon[i,NegT],select=c(Patient,FASTASampleCode,Died))
					tmp2			<- rep(df.serocon[i,FASTASampleCode],nrow(tmp))
					tmp[,query.FASTASampleCode:= tmp2]
					setkey(tmp,"FASTASampleCode")
					tmp					
				})
		names(unlinked.bytime)	<- df.serocon[,FASTASampleCode]			
		
		#need also NegT for each seq			
		df.seqinfo				<- merge(df.seqinfo, df.serocon[,list(NegT=min(NegT)),by=Patient], all.x=1, by="Patient")
		
		#extract list of truly unlinked seqs -- geographical separation
		unlinked.byspace		<- data.table(FASTASampleCode=rownames(seq.PROT.RT.TN), key="FASTASampleCode")
						
		outdir					<- paste(dir.name,"tmp",sep='/')
		outfile					<- "ATHENA_2013_03_Unlinked_and_Linked"
		outsignat				<- "Sat_Jun_16_17/23/46_2013"
		file					<- paste(outdir,'/',outfile,'_',gsub('/',':',insignat),".R",sep='')
		if(verbose) cat(paste("\nwrite linked and unlinked to file",file))
		save(unlinked.byspace,unlinked.bytime,linked.bypatient,df.seqinfo,file=file)
		
		#
		#some quick statistics
		#
		if(0)
		{
			#number of unlinked HIV+ by date (this date is when someone else is still seroneg)
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
		quit("no")
	}
	if(0)	#extract unlinked pairs by temporal separation -- this is version Fri_May_24_12/59/06_2013
	{
		verbose				<- 1
		unlinked.closest.n	<- NA
		indir				<- paste(dir.name,"derived",sep='/')
		outdir				<- paste(dir.name,"derived",sep='/')
		infile				<- "ATHENA_2013_03_All1stPatientCovariates.R"
		outfile				<- "ATHENA_2013_03_Unlinked_SeroConv_Dead"
		outsignat			<- "Fri_May_24_12/59/06_2013"
		file				<- paste(indir,infile,sep='/')
		if(verbose)
		{
			cat(paste("\nunlinked.closest.n",unlinked.closest.n))			
		}
		
		if(verbose)	cat(paste("\nread file",file))
		load(file)
		if(verbose) str(df.all)
		
		df.dead						<- subset(df.all, !is.na(Died), c(Patient,Died))
		setkey(df.dead,Died)
		df.serocon.acc				<- subset(df.all, NegT_Acc=="Yes" & NegT>=df.dead[1,Died],)
		if(1)
		{
			df.serocon.nacc				<- subset(df.all, NegT_Acc=="No" & !is.na(NegT) & NegT>=df.dead[1,Died], )
			#for inaccurate days, we (conservatively) assume the patient was only seronegative at the start of the month
			df.serocon.nacc.dy			<- subset(df.serocon.nacc, as.POSIXlt(NegT)$mday==15, )
			tmp							<- as.POSIXlt(df.serocon.nacc.dy[,NegT] )
			tmp$mday					<- 1
			df.serocon.nacc.dy[,NegT:=as.Date(tmp)]
			#for inaccurate months and days, we (conservatively) assume the patient was only seronegative at the start of the year
			df.serocon.nacc.mody		<- subset(df.serocon.nacc, as.POSIXlt(NegT)$mon==6 & as.POSIXlt(NegT)$mday==1, )
			tmp							<- as.POSIXlt(df.serocon.nacc.mody[,NegT] )
			tmp$mon						<- 0
			df.serocon.nacc.mody[,NegT:=as.Date(tmp)]
			#merge all
			df.serocon					<- rbind(df.serocon.acc, df.serocon.nacc.dy, df.serocon.nacc.mody)
		}
		else
			df.serocon					<- df.serocon.acc
		
		if(verbose) cat(paste("\nnumber of seroconverters with at least 1 preceeding dead HIV+",nrow(df.serocon)))
		#for each accurate seroconverter, extract HIV+ that are dead before seroconversion
		if( any(as.logical(df.serocon[,is.na(NegT)])) )	warning("Found accurate seroconverters with missing NegT")		
		setkey(df.serocon,NegT)
		unlinked.bytime				<- lapply(seq_len(nrow(df.serocon)), function(i)
				{
					tmp<- subset(df.dead, Died<=df.serocon[i,NegT],)											
					#tmp<- tmp[,Patient]
					#if(length(tmp)<unlinked.closest.n)	
					#	tmp<- c(rep(NA,unlinked.closest.n-length(tmp)),tmp)
					#rev(tmp)[1:unlinked.closest.n]
					if(is.na(unlinked.closest.n))	return( tmp[,Patient] )
					else							return( rev(tmp)[1:min(length(tmp),unlinked.closest.n)] )
				})
		names(unlinked.bytime)	<- df.serocon[,Patient]			
		#
		#some quick statistics
		#
		if(1)
		{
			#number of unlinked HIV+ by date (this date is when someone else is still seroneg)
			y		<- sapply(unlinked.bytime, function(x) length(x) )
			x		<- df.serocon[,NegT]
			xlim	<- range(x)
			tmp		<- as.POSIXlt(xlim[1])
			tmp$mday<- 1
			tmp$mon	<- 1
			xlim[1]	<- as.Date(tmp)
			plot(1,1,type='n',bty='n',xlab="time of seroconversion",ylab="number unlinked",xaxt='n',xlim=xlim,ylim=range(y))
			axis.Date(1,Year,at=seq(xlim[1], xlim[2],by="12 months"))
			polygon(c(x,c(x[length(x)],x[1])), c(y,0,0), border=NA, col="grey60" )
			
			#average PosT for all unlinked patients by date (this date is when someone else is still seroneg)
			unlinked.PosT	<- sapply(seq_along(unlinked.bytime), function(i)
					{
						tmp<- data.table(Patient=unlinked.bytime[[i]], key="Patient")
						tmp<- subset(df.all[tmp], select=c(Patient, PosT))
						tmp<- tmp[,mean(as.numeric(difftime(df.serocon[i,NegT],PosT,units="weeks"))/52,na.rm=1)]
						tmp
					})
			y		<- unlinked.PosT				
			par(mar=c(5,6,1,1))
			xlim[1]	<- as.Date("2000-01-01")
			plot(1,1,type='n',bty='n',xlab="time of seroconversion",ylab="median time difference\nbetween unlinked sequences [yrs]",xaxt='n',xlim=xlim,ylim=range(y))
			axis.Date(1,Year,at=seq(xlim[1], xlim[2],by="12 months"))
			lines(x,y)
		}
		#
		#print(unlinked.bytime.n)
		#print(any(diff(unlinked.bytime.n)<0))
		#
		#save
		#
		if(is.na(unlinked.closest.n))			
			file						<- paste(outdir,paste(outfile,"_UnlinkedAll_",gsub('/',':',outsignat),".R",sep=''),sep='/')
		else
			file						<- paste(outdir,paste(outfile,"_Unlinked",unlinked.closest.n,'_',gsub('/',':',outsignat),".R",sep=''),sep='/')
		if(verbose) cat(paste("\nwrite unlinked pairs to file",file))
		save(unlinked.bytime, df.serocon, df.all, file=file)		
		quit("no")
	}
}
######################################################################################
project.hivc.clustering.compare.NoDR.to.NoRecombNoDR<- function()
{	
	verbose		<- 1
	resume		<- 1
	patient.n	<- 15700; 	thresh.brl<- 0.096; 	thresh.bs<- 0.8
	indircov	<- paste(DATA,"derived",sep='/')
	infilecov	<- "ATHENA_2013_03_AllSeqPatientCovariates"									
	#
	# get clusters for No Drug resistance mutations, single linkage criterion		
	#		
	indir			<- paste(DATA,"tmp",sep='/')		
	infile			<- "ATHENA_2013_03_NoDRAll+LANL_Sequences_examlbs500"			
	insignat		<- "Thu_Aug_01_17/05/23_2013"
	#
	argv			<<- hivc.cmd.preclustering(indir, infile, insignat, indircov, infilecov, resume=resume)				 
	argv			<<- unlist(strsplit(argv,' '))
	ndr.clu.pre		<- hivc.prog.get.clustering.precompute()
	#
	argv			<<- hivc.cmd.clustering.tptn(indir, infile, insignat, indircov, infilecov, opt.brl="dist.brl.casc", patient.n=patient.n, resume=resume)
	argv			<<- unlist(strsplit(argv,' '))
	ndr.clu.tptn	<- hivc.prog.get.clustering.TPTN(clu.pre=clu.pre)
	#			
	argv			<<- hivc.cmd.clustering(indir, infile, insignat, opt.brl="dist.brl.casc", thresh.brl, thresh.bs, resume=resume)				 
	argv			<<- unlist(strsplit(argv,' '))
	ndr.clu			<- hivc.prog.get.clustering()
	#
	# get clusters for No Recombination + No Drug resistance mutations, single linkage criterion		
	#						
	infile			<- "ATHENA_2013_03_NoRCDRAll+LANL_Sequences_examlbs500"			
	insignat		<- "Fri_Nov_01_16/07/23_2013"
	#
	argv			<<- hivc.cmd.preclustering(indir, infile, insignat, indircov, infilecov, resume=resume)				 
	argv			<<- unlist(strsplit(argv,' '))
	nrc.clu.pre		<- hivc.prog.get.clustering.precompute()
	#
	argv			<<- hivc.cmd.clustering.tptn(indir, infile, insignat, indircov, infilecov, opt.brl="dist.brl.casc", patient.n=patient.n, resume=resume)
	argv			<<- unlist(strsplit(argv,' '))
	nrc.clu.tptn	<- hivc.prog.get.clustering.TPTN(clu.pre=clu.pre)
	#			
	argv			<<- hivc.cmd.clustering(indir, infile, insignat, opt.brl="dist.brl.casc", thresh.brl, thresh.bs, resume=resume)				 
	argv			<<- unlist(strsplit(argv,' '))
	nrc.clu			<- hivc.prog.get.clustering()
	#
	#	compare distribution of bootstrap values
	#
	hist( ndr.clu.pre$ph.node.bs )
	hist( nrc.clu.pre$ph.node.bs, border="blue", add=1 )
	#
	#	compare TP
	#
	ndr.clu.tptn$tp.by.all
	nrc.clu.tptn$tp.by.all
	#
	#	compare FP
	#
	ndr.clu.tptn$fpn.by.sum
	nrc.clu.tptn$fpn.by.sum		
	#
}
######################################################################################
project.hivc.clustering.computeclusterstatistics.fordistbrl<- function(thresh, clusters, with.withinpatientclu=0, dist.brl="dist.brl.med")
{
	#select clusters for analysis
	clusters			<- clusters[[dist.brl]]		
	patient.nNL			<- 15700	#estimate from 2012 annual report	
	#length(df.seqinfo[,unique(Patient)])		
	patient.n			<- patient.nNL	 + length(which(substr(ph$tip.label,1,2)=="TN")) + length(which(substr(ph$tip.label,1,8)=="PROT+P51"))
	
	clusters.pairs		<- sapply(clusters,function(x){	table(x$size.tips)[1]  	})
	clusters.seqinclu	<- sapply(clusters,function(x){	sum(x$size.tips) / patient.n 	})
	clusters.NLseqinclu	<- sapply(clusters,function(x){	sum(x$size.tips) / patient.nNL 	})
	clusters.totalclu	<- sapply(clusters,function(x){	length(na.omit( unique(x$clu.mem) ))  	})	
	clusters.maxsize		<- sapply(clusters,function(x){	max(x$size.tips)  	})
	clusters.info		<- as.data.table( cbind(thresh,clusters.totalclu,clusters.pairs,clusters.maxsize,clusters.seqinclu,clusters.NLseqinclu) )
	clusters.info		<- subset(clusters.info, brl>0.025 & bs>0.6)	
	clusters.info.bs.sd	<- clusters.info[, list(sd.totalclu=sd(clusters.totalclu), sd.pairs=sd(clusters.pairs), sd.maxsize=sd(clusters.maxsize), sd.seqinclu=sd(clusters.seqinclu), sd.NLseqinclu=sd(clusters.seqinclu)), by=bs]
	clusters.info.brl.sd<- clusters.info[, list(sd.totalclu=sd(clusters.totalclu), sd.pairs=sd(clusters.pairs), sd.maxsize=sd(clusters.maxsize), sd.seqinclu=sd(clusters.seqinclu), sd.NLseqinclu=sd(clusters.seqinclu)), by=brl]
	print(clusters.info)
	print(clusters.info.bs.sd)
	print(clusters.info.brl.sd)	
}		
######################################################################################
hivc.prog.get.clustering.MSM<- function(clu.pre= NULL)
{
	verbose		<- 1
	resume		<- 1
	indir		<- paste(DATA,"tmp",sep='/')		
	infile		<- "ATHENA_2013_03_NoDRAll+LANL_Sequences_examlbs100"			
	insignat	<- "Thu_Aug_01_17/05/23_2013"
	indircov	<- paste(DATA,"derived",sep='/')
	infilecov	<- "ATHENA_2013_03_AllSeqPatientCovariates"
	opt.brl		<- "dist.brl.casc" 
	thresh.brl	<- 0.096
	thresh.bs	<- 0.8
	
	if(exists("argv"))
	{
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,6),
									indir= return(substr(arg,8,nchar(arg))),NA)	}))
		if(length(tmp)>0) indir<- tmp[1]		
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,7),
									infile= return(substr(arg,9,nchar(arg))),NA)	}))
		if(length(tmp)>0) infile<- tmp[1]				
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,9),
									insignat= return(substr(arg,11,nchar(arg))),NA)	}))
		if(length(tmp)>0) insignat<- tmp[1]		
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,9),
									indircov= return(substr(arg,11,nchar(arg))),NA)	}))
		if(length(tmp)>0) indircov<- tmp[1]		
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,10),
									infilecov= return(substr(arg,12,nchar(arg))),NA)	}))
		if(length(tmp)>0) infilecov<- tmp[1]		
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,7),
									resume= return(as.numeric(substr(arg,9,nchar(arg)))),NA)	}))
		if(length(tmp)>0) resume<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,2),
									v= return(as.numeric(substr(arg,4,nchar(arg)))),NA)	}))
		if(length(tmp)>0) verbose<- tmp[1]
		
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,8),
									opt.brl= return(substr(arg,10,nchar(arg))),NA)	}))
		if(length(tmp)>0) opt.brl<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,10),
									thresh.bs= return(as.numeric(substr(arg,12,nchar(arg)))),NA)	}))
		if(length(tmp)>0) thresh.bs<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,11),
									thresh.brl= return(as.numeric(substr(arg,13,nchar(arg)))),NA)	}))
		if(length(tmp)>0) thresh.brl<- tmp[1]		
	}	
	if(verbose)
	{
		print(indir)
		print(infile)
		print(insignat)
		print(indircov)
		print(infilecov)
		print(resume)
		print(opt.brl)
		print(thresh.brl)
		print(thresh.bs)		
	}	
	#stop()
	if(resume)												#//load if there is R Master data.table
	{
		outdir			<- indir
		outfile			<- paste(infile,"_clust_",opt.brl,"_bs",thresh.bs*100,"_brl",thresh.brl*100,'_',"msmexpgr",sep='')
		outsignat		<- insignat	
		file			<- paste(outdir,'/',outfile,'_',gsub('/',':',outsignat),".R",sep='')		
		options(show.error.messages = FALSE)		
		readAttempt<-try(suppressWarnings(load(file)))
		if(!inherits(readAttempt, "try-error"))	cat(paste("\nresumed file",file))			
		options(show.error.messages = TRUE)		
	}
	if(!resume || inherits(readAttempt, "try-error"))
	{			
		#
		# precompute clustering stuff		
		#
		if(is.null(clu.pre))
		{
			argv		<<- hivc.cmd.preclustering(indir, infile, insignat, indircov, infilecov, resume=resume)				 
			argv		<<- unlist(strsplit(argv,' '))
			clu.pre		<- hivc.prog.get.clustering.precompute()
		}
		if(0)
		{
			clu.pre$df.seqinfo<- merge(df.all, subset(clu.pre$df.seqinfo,select=c(FASTASampleCode, Node)), all.y=1, by="FASTASampleCode")
			ph<- clu.pre$ph; dist.brl.max<- clu.pre$dist.brl.max; dist.brl.med<- clu.pre$dist.brl.med; dist.brl.casc<- clu.pre$dist.brl.casc; ph.node.bs<- clu.pre$ph.node.bs; ph.linked<- clu.pre$ph.linked; ph.unlinked.info<- clu.pre$ph.unlinked.info; ph.unlinked<- clu.pre$ph.unlinked; df.seqinfo<- clu.pre$df.seqinfo; unlinked.byspace<- clu.pre$unlinked.byspace; unlinked.bytime<- clu.pre$unlinked.bytime; linked.bypatient<- clu.pre$linked.bypatient	
			save(ph, dist.brl.max, dist.brl.med, dist.brl.casc, ph.node.bs, ph.linked, ph.unlinked.info, ph.unlinked, df.seqinfo, unlinked.byspace, unlinked.bytime, linked.bypatient, file="/Users/Oliver/duke/2013_HIV_NL/ATHENA_2013/data/tmp/ATHENA_2013_03_NoDRAll+LANL_Sequences_examlbs100_preclust_Thu_Aug_01_17:05:23_2013.R")
		}
		#
		#precompute clustering stuff for particular thresholds etc	
		argv		<<- hivc.cmd.clustering(indir, infile, insignat, opt.brl, thresh.brl, thresh.bs, resume=resume)				 
		argv		<<- unlist(strsplit(argv,' '))
		clu			<- hivc.prog.get.clustering()
		if(0)
		{
			clu$df.seqinfo<- merge(clu.pre$df.seqinfo,subset(clu$df.seqinfo,select=c(FASTASampleCode, cluster)), all.y=1, by="FASTASampleCode")
			set(clu$df.seqinfo, which( clu$df.seqinfo[,substr(FASTASampleCode,1,2)=="TN"] ), "CountryInfection", "FRGNTN")			
			set(clu$df.seqinfo, which( clu$df.seqinfo[,substr(FASTASampleCode,1,8)=="PROT+P51"] ), "CountryInfection", "FRGN")
			df.seqinfo<- clu$df.seqinfo; clustering<- clu$clustering; clusters.tp<- clu$clusters.tp; clusters.tn<- clu$clusters.tn
			save(df.seqinfo, clustering, clusters.tp, clusters.tn,file="/Users/Oliver/duke/2013_HIV_NL/ATHENA_2013/data/tmp/ATHENA_2013_03_NoDRAll+LANL_Sequences_examlbs100_clust_dist.brl.casc_bs80_brl9.6_Thu_Aug_01_17:05:23_2013.R")		
		}
		#
		# remove singletons
		#
		if(verbose) cat(paste("\nnumber of seq in tree is n=", nrow(clu$df.cluinfo)))
		df.cluinfo	<- subset(clu$df.seqinfo, !is.na(cluster) )
		if(verbose) cat(paste("\nnumber of seq in clusters is n=", nrow(df.cluinfo)))
		if(verbose) cat(paste("\nnumber of clusters is n=", length(unique(df.cluinfo[,cluster]))))
		#
		# remove within patient clusters
		#
		tmp			<- subset(df.cluinfo[,list(clu.is.bwpat=length(unique(Patient))>1),by="cluster"], clu.is.bwpat, cluster )
		df.cluinfo	<- merge(tmp, df.cluinfo, by="cluster", all.x=1)
		if(verbose) cat(paste("\nnumber of seq in clusters between patients is n=", nrow(df.cluinfo)))
		if(verbose) cat(paste("\nnumber of clusters between patients is n=", length(unique(df.cluinfo[,cluster]))))				
		if(0)	#plot clusters that have multiple foreign infections
		{
			outdir			<- indir
			outfile			<- paste(infile,"_clust_",opt.brl,"_bs",thresh.bs*100,"_brl",thresh.brl*100,'_',"multifrgninfection",sep='')
			outsignat		<- insignat									
			hivc.clu.getplot.multifrgninfection(clu.pre$ph, clu$clustering, df.cluinfo, plot.file=paste(outdir,'/',outfile,'_',gsub('/',':',outsignat),".pdf",sep=''))
		}					
		if(0)	#plot clusters with mixed exposure group
		{			
			outdir		<- indir
			outfile		<- paste(infile,"_clust_",opt.brl,"_bs",thresh.bs*100,"_brl",thresh.brl*100,'_',"mixedexpgr",sep='')
			outsignat	<- insignat		
			hivc.clu.getplot.mixedexposuregroup( clu.pre$ph, clu$clustering, df.cluinfo, plot.file=paste(outdir,'/',outfile,'_',gsub('/',':',outsignat),".pdf",sep='') )
		}	
		#
		# remove clusters with all seq coming from frgn infection
		#
		tmp			<- hivc.clu.getplot.excludeallmultifrgninfection(clu.pre$ph, clu$clustering, df.cluinfo )
		ph			<- tmp$cluphy
		df.cluinfo	<- tmp$cluphy.df
		clustering	<- tmp$cluphy.clustering		
		#
		#get in-country clusters. this splits clusters with a foreign sequence
		#		
		outdir			<- indir
		outfile			<- paste(infile,"_clust_",opt.brl,"_bs",thresh.bs*100,"_brl",thresh.brl*100,'_',"incountry",sep='')
		outsignat		<- insignat									
		#ph			<- clu.pre$ph; clustering	<- clu$clustering; plot.file=paste(outdir,'/',outfile,'_',gsub('/',':',outsignat),".pdf",sep=''); char.frgn  	='CountryInfection=="FRGN"'; char.frgntn	='CountryInfection=="FRGNTN"'; 
		incountry		<- hivc.clu.getplot.incountry(ph, clustering, df.cluinfo, plot.file=paste(outdir,'/',outfile,'_',gsub('/',':',outsignat),".pdf",sep=''))
		ph				<- incountry$cluphy
		clustering		<- incountry$clustering
		df.cluinfo		<- incountry$df.cluinfo		
		#
		# get msm exposure group clusters. this splits clusters with HET-F
		#
		set(df.cluinfo, which( df.cluinfo[,Trm%in%c("BLOOD","BREAST","PREG","NEEACC")] ), "Trm", "OTH" )
		set(df.cluinfo, which( df.cluinfo[,Trm=="HETfa"] ), "Trm", "HET" )		
		set(df.cluinfo, NULL, "Trm", factor(df.cluinfo[,Trm]) )		
		#ph<- incountry$cluphy; 		plot.file	<- paste(outdir,'/',outfile,'_',gsub('/',':',outsignat),".pdf",sep=''); levels.msm=c("BI","MSM","IDU","NA"); levels.het=c("BI","HET","IDU","NA"); levels.mixed=c("BI","MSM","HET","IDU","NA"); levels.oth="OTH"
		outdir			<- indir
		outfile			<- paste(infile,"_clust_",opt.brl,"_bs",thresh.bs*100,"_brl",thresh.brl*100,'_',"msmexpgr",sep='')
		outsignat		<- insignat							
		msm				<- hivc.clu.getplot.msmexposuregroup(ph, clustering, df.cluinfo, plot.file=paste(outdir,'/',outfile,'_',gsub('/',':',outsignat),".pdf",sep=''))		
		#
		# collapse within patient subclades in each subtree
		#
		msm				<- hivc.clu.collapse.monophyletic.withinpatientseq(msm$cluphy.subtrees, msm$cluphy.df )
		ph				<- msm$cluphy
		df.cluinfo		<- msm$cluphy.df
		cluphy.subtrees	<- msm$cluphy.subtrees		
		clustering		<- msm$cluphy.clustering
		#
		# get statistics of clusters that describe acute/early infection, branch lengths (#mutations)
		#				
		msm.brl.bwpat				<- hivc.clu.brl.bwpat(cluphy.subtrees, df.cluinfo)
		tmp							<- sapply(msm.brl.bwpat, function(x) c(median(na.omit(x)), sd(na.omit(x))))
		if(any(is.na(tmp[1,])))	stop("unexpected NA in median branch length - all within patient clusters removed ?")
		msm.clusu.df				<- data.table(cluster= as.numeric( names(cluphy.subtrees) ), clu.bwpat.medbrl= tmp[1,])		
		# add number of patients in cluster to 'msm$cluphy.df'
		tmp							<- df.cluinfo[, list(	cluster			= unique(cluster),
															nseq			= length(FASTASampleCode),
															FrgnInfection	= (!is.na(CountryInfection) & CountryInfection!="NL")[1],
															PossAcute		= (isAcute%in%c("Yes","Maybe"))[1],
															AnyPos_T1		= AnyPos_T1[1]
													), by="Patient"]
		tmp							<- tmp[, list(		clu.npat			= length(Patient), 
														clu.ntip			= sum(nseq),
														clu.nFrgnInfection	= length(which(FrgnInfection)),
														clu.fPossAcute		= length(which(PossAcute)) / length(Patient),
														clu.AnyPos_T1		= min(AnyPos_T1)						
												), by="cluster"]
		msm.clusu.df				<- merge(tmp, msm.clusu.df, by="cluster")
		#
		# re-order cluphy by branch length and plot
		#
		df.cluinfo				<- merge( df.cluinfo, msm.clusu.df, by="cluster" )
		setkey(df.cluinfo,clu.bwpat.medbrl)											#sort clusters by median branch length
		cluphy.subtrees			<- lapply( as.character(unique(df.cluinfo[,cluster])), function(name) cluphy.subtrees[[name]] )
		names(cluphy.subtrees)	<- as.character(unique(df.cluinfo[,cluster]))
		#
		# plot
		#		
		outdir			<- indir
		outfile			<- paste(infile,"_clust_",opt.brl,"_bs",thresh.bs*100,"_brl",thresh.brl*100,'_',"msmexpgr_bybwpatmedbrl",sep='')
		outsignat		<- insignat									
		tmp				<- hivc.clu.polyphyletic.clusters(df.cluinfo, cluphy.subtrees=cluphy.subtrees, plot.file=paste(outdir,'/',outfile,'_',gsub('/',':',outsignat),".pdf",sep=''), pdf.scaley=35, adj.tiplabel= c(-0.05,0.5), cex.tiplabel=0.3, pdf.xlim=0.36)
		cluphy			<- tmp$cluphy
		#
		# save
		#
		outdir			<- indir
		outfile			<- paste(infile,"_clust_",opt.brl,"_bs",thresh.bs*100,"_brl",thresh.brl*100,'_',"msmexpgr",sep='')
		outsignat		<- insignat	
		file			<- paste(outdir,'/',outfile,'_',gsub('/',':',outsignat),".R",sep='')
		if(verbose)	cat(paste("\nsave msm output to",file))
		save(df.cluinfo, cluphy.subtrees, clustering, cluphy, file=file)
	}
	list(df.cluinfo=df.cluinfo, cluphy=cluphy, cluphy.subtrees=cluphy.subtrees, clustering=clustering)
}
######################################################################################
hivc.prog.get.clustering.TPTN<- function(clu.pre= NULL)
{
	require(RColorBrewer)
	require(colorspace)
	indir		<- paste(DATA,"tmp",sep='/')		
	infile		<- "ATHENA_2013_03_NoDRAll+LANL_Sequences_examlbs100"			
	insignat	<- "Thu_Aug_01_17/05/23_2013"
	indircov	<- paste(DATA,"derived",sep='/')
	infilecov	<- "ATHENA_2013_03_AllSeqPatientCovariates"							
	
	patient.n	<- 15700
	opt.brl		<- "dist.brl.casc"
	thresh.brl	<- c(seq(0.02,0.05,0.01),seq(0.06,0.12,0.02),seq(0.16,0.24,0.04))
	thresh.bs	<- c(0.7,0.75,0.8,0.85,0.9,0.95)
	resume		<- 1
	verbose		<- 1
	
	if(exists("argv"))
	{
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,6),
									indir= return(substr(arg,8,nchar(arg))),NA)	}))
		if(length(tmp)>0) indir<- tmp[1]		
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,7),
									infile= return(substr(arg,9,nchar(arg))),NA)	}))
		if(length(tmp)>0) infile<- tmp[1]				
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,9),
									insignat= return(substr(arg,11,nchar(arg))),NA)	}))
		if(length(tmp)>0) insignat<- tmp[1]		
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,9),
									indircov= return(substr(arg,11,nchar(arg))),NA)	}))
		if(length(tmp)>0) indircov<- tmp[1]		
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,10),
									infilecov= return(substr(arg,12,nchar(arg))),NA)	}))
		if(length(tmp)>0) infilecov<- tmp[1]		
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,7),
									resume= return(as.numeric(substr(arg,9,nchar(arg)))),NA)	}))
		if(length(tmp)>0) resume<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,2),
									v= return(as.numeric(substr(arg,4,nchar(arg)))),NA)	}))
		if(length(tmp)>0) verbose<- tmp[1]
		
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,8),
									opt.brl= return(substr(arg,10,nchar(arg))),NA)	}))
		if(length(tmp)>0) opt.brl<- tmp[1]
		
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,10),
									patient.n= return(as.numeric(substr(arg,12,nchar(arg)))),NA)	}))
		if(length(tmp)>0) patient.n<- tmp[1]
	}	
	outdir			<- indir
	outsignat		<- insignat
	outfile			<- paste(infile,"_clusttptn_",opt.brl,sep='')	
	if(verbose)
	{
		print(indir)
		print(infile)
		print(insignat)
		print(indircov)
		print(infilecov)
		print(resume)
		print(opt.brl)
		print(patient.n)
	}
	#precompute clustering stuff	
	if(is.null(clu.pre))
	{
		argv	<<- hivc.cmd.preclustering(indir, infile, insignat, indircov, infilecov, resume=resume)				 
		argv	<<- unlist(strsplit(argv,' '))
		clu.pre	<- hivc.prog.get.clustering.precompute()
	}
	dist.brl	<- switch(	opt.brl, 
							"dist.brl.max"		= clu.pre$dist.brl.max,
							"dist.brl.med"		= clu.pre$dist.brl.med,
							"dist.brl.casc"		= clu.pre$dist.brl.casc,
							NA)
	if(any(is.na(dist.brl)))	stop("unexpected NA in dist.brl")
	#
	# evaluate TP and TN for several clustering thresholds	 
	#
	file		<- paste(outdir, '/', outfile, '_', gsub('/',':',outsignat),".R", sep='' )
	if(resume)												#//load if there is R Master data.table
	{
		options(show.error.messages = FALSE)		
		readAttempt<-try(suppressWarnings(load(file)))
		if(!inherits(readAttempt, "try-error"))	cat(paste("\nresumed file",file))			
		options(show.error.messages = TRUE)		
	}
	if(!resume || inherits(readAttempt, "try-error"))		#else generate R Master data.table
	{
		cat(paste("\ncreate file",file))
		clusters		<- lapply(thresh.bs, function(bs)
							{
								clusters	<- lapply(thresh.brl, function(brl)
										{
											hivc.clu.clusterbythresh(clu.pre$ph,thresh.brl=brl,dist.brl=dist.brl,thresh.nodesupport=bs,nodesupport=clu.pre$ph.node.bs,retval="all")
										})
								clusters.tp	<- lapply(clusters,function(x)
										{
											hivc.clu.truepos(x, clu.pre$ph.linked, Ntip(clu.pre$ph), verbose=0)
										})		
								clusters.fp	<- lapply(clusters,function(x)
										{
											hivc.clu.trueneg(x, clu.pre$ph.unlinked.info, clu.pre$ph.unlinked, Ntip(clu.pre$ph), verbose=0)
										})
								list(clu= clusters, tp= clusters.tp, fp= clusters.fp)							
							})	
		names(clusters)	<- thresh.bs	
		#
		# get number of patients in clustering	 
		#
		df.cluinfo						<- copy(clu.pre$df.seqinfo)
		df.cluinfo[,"cluster":= NA]
		setkey(df.cluinfo, Node)
		clusters.nbwpat					<- t( sapply(seq_along(clusters),function(i)
											{				
												sapply(seq_along(clusters[[i]][["clu"]]), function(j)
														{
															clustering			<- clusters[[i]][["clu"]][[j]]
															set(df.cluinfo, NULL, "cluster", clustering[["clu.mem"]][seq_len(Ntip(clu.pre$ph))] )
															df.bwpatclu			<- subset(df.cluinfo, !is.na(cluster) )[ , list(n.pat= length(unique(na.omit(Patient))) ), by="cluster"]
															df.bwpatclu			<- subset(df.bwpatclu, n.pat>1)
															sum(df.bwpatclu[,n.pat]) 							
														})
											}) )
		rownames(clusters.nbwpat)		<- thresh.bs
		colnames(clusters.nbwpat)		<- thresh.brl	
		#
		# get number of sequences in clustering	 
		#	
		clusters.nseq					<- t( sapply(seq_along(clusters),function(i)
											{				
												sapply(seq_along(clusters[[i]][["clu"]]), function(j)
														{
															clustering			<- clusters[[i]][["clu"]][[j]]
															set(df.cluinfo, NULL, "cluster", clustering[["clu.mem"]][seq_len(Ntip(clu.pre$ph))] )
															nrow( subset( df.cluinfo, 	!is.na(cluster) & substr(FASTASampleCode, 1, 2)!="TN" & substr(FASTASampleCode, 1, 8)!="PROT+P51" ) ) 							
														})
											}) )
		rownames(clusters.nseq)		<- thresh.bs
		colnames(clusters.nseq)		<- thresh.brl	
		#
		# get distribution of patients in clustering	 
		#
		clusters.dbwpat				<- lapply(seq_along(clusters),function(i)
										{				
											tmp			<- lapply(seq_along(clusters[[i]][["clu"]]), function(j)
															{
																clustering			<- clusters[[i]][["clu"]][[j]]
																set(df.cluinfo, NULL, "cluster", clustering[["clu.mem"]][seq_len(Ntip(clu.pre$ph))] )
																df.bwpatclu			<- subset(df.cluinfo, !is.na(cluster) )[ , list(n.pat= length(unique(na.omit(Patient))), size=length(FASTASampleCode) ), by="cluster"]
																
																df.bwpatclu			<- df.cluinfo[ , list(n.pat= length(unique(na.omit(Patient))), size=length(FASTASampleCode) ), by="cluster"]													
																t.bwpatclu			<- table( subset(df.bwpatclu, !is.na(cluster) & n.pat>0)[, n.pat] )		#n.pat>0 excludes clusters of foreign seq
																#add non-clustering ATHENA seqs to '1'; these have cluster==NA and Patient!=NA 
																t.bwpatclu[1]		<- t.bwpatclu[1] + subset(df.bwpatclu, is.na(cluster))[,n.pat]
																t.bwpatclu
															})
											names(tmp)	<- thresh.brl
											tmp
										})
		names(clusters.dbwpat)		<- thresh.bs
	
		file			<- paste(outdir, '/', outfile, '_', gsub('/',':',outsignat),".R", sep='')
		save(clusters,clusters.nbwpat,clusters.dbwpat,clusters.nseq,file=file)
	}
	clusters.cov.epidemic	<- clusters.nbwpat / patient.n
	clusters.cov.patindb	<- clusters.nbwpat / length(unique(subset(df.cluinfo,!is.na(Patient))[,Patient]))	
	clusters.cov.seqindb	<- clusters.nseq / nrow( subset( df.cluinfo, 	substr(FASTASampleCode, 1, 2)!="TN" & substr(FASTASampleCode, 1, 8)!="PROT+P51" ) )
	#
	# plot cluster size distributions per bootstrap
	#
	with.singleton			<- 0
	thresh.brl.select		<- which(thresh.brl %in% c(0.02, 0.04, 0.06, 0.1, 0.16) )
	thresh.bs.select		<- which(thresh.bs %in% c(0.7, 0.8, 0.95) )
	pch						<- 15 + seq_along(thresh.bs.select)
	lapply(thresh.bs.select,function(i)
							{				
								file				<- paste(outdir, paste(infile,"_clustpdf_",opt.brl,"_bs",thresh.bs[i],'_',gsub('/',':',outsignat),".pdf",sep=''),sep='/')				
								pdf(file=file, width=5,height=5)				
								if( with.singleton )
								{
									x					<- lapply(thresh.brl.select, function(j){	as.numeric(names(clusters.dbwpat[[i]][[j]]))	})
									y					<- lapply(thresh.brl.select, function(j){	clusters.dbwpat[[i]][[j]]	})					
								}
								else
								{
									x					<- lapply(thresh.brl.select, function(j){	as.numeric(names(clusters.dbwpat[[i]][[j]]))[-1]	})
									y					<- lapply(thresh.brl.select, function(j){	clusters.dbwpat[[i]][[j]][-1]	})
								}
								xlim				<- range( c(70, unlist(x) ) )
								cols				<- brewer.pal( length(x), "Accent")
								plot(1,1,bty='n',type='n',xlim=xlim, ylim=range(unlist(y)), xlab="patients", ylab="frequency", log='y')
								dummy<- lapply(seq_along(x),function(j)
										{
											points(x[[j]],y[[j]],type='l',col=cols[j],pch=pch[i])
											points(x[[j]],y[[j]],col=cols[j],pch=pch[i], cex=0.5)
										})				
								legend("topright", 		fill=cols, 	legend= paste("BRL", thresh.brl[ thresh.brl.select ]), bty='n', border=NA)
								legend("bottomright", 	pch=pch[i], legend= paste("BS", thresh.bs[i] ), bty='n', border=NA)
								dev.off()				
							})
	#
	# plot cluster size distributions per brl
	#
	clusters.dbwpat			<- lapply(seq_along(clusters.dbwpat[[1]]), function(i)
									{
										tmp			<- lapply(seq_along(clusters.dbwpat), function(j) 		clusters.dbwpat[[j]][[i]]		)
										names(tmp)	<- thresh.bs
										tmp
									})
	names(clusters.dbwpat)	<- thresh.brl
	thresh.brl.select	<- which(thresh.brl %in% c(0.06, 0.08, 0.1) )
	thresh.bs.select	<- which(thresh.bs %in% c(0.7, 0.75, 0.8, 0.85, 0.9, 0.95) )
	pch					<- 15 + seq_along(thresh.brl.select)
	lapply(thresh.brl.select,function(i)
							{												
								file				<- paste(outdir, paste(infile,"_clustpdf_",opt.brl,"_brl",thresh.brl[i],'_',gsub('/',':',outsignat),".pdf",sep=''),sep='/')				
								pdf(file=file, width=5,height=5)				
								if( with.singleton )
								{
									x					<- lapply(thresh.bs.select, function(j){	as.numeric(names(clusters.dbwpat[[i]][[j]]))	})
									y					<- lapply(thresh.bs.select, function(j){	clusters.dbwpat[[i]][[j]]	})					
								}
								else
								{
									x					<- lapply(thresh.bs.select, function(j){	as.numeric(names(clusters.dbwpat[[i]][[j]]))[-1]	})
									y					<- lapply(thresh.bs.select, function(j){	clusters.dbwpat[[i]][[j]][-1]	})
								}
								xlim				<- range( c(70, unlist(x) ) )
								cols				<- brewer.pal( length(x), "Accent")
								plot(1,1,bty='n',type='n',xlim=xlim, ylim=range(unlist(y)), xlab="patients", ylab="frequency", log='y')
								dummy<- lapply(seq_along(x),function(j)
										{
											points(x[[j]],y[[j]],type='l',col=cols[j],pch=pch[i])
											points(x[[j]],y[[j]],col=cols[j], pch=pch[i], cex=0.5)
										})				
								legend("topright", 		pch=pch[i],	legend= paste("BRL", thresh.brl[ i ]), bty='n', border=NA)
								legend("bottomright", 	fill=cols,	legend= paste("BS", thresh.bs[ thresh.bs.select ] ), bty='n', border=NA)
								dev.off()				
							})					
	#
	# get TP and TN for several clustering thresholds
	#
	tp.by.all				<- t( sapply(seq_along(clusters),function(i){				sapply(clusters[[i]][["tp"]], function(x) x["tp.by.all"] )		}))
	rownames(tp.by.all)		<- thresh.bs
	colnames(tp.by.all)		<- thresh.brl	
	tpclu.by.all			<- t( sapply(seq_along(clusters),function(i){				sapply(clusters[[i]][["tp"]], function(x) x["tpclu.by.all"] )		}))
	rownames(tpclu.by.all)	<- thresh.bs
	colnames(tpclu.by.all)	<- thresh.brl
	fpn.by.sum				<- t( sapply(seq_along(clusters),function(i){		sapply(clusters[[i]][["fp"]], function(x) x["fpn.by.sum"] )	}) )
	rownames(fpn.by.sum)	<- thresh.bs
	colnames(fpn.by.sum)	<- thresh.brl
	fp.by.all				<- t( sapply(seq_along(clusters),function(i){		sapply(clusters[[i]][["fp"]], function(x) x["fp.by.all"] )	}) )
	rownames(fp.by.all)		<- thresh.bs
	#
	# plot TP and TN for several clustering thresholds
	#	
	colnames(fp.by.all)		<- thresh.brl	
	cols					<- diverge_hcl(nrow(fpn.by.sum), h = c(246, 40), c = 96, l = c(65, 90))
	names(cols)				<- thresh.bs
	#	 
	hivc.clu.plot.tptn(	fpn.by.sum, tp.by.all, paste(outdir, paste(outfile,"_tpbyall_fpnbysum_",gsub('/',':',outsignat),".pdf",sep=''),sep='/'), 
						cols, xlab= "#FP (among all)", ylab= "%TP (among all)", xlim= range(c(0,16,fpn.by.sum)))
	#
	hivc.clu.plot.tptn(	fp.by.all, tp.by.all, paste(outdir, paste(outfile,"_tpbyall_fpbyall_",gsub('/',':',outsignat),".pdf",sep=''),sep='/'), 
						cols, xlab= "%FP (among all)", ylab= "%TP (among all)", xlim= range(c(0,0.01,fp.by.all)))
	#
	hivc.clu.plot.tptn(	fpn.by.sum, tpclu.by.all, paste(outdir, paste(outfile,"_tpclubyall_fpnbysum_",gsub('/',':',outsignat),".pdf",sep=''),sep='/'), 
						cols, xlab= "#FP (among all)", ylab= "%TP (among clu)", xlim= range(c(0,16,fpn.by.sum)))
	#
	hivc.clu.plot.tptn(	fp.by.all, tpclu.by.all, paste(outdir, paste(outfile,"_tpclubyall_fpbyall_",gsub('/',':',outsignat),".pdf",sep=''),sep='/'), 
						cols, xlab= "%FP (among all)", ylab= "%TP (among clu)", xlim= range(c(0,0.01,fp.by.all)))
	#
	hivc.clu.plot.tptn(	fpn.by.sum, clusters.cov.epidemic, paste(outdir, paste(outfile,"_covepi_fpnbysum_",gsub('/',':',outsignat),".pdf",sep=''),sep='/'), 
						cols, xlab= "#FP (among all)", ylab= "%coverage (of epi)", xlim= range(c(0,16,fpn.by.sum)))
	#
	hivc.clu.plot.tptn(	fp.by.all, clusters.cov.epidemic, paste(outdir, paste(outfile,"_covepi_fpbyall_",gsub('/',':',outsignat),".pdf",sep=''),sep='/'), 
						cols, xlab= "%FP (among all)", ylab= "%coverage (of epi)", xlim= range(c(0,0.01,fp.by.all)))
	#
	hivc.clu.plot.tptn(	fpn.by.sum, clusters.cov.patindb, paste(outdir, paste(outfile,"_covpat_fpnbysum_",gsub('/',':',outsignat),".pdf",sep=''),sep='/'), 
						cols, xlab= "#FP (among all)", ylab= "%coverage (of db)", xlim= range(c(0,16,fpn.by.sum)))
	#
	hivc.clu.plot.tptn(	fp.by.all, clusters.cov.patindb, paste(outdir, paste(outfile,"_covpat_fpbyall_",gsub('/',':',outsignat),".pdf",sep=''),sep='/'), 
						cols, xlab= "%FP (among all)", ylab= "%coverage (of db)", xlim= range(c(0,0.01,fp.by.all)))
	#
	hivc.clu.plot.tptn(	fpn.by.sum, clusters.cov.seqindb, paste(outdir, paste(outfile,"_covseq_fpnbysum_",gsub('/',':',outsignat),".pdf",sep=''),sep='/'), 
						cols, xlab= "#FP (among all)", ylab= "%seq coverage (of db)", xlim= range(c(0,16,fpn.by.sum)))
	#
	hivc.clu.plot.tptn(	fp.by.all, clusters.cov.seqindb, paste(outdir, paste(outfile,"_covseq_fpbyall_",gsub('/',':',outsignat),".pdf",sep=''),sep='/'), 
						cols, xlab= "%FP (among all)", ylab= "%seq coverage (of db)", xlim= range(c(0,0.01,fp.by.all)))

	
	list(	clusters				= clusters,
			clusters.cov.epidemic	= clusters.cov.epidemic, 
			clusters.cov.patindb	= clusters.cov.patindb, 
			clusters.cov.seqindb	= clusters.cov.seqindb,
			clusters.dbwpat			= clusters.dbwpat,
			fpn.by.sum				= fpn.by.sum,
			fp.by.all				= fp.by.all,
			tp.by.all				= tp.by.all, 
			tpclu.by.all			= tpclu.by.all		
			)
}
######################################################################################
project.hivc.beast<- function(dir.name= DATA)
{
	require(ape)
	require(data.table)
	require(RColorBrewer)
	require(XML)
	if(0)	#get BEAST nexus file for seroconverters
	{		
		indir			<- paste(DATA,"tmp",sep='/')		
		infile			<- "ATHENA_2013_03_NoDRAll+LANL_Sequences"		
		insignat		<- "Thu_Aug_01_17/05/23_2013"
		indircov		<- paste(DATA,"derived",sep='/')
		infilecov		<- "ATHENA_2013_03_AllSeqPatientCovariates"
		infiletree		<- paste(infile,"examlbs100",sep="_")
		infilexml		<- paste(infile,'_',"beast",'_',"seroneg",sep='')
		
		outdir			<- indir
		outsignat		<- "Tue_Aug_26_09/13/47_2013"
		
		opt.brl			<- "dist.brl.casc" 
		thresh.brl		<- 0.096
		thresh.bs		<- 0.8
		pool.ntip		<- 130
		infilexml.opt	<- "txs4clu"
		resume			<- 1
		verbose			<- 1
		
		argv		<<- hivc.cmd.beast.poolrunxml(indir, infile, insignat, indircov, infilecov, infiletree, infilexml, outsignat, pool.ntip, infilexml.opt=infilexml.opt, opt.brl=opt.brl, thresh.brl=thresh.brl, thresh.bs=thresh.bs, resume=resume, verbose=1)
		cat(argv)
		argv		<<- unlist(strsplit(argv,' '))
		hivc.prog.BEAST.generate.xml()		
	}
	if(0)	#get distribution of TMRCAs of tip nodes for a particular BEAST run
	{
		indir				<- paste(DATA,"beast/beast_131011",sep='/')		
		infile				<- "ATHENA_2013_03_NoDRAll+LANL_Sequences_beast_seroneg"
		insignat			<- "Tue_Aug_26_09/13/47_2013"
		infilexml.opt		<- "mph4clutx4tipLdTd"
		infilexml.template	<- "um232rhU2045"
		verbose				<- 1
		
		
		tmp			<- list.files(indir, pattern=paste(".log$",sep=''))
		files		<- tmp[ grepl(paste('_',infilexml.template,'_',sep=''), tmp) & grepl(paste('_',infilexml.opt,'_',sep=''), tmp) & grepl(gsub('/',':',insignat), tmp) ]
		if(verbose)	cat(paste("\nFound files matching input args, n=", length(files)))
		
		#	read length of tip stems
		df.tstem	<- lapply( seq_along(files), function(i)
					{
						file.log	<- files[i]
						if(verbose)	cat(paste("\nprocess file,", file.log))
						file.log	<- paste(indir,file.log,sep='/')
						file.xml	<- paste(substr(file.log,1,nchar(file.log)-3), "xml", sep='')
						if(verbose)	cat(paste("\nand file,", file.xml))
						hivc.beast.read.log2tstem(file.log, file.xml, beastlabel.idx.samplecode=6, burn.in= 5e6, breaks.n= 30, verbose=1)					
					})
		df.tstem<- rbindlist( df.tstem )		
	}
	if(1)	#plot BEAST clusters
	{					
		indir				<- paste(DATA,"beast/beast_130908",sep='/')		
		infile				<- "ATHENA_2013_03_NoDRAll+LANL_Sequences"		
		insignat			<- "Tue_Aug_26_09/13/47_2013"
		indircov			<- paste(DATA,"derived",sep='/')
		infilecov			<- "ATHENA_2013_03_AllSeqPatientCovariates"
		infiletree			<- paste(infile,"examlbs100",sep="_")
		infilexml			<- paste(infile,'_',"beast",'_',"seroneg",sep='')
		#infilexml.template	<- "um22rhU2050"
		infilexml.template	<- "um22rhG202018"
		#infilexml.template	<- "rhU65rho753"
		#infilexml.template	<- "rhU65rho903"
		#infilexml.template	<- "rhU65rho906"
		#infilexml.template	<- "rhU65rho909"	
		#infilexml.template	<- "um181rhU2045"
		#infilexml.template	<- "um182rhU2045"
		#infilexml.template	<- "um183rhU2045"
		#infilexml.template	<- "um182us45"
		#infilexml.template	<- "um182us60"
		infilexml.opt		<- "txs4clu"
		#infilexml.opt		<- "txs4clufx03"
		#infilexml.opt		<- "mph4clu"
		#infilexml.opt		<- "mph4clumph4tu"
		#infilexml.opt		<- "mph4clufx03"
	
		argv				<<- hivc.cmd.beast.evalrun(indir, infilexml, insignat, infilexml.opt, infilexml.template, pool.n, verbose=verbose)
		argv				<<- unlist(strsplit(argv,' '))
		hivc.prog.BEAST.evalpoolrun()
	}
}
######################################################################################
project.ukca.TPTN.bootstrapvalues<- function(dir.name= DATA)
{
	require(ape)
	require(data.table)
	require(RColorBrewer)
	
	indir				<- paste(DATA,"tmp",sep='/')
	infile				<- "UKCA_2013_07_TNTPHIVnGTR_examlbs100"
	insignat			<- "Mon_Sep_22_17/23/46_2013"		
	file				<- paste(indir,'/',infile,'_',gsub('/',':',insignat),".newick",sep='')
	plot.file			<- paste(indir,'/',infile,"_preclust_distbs_",gsub('/',':',insignat),".pdf", sep='')
	ph					<- read.tree(file)
	ph$node.label		<- as.numeric(ph$node.label)/100
	ph$nodel.label[1]	<- 0
	ph$node.label[3]	<- 0.01		#little phylo cleaning hack ;-)
	ph.mrca				<- mrca(ph)
	dist.brl.casc		<- hivc.clu.brdist.stats(ph, eval.dist.btw="leaf", stat.fun=hivc.clu.min.transmission.cascade)
	
	infile				<- "UKCA_2013_07_B_true_pos"
	file				<- paste(indir,'/',infile,'_',gsub('/',':',insignat),".csv",sep='')
	df.tp				<- as.data.table(read.csv(file, header=F, sep=' ', col.names=c("tip1","tip2")))		
	set(df.tp,NULL,"tip1", as.character(df.tp[,tip1]))
	set(df.tp,NULL,"tip2", as.character(df.tp[,tip2]))
	if(verbose)	cat(paste("\nNumber of TP pairs read, n=",nrow(df.tp)))
	df.tp				<- subset(df.tp, df.tp[,tip1]%in%rownames(ph.mrca)  &  df.tp[,tip2]%in%rownames(ph.mrca))
	if(verbose)	cat(paste("\nNumber of TP pairs with tips in ph, n=",nrow(df.tp)))
	
	infile				<- "UKCA_2013_07_B_true_neg"
	file				<- paste(indir,'/',infile,'_',gsub('/',':',insignat),".csv",sep='')
	df.tn				<- as.data.table(read.csv(file, header=F, sep=' ', col.names=c("tip1","tip2")))		
	set(df.tn,NULL,"tip1", as.character(df.tn[,tip1]))
	set(df.tn,NULL,"tip2", as.character(df.tn[,tip2]))
	if(verbose)	cat(paste("\nNumber of TN pairs read, n=",nrow(df.tn)))		
	df.tn				<- subset(df.tn, df.tn[,tip1]%in%rownames(ph.mrca)  &  df.tn[,tip2]%in%rownames(ph.mrca))
	if(verbose)	cat(paste("\nNumber of TN pairs with tips in ph, n=",nrow(df.tn)))
	
	#set MRCAs
	df.tn[,dummy:=seq_len(nrow(df.tn))]
	df.tp[,dummy:=seq_len(nrow(df.tp))]
	df.tn				<- df.tn[, list(tip1=tip1, tip2=tip2, mrca= ph.mrca[tip1,tip2]),by="dummy"]
	df.tp				<- df.tp[, list(tip1=tip1, tip2=tip2, mrca= ph.mrca[tip1,tip2]),by="dummy"]
	
	bs.linked.bypatient	<- df.tp
	bs.unlinkedpairs	<- df.tn		
	
	hivc.phy.get.TP.and.TN.bootstrapvalues(ph, bs.linked.bypatient, ph.mrca=ph.mrca ,df.seqinfo=NULL, bs.unlinkedpairs=bs.unlinkedpairs, bs.unlinked.byspace=NULL, dist.brl=NULL, thresh.brl=0.096, plot.file=plot.file, verbose= 1)	
}
######################################################################################
project.athena.TPTN.bootstrapvalues<- function(dir.name= DATA)
{	
	require(ape)
	require(data.table)
	require(RColorBrewer)
	
	verbose		<- 1
	resume		<- 1
	
	# load sequences
	indir								<- paste(dir.name,"tmp",sep='/')
	infile								<- "ATHENA_2013_03_CurAll+LANL_Sequences"			
	insignat							<- "Sat_Jun_16_17/23/46_2013"
	file								<- paste(indir,'/',infile,'_',gsub('/',':',insignat),".R",sep='')			
	load(file)		
	#
	# precompute clustering stuff		
	#
	patient.n	<- 15700
	indir		<- paste(DATA,"tmp",sep='/')		
	infile		<- "ATHENA_2013_03_CurAll+LANL_Sequences_examlbs100"
	insignat	<- "Sat_Jun_16_17/23/46_2013"
	infile		<- "ATHENA_2013_03_NoDRAll+LANL_Sequences_examlbs100"			
	insignat	<- "Thu_Aug_01_17/05/23_2013"
	indircov	<- paste(DATA,"derived",sep='/')
	infilecov	<- "ATHENA_2013_03_AllSeqPatientCovariates"							
	argv		<<- hivc.cmd.preclustering(indir, infile, insignat, indircov, infilecov, resume=resume)				 
	argv		<<- unlist(strsplit(argv,' '))
	clu.pre		<- hivc.prog.get.clustering.precompute()
	
	load(paste(indir,'/',"mrcas.R",sep=''))
	
	outfile				<- paste(infile,"preclust",sep='_')
	
	
	dist.brl			<- clu.pre$dist.brl.casc
	linked.bypatient	<- clu.pre$linked.bypatient
	unlinked			<- clu.pre$ph.unlinked
	unlinked.byspace	<- clu.pre$unlinked.byspace
	unlinked.bytime		<- clu.pre$unlinked.bytime
	unlinked.info		<- clu.pre$ph.unlinked.info
	ph					<- clu.pre$ph
	df.seqinfo			<- clu.pre$df.seqinfo
	#ph.mrca			<- mrca(ph)		#probably fastest
	thresh.brl			<- 0.096
	plot.file			<- paste(indir,'/',outfile,"_distbs_",gsub('/',':',insignat),".pdf", sep='')
	
	#prepare data.tables with mrca		
	bs.unlinkedpairs	<- lapply(unlinked.bytime, function(x)
			{						
				set(x, NULL, "query.FASTASampleCode", as.character(x[,query.FASTASampleCode]))
				set(x, NULL, "FASTASampleCode", as.character(x[,FASTASampleCode]))
				ans					<- merge(x, x[, list(mrca= ph.mrca[query.FASTASampleCode,FASTASampleCode]) ,by=FASTASampleCode],by="FASTASampleCode")
				setnames(ans, c("FASTASampleCode","query.FASTASampleCode"), c("tip2","tip1"))
				subset(ans, select=c(tip1, tip2, mrca))
			})
	bs.unlinkedpairs	<- rbindlist(bs.unlinkedpairs)
	#
	unlinked.byspace[,dummy:=seq_len(nrow(unlinked.byspace))]
	set(unlinked.byspace, NULL, "FASTASampleCode", as.character(unlinked.byspace[,FASTASampleCode]))	
	seq.indb			<- colnames(ph.mrca)[ which( substr(colnames(ph.mrca),1,2)!="TN" & substr(colnames(ph.mrca),1,8)!="PROT+P51" ) ]
	bs.unlinked.byspace	<- unlinked.byspace[,	list(tip1=FASTASampleCode,tip2=seq.indb, mrca= ph.mrca[seq.indb,FASTASampleCode]), by="dummy"]
	#
	setkey(linked.bypatient,Patient)
	bs.linked.bypatient	<- linked.bypatient[, {
				tmp					<- match(FASTASampleCode, ph$tip.label)
				ans					<- t(combn(tmp, 2 ))
				ans					<- cbind(ans, apply(ans, 1, function(z)  ph.mrca[z[1],z[2]]))
				data.table(tip1=ans[,1], tip2=ans[,2], mrca=ans[,3])
			}, by="Patient"]
	#	get bootstrap values
	ph.mrca				<- mrca(ph)
	tmp					<- hivc.phy.get.TP.and.TN.bootstrapvalues(ph,  bs.linked.bypatient, ph.mrca=ph.mrca, clu.pre$df.seqinfo, bs.unlinkedpairs=bs.unlinkedpairs, bs.unlinked.byspace=bs.unlinked.byspace, dist.brl=clu.pre$dist.brl.casc, thresh.brl=0.096, plot.file=paste(indir,'/',outfile,"_distbs_",gsub('/',':',insignat),".pdf"), verbose= 1)	
	#	further analysis of those with pairs with PosSeq.diff==0		
	if(0)
	{
		bs.linked.bypatient.eqPosSeqT		<- subset(bs.linked.bypatient, PosSeqT.diff==0)
		# compute raw genetic distance between sequences
		indir								<- paste(dir.name,"tmp",sep='/')
		infile								<- "ATHENA_2013_03_CurAll+LANL_Sequences"			
		insignat							<- "Sat_Jun_16_17/23/46_2013"
		file								<- paste(indir,'/',infile,'_',gsub('/',':',insignat),".R",sep='')			
		load(file)
		#dist.dna( rbind( seq.PROT.RT[bs.linked.bypatient.eqPosSeqT[i,tip1], ], seq.PROT.RT[bs.linked.bypatient.eqPosSeqT[i,tip2], ]), model="raw", as.matrix=1)
		#tmp		<- hivc.seq.dist( rbind( seq.PROT.RT[bs.linked.bypatient.eqPosSeqT[i,tip1], ], seq.PROT.RT[bs.linked.bypatient.eqPosSeqT[i,tip2], ]) )
		dummy	<- 0
		tmp		<- sapply(seq_len(nrow(bs.linked.bypatient.eqPosSeqT)),function(i)
				{
					1 - .C("hivc_dist_ambiguous_dna", seq.PROT.RT[bs.linked.bypatient.eqPosSeqT[i,tip1], ], seq.PROT.RT[bs.linked.bypatient.eqPosSeqT[i,tip2], ], ncol(seq.PROT.RT), dummy )[[4]]
				})
		bs.linked.bypatient.eqPosSeqT[,dist:=tmp]
		#
		#	conclusion: despite identical sequences, BS can be suprisingly low! Unclear if this improves with more BS replications
		#
		save(bs.linked.bypatient.eqPosSeqT, file="/Users/Oliver/duke/2013_HIV_NL/ATHENA_2013/notes/20130917_ATHENA_2013_03_NoDRAll+LANL_Sequences_examlbs100_poorBSvalues.R")
	}
}
######################################################################################
project.hivc.clustering<- function(dir.name= DATA)
{
	require(ape)
	require(data.table)
	require(RColorBrewer)
		
	if(0)	#plot composition of selected MSM clusters
	{
		hivc.prog.eval.clustering.bias()
	}	
	if(0)
	{
		hivc.prog.recombination.checkcandidates()
	}
	if(1)
	{
		project.hivc.clustering.compare.NoDR.to.NoRecombNoDR()
		project.hivc.seq.dataset.mDR.mRC.mSH.pLANL()
	}
	if(0)	#min brl to get a transmission cascade from brl matrix
	{
		brlmat	<- matrix( c(0,0.1,0.1, 	0.1,0,0.2,	0.1,0.2,0), ncol=3, nrow=3)
		brlmat	<- matrix( c(0,0.1,0.1,0.2, 	0.1,0,0.2,0.3,	0.1,0.2,0,0.1,	0.2,0.3,0.1,0), ncol=4, nrow=4)
		#brlmat	<- matrix( c(0,0.1,0.5,0.5, 	0.1,0,0.5,0.5,	0.5,0.5,0,0.1,	0.5,0.5,0.1,0), ncol=4, nrow=4)
		#brlmat	<- matrix( c(0,0.5,0.5, 	0.5,0,0.1,	0.5,0.1,0), ncol=3, nrow=3)
		#brlmat	<- matrix( c(0,0.5,0.1, 	0.5,0,0.5,	0.1,0.5,0), ncol=3, nrow=3)
		#brlmat	<- 0.5
		#brlmat	<- matrix(2, ncol=1, nrow=1)
		#brlmat	<- brlmat[upper.tri(brlmat)]
		brl		<- hivc.clu.min.transmission.cascade(brlmat)
		str(brl)
		stop()
	}
	if(0)	#test clustering on simple test tree
	{		
		ph	<- "(Wulfeniopsis:0.196108,(((TN_alpinus:0.459325,TN_grandiflora:0.259364)1.00:0.313204,uniflora:1.155678)1.00:0.160549,(((TP_angustibracteata:0.054609,(TN_brevituba:0.085086,TP_stolonifera:0.086001)0.76:0.035958)1.00:0.231339,(((axillare:0.017540,liukiuense:0.018503)0.96:0.038019,stenostachyum:0.049803)1.00:0.083104,virginicum:0.073686)1.00:0.103843)1.00:0.086965,(carinthiaca:0.018150,orientalis:0.019697)1.00:0.194784)1.00:0.077110)1.00:0.199516,(((((abyssinica:0.077714,glandulosa:0.063758)1.00:0.152861,((((allionii:0.067154,(morrisonicola:0.033595,officinalis:0.067266)1.00:0.055175)1.00:0.090694,(alpina:0.051894,baumgartenii:0.024152,(bellidioides:0.016996,nutans:0.063292)0.68:0.031661,urticifolia:0.032044)0.96:0.036973,aphylla:0.117223)0.67:0.033757,(((japonensis:0.018053,miqueliana:0.033676)1.00:0.160576,vandellioides:0.099761)0.69:0.036188,montana:0.050690)1.00:0.058380)1.00:0.115874,scutellata:0.232093)0.99:0.055014)1.00:0.209754,((((((acinifolia:0.112279,reuterana:0.108698)0.94:0.055829,pusilla:0.110550)1.00:0.230282,((davisii:0.053261,serpyllifolia:0.087290)0.89:0.036820,(gentianoides:0.035798,schistosa:0.038522)0.95:0.039292)1.00:0.092830)1.00:0.169662,(((anagalloides:0.018007,scardica:0.017167)1.00:0.135357,peregrina:0.120179)1.00:0.098045,beccabunga:0.069515)1.00:0.103473)1.00:0.287909,(((((((((((agrestis:0.017079,filiformis:0.018923)0.94:0.041802,ceratocarpa:0.111521)1.00:0.072991,amoena:0.229452,(((argute_serrata:0.017952,campylopoda:0.075210)0.64:0.034411,capillipes:0.022412)0.59:0.034547,biloba:0.037143)1.00:0.141513,intercedens:0.339760,((opaca:0.019779,persica:0.035744)0.94:0.038558,polita:0.036762)1.00:0.108620,rubrifolia:0.186799)1.00:0.144789,(((bombycina_11:0.033926,bombycina_bol:0.035290,cuneifolia:0.017300,jacquinii:0.054249,oltensis:0.045755,paederotae:0.051579,turrilliana:0.017117)0.85:0.049052,czerniakowskiana:0.089983)0.93:0.051111,farinosa:0.138075)1.00:0.080565)1.00:0.104525,((albiflora:0.017984,ciliata_Anna:0.032685,vandewateri:0.017610)0.97:0.045649,arguta:0.063057,(catarractae:0.022789,decora:0.049785)0.96:0.048220,((cheesemanii:0.040125,cupressoides:0.146538)1.00:0.067761,macrantha:0.038130)1.00:0.088158,(densifolia:0.090044,formosa:0.116180)0.71:0.046353,(elliptica:0.038650,(odora:0.019325,salicornioides:0.021228)0.94:0.042950,salicifolia:0.020829)0.92:0.043978,(nivea:0.070429,(papuana:0.035003,tubata:0.031140)0.98:0.064379)0.93:0.065336,raoulii:0.109101)0.97:0.076607)0.93:0.085835,chamaepithyoides:0.485601)0.57:0.072713,(ciliata_157:0.069943,lanuginosa:0.052833)1.00:0.098638,(densiflora:0.069429,macrostemon:0.118926)0.92:0.124911,(fruticulosa:0.086891,saturejoides:0.041181)0.94:0.086148,kellererii:0.083762,lanosa:0.263033,mampodrensis:0.103384,nummularia:0.191180,pontica:0.128944,thessalica:0.129197)0.65:0.031006,(arvensis:0.342138,(((((chamaedrys:0.043720,micans:0.032021,vindobonensis:0.033309)0.51:0.034053,micrantha:0.019084)0.64:0.037906,krumovii:0.020175)1.00:0.103875,verna:0.254017)0.81:0.057105,magna:0.112657)1.00:0.104070)1.00:0.101845)1.00:0.149208,(((aznavourii:0.664103,glauca:0.405588)0.85:0.209945,praecox:0.447238)1.00:0.185614,(donii:0.260827,triphyllos:0.176032)1.00:0.194928)1.00:0.611079)0.74:0.055152,((crista:0.591702,(((cymbalaria_Avlan:0.017401,panormitana:0.017609)1.00:0.229508,((cymbalaria_Istanbul:0.028379,trichadena_332:0.016891,trichadena_Mugla:0.019131)1.00:0.196417,lycica_333:0.146772)1.00:0.097646,lycica_192:0.154877)1.00:0.234748,(((hederifolia:0.018068,triloba:0.075784)1.00:0.084865,(sibthorpioides:0.122542,sublobata:0.136951)1.00:0.074683)0.89:0.043623,stewartii:0.040679)1.00:0.596859)1.00:0.237324)0.58:0.057120,javanica:0.133802)1.00:0.137214)1.00:0.269201,(missurica:0.016685,rubra:0.019696)1.00:0.351184)0.54:0.058275)0.52:0.062485,((dahurica:0.023542,longifolia:0.016484,spicata:0.018125)0.95:0.042294,(nakaiana:0.016270,schmidtiana:0.058451)0.88:0.037207)1.00:0.261643)0.55:0.056458)1.00:0.229509,kurrooa:0.100611)0.74:0.068198,(bonarota:0.040842,lutea:0.115316)1.00:0.241657)0.99:0.085772);"
		ph <- ladderize( read.tree(text = ph) )				
		#read bootstrap support values
		thresh.bs						<- 0.9
		ph.node.bs						<- as.numeric( ph$node.label )		
		ph.node.bs[is.na(ph.node.bs)]	<- 0
		ph.node.bs[c(13,15,27,41,43)]	<- thresh.bs-0.05*seq(0.01,length.out=5)
		ph.node.bs[ph.node.bs==1]		<- seq_along(which(ph.node.bs==1))*0.005 + 0.7
		ph$node.label					<- ph.node.bs
		#read patristic distances
		stat.fun						<- hivc.clu.min.transmission.cascade
		#stat.fun						<- max
		dist.brl						<- hivc.clu.brdist.stats(ph, eval.dist.btw="leaf", stat.fun=stat.fun)
		print(dist.brl)
		thresh.brl						<- quantile(dist.brl,seq(0.1,1,by=0.05))["100%"]
		print(quantile(dist.brl,seq(0.1,0.5,by=0.05)))
		print(thresh.brl)
		#produce clustering 
		clustering	<- hivc.clu.clusterbythresh(ph, thresh.nodesupport=thresh.bs, thresh.brl=thresh.brl, dist.brl=dist.brl, nodesupport=ph.node.bs,retval="all")
		print(clustering)		
		hivc.clu.plot(ph, clustering[["clu.mem"]], highlight.edge.of.tiplabel=c("TN_","TP_"), highlight.edge.of.tiplabel.col= c("red","blue") )
		#produce some tip states
		ph.tip.state<- rep(1:20, each=ceiling( length( ph$tip.label )/20 ))[seq_len(Ntip(ph))]
		states		<- data.table(state.text=1:20, state.col=rainbow(20))
		states		<- states[ph.tip.state,]
		hivc.clu.plot.tiplabels( seq_len(Ntip(ph)), matrix(states[,state.text], nrow=1,ncol=nrow(states)), matrix(states[,state.col], nrow=1,ncol=nrow(states)), cex=0.4, adj=c(-1,0.5) )
				
		#get nodes in each cluster
		#clu.nodes		<- sapply(unique(clustering[["clu.mem"]])[-1],function(clu){		which(clustering[["clu.mem"]]==clu)	})
		#get boostrap support of index cases
		clu.idx			<- clustering[["clu.idx"]]-Ntip(ph)
		print(clu.idx)
		print(ph.node.bs[clu.idx])
		stop()
	}
	if(0)
	{
		#plot properties of concentrated clusters by region
		verbose							<- 1
		thresh.bs						<- 0.85
		thresh.brl						<- 0.06				
		
		infile		<- paste("ATHENA_2013_03_CurAll+LANL_Sequences_examlbs100_withinpatientclusters_bs",thresh.bs*100,"_brlmed",thresh.brl*100,"_clusterinfo",sep='')
		insignat	<- "Sat_Jun_16_17/23/46_2013"
		file		<- paste(dir.name,"tmp",paste(infile,'_',gsub('/',':',insignat),".R",sep=''),sep='/')
		if(verbose)	cat(paste("read cluster info from file",file))
		load(file)
		if(verbose)	cat(paste("found seq, n=",nrow(df.seqinfo)))		
		print(df.seqinfo)
		df.seqinfo	<- subset(df.seqinfo, !is.na(cluster) )
		if(verbose)	cat(paste("found seq in clu, n=",nrow(df.seqinfo)))
		
		#remove within patient clusters before out-of-NL taken out		
		set(df.seqinfo, which( df.seqinfo[,is.na(Patient)] ), "Patient", "NA")
		tmp			<- df.seqinfo[, list(isWithinCluster= length(unique(Patient))==1), by=cluster]
		df.seqinfo	<- merge( subset(tmp,!isWithinCluster), df.seqinfo, by="cluster" )
		#from between patient clusters, remove within patient seqs and foreign seqs
		df.pat		<- df.seqinfo[,list(CountryInfection=CountryInfection[1], Trm=Trm[1], Sex=Sex[1], NegT=NegT[1], AnyPos_T1=AnyPos_T1[1], RegionHospital=RegionHospital[1]), by=Patient]			
		df.clu		<- df.seqinfo[, list(Patient=unique(Patient)), by=cluster]		
		df.clu		<- merge(df.clu, df.pat, all.x=1, by="Patient")
		
		
		#determine if cluster concentrated in region, and where
		setkey(df.clu, cluster)
		if(verbose)	cat(paste("found pat in clu, n=",nrow(df.clu)))
		clu.reg		<- table( df.clu[,cluster,RegionHospital] )
		clu.reg		<- clu.reg / matrix(apply(clu.reg, 2, sum), nrow=nrow(clu.reg), ncol=ncol(clu.reg), byrow=1)				
		clu.regconc	<- apply( clu.reg, 2, function(x)	any(x>0.4) && length(which(x!=0))>2 ) |
					   apply( clu.reg, 2, function(x)	any(x>0.5) && length(which(x!=0))<=2 )
		if(verbose)	cat(paste("number of clusters, n=",length(clu.regconc)))
		if(verbose)	cat(paste("number of spatially concentrated clusters, n=",length(which(clu.regconc))))
		tmp			<- rownames(clu.reg)
		clu.regconc2<- sapply( seq_len(ncol(clu.reg)), function(j)
						{
							ifelse(!clu.regconc[j], NA, tmp[ which.max(clu.reg[,j]) ] )						
						})
		df.clureg	<- data.table(cluster= as.numeric(colnames( clu.reg )), IsRegionConcentrated= clu.regconc, RegionConcentrated=clu.regconc2, key="cluster" )
		
		#determine median AnyPos_T1 time for cluster and add
		tmp			<- df.clu[,list(cluster.PosT=mean(AnyPos_T1, na.rm=T)),by=cluster]				
		df.clureg	<- merge(df.clureg, tmp, all.x=1, by="cluster")

		#merge all
		df.seqinfo	<- merge( df.seqinfo, df.clureg, all.x=1, by="cluster" )
						
		
		#plot Bezemer style
		file				<- paste(dir.name,"tmp",paste(infile,"_spatialvariation_",gsub('/',':',insignat),".pdf",sep=''),sep='/')
		df					<- df.seqinfo
		df[,time:=AnyPos_T1 ]		
		set(df, which(df[,!IsRegionConcentrated]), "RegionConcentrated", "Bridge")
		tmp					<- numeric(nrow(df))
		tmp2				<- c("N","E","S","W","Amst","Bridge")
		for( i in seq_along(tmp2))
			tmp[which(df[,RegionConcentrated==tmp2[i]])]	<- i
		df[,sort1:=factor(tmp, levels=1:6, labels=tmp2) ]
		df[,sort2:=cluster.PosT]		
		tmp					<- numeric(nrow(df))		
		tmp2				<- c("N","E","S","W","Amst")
		for( i in seq_along(tmp2))
			tmp[which(df[,RegionHospital==tmp2[i]])]	<- i
		df[,covariate:=factor(tmp, levels=1:5, labels=tmp2) ]
		xlab				<- "first HIV+ event"
		ylab				<- paste("clusters BS=",thresh.bs*100," BRL.med=",thresh.brl*100," version=",gsub('/',':',insignat),sep="")
		cex.points			<- 0.5
		
		#extract clusters one by one in desired order
		tmp					<- df[,list(sort1=sort1[1], sort2=sort2[1]),by=cluster]
		clusters.sortby1	<- as.numeric( tmp[,sort1] )
		clusters.sortby2	<- as.numeric( tmp[,sort2] )
		clusters.sortby		<- order(clusters.sortby1,clusters.sortby2, decreasing=F)
		clusters			<- tmp[clusters.sortby,cluster]								
		clusters.levels		<- levels(df[,covariate])					
		clusters			<- lapply(clusters,function(x)	subset(df,cluster==x,c(time,covariate))		)
		
		ylim				<- c(-2,length(clusters))
		cols				<- brewer.pal(length(clusters.levels),"Set1")	#[c(2,1,3,4,5)]		
#cols			<- c("green","red","blue","grey30")
		cols				<- sapply(cols, function(x) my.fade.col(x,0.7))
		pch					<- c(rep(19,length(clusters.levels)))
		xlim				<- range( df[,time], na.rm=1 )
		xlim[1]				<- xlim[1]-700
		
		cat(paste("\nwrite plot to",file))
		pdf(file,width=7,height=14)
		plot(1,1,type='n',bty='n',xlim=xlim,ylim=ylim,ylab=ylab, xaxt='n',yaxt='n',xlab=xlab)
		axis.Date(1, seq.Date(xlim[1],xlim[2],by="year"), labels = TRUE )
		dummy<- sapply(seq_along(clusters),function(i)
				{							
					cluster.cex	<- cex.points	#(1 - nrow(subset(clusters[[i]],hregion=="unknown")) / nrow(clusters[[i]])) * cex.points
					cluster.ix	<- order(as.numeric(clusters[[i]][,time]))
					cluster.x	<- as.numeric(clusters[[i]][,time])[cluster.ix]										
					cluster.y	<- rep(i,nrow(clusters[[i]]))
					cluster.z	<- as.numeric(clusters[[i]][,covariate])[cluster.ix]
					#print(clusters[[i]][,covariate])
					#if(cluster.z[1]!=4)
					#lines( c(cluster.x[1],xlim[2]), rep(i,2), col=cols[ cluster.z[1] ] )					
					points( cluster.x, cluster.y, col=cols[ cluster.z ], pch=pch[1], cex=cluster.cex )										
				})	
		pch[pch==19]	<- 21
		legend("topleft",bty='n',pt.bg=cols,pch=pch,legend=levels(df[,covariate]), col=rep("transparent",length(clusters.levels)))				
		dev.off()	
		stop()
	}
	if(0)
	{
		#check BEEHIVE sequences
		indircov	<- paste(DATA,"derived",sep='/')
		infilecov	<- "ATHENA_2013_03_AllSeqPatientCovariates"							
		load(paste(indircov,'/',infilecov,".R",sep=''))
		
		df.bee		<- as.data.table(read.csv("~/duke/2013_HIV_NL/ATHENA_2013/data/BEEHIVE1_data.csv", stringsAsFactors=0))
		setnames(df.bee,"mcode","Patient")
		df.bee		<- merge(df.all, subset(df.bee,select=Patient), by="Patient")
		save(df.bee, file="~/duke/2013_HIV_NL/ATHENA_2013/data/BEEHIVE1_data.R")
		
	}
	if(0)
	{		
		verbose		<- 1
		resume		<- 1
		#
		# precompute clustering stuff		
		#
		patient.n	<- 15700
		indir		<- paste(DATA,"tmp",sep='/')		
		infile		<- "ATHENA_2013_03_CurAll+LANL_Sequences_examlbs100"
		insignat	<- "Sat_Jun_16_17/23/46_2013"
		infile		<- "ATHENA_2013_03_NoDRAll+LANL_Sequences_examlbs100"			
		insignat	<- "Thu_Aug_01_17/05/23_2013"
		indircov	<- paste(DATA,"derived",sep='/')
		infilecov	<- "ATHENA_2013_03_AllSeqPatientCovariates"							
		argv		<<- hivc.cmd.preclustering(indir, infile, insignat, indircov, infilecov, resume=resume)				 
		argv		<<- unlist(strsplit(argv,' '))
		clu.pre		<- hivc.prog.get.clustering.precompute()
		#
		# evaluate TPTN for various thresholds
		#
		if(verbose) cat(paste("compute TPTN for dist.brl.casc"))
		argv			<<- hivc.cmd.clustering.tptn(indir, infile, insignat, indircov, infilecov, opt.brl="dist.brl.casc", patient.n=patient.n, resume=resume)
		argv			<<- unlist(strsplit(argv,' '))
		clu.tptn.casc	<- hivc.prog.get.clustering.TPTN(clu.pre=clu.pre)
		if(verbose) cat(paste("compute TPTN for dist.brl.max"))
		argv			<<- hivc.cmd.clustering.tptn(indir, infile, insignat, indircov, infilecov, opt.brl="dist.brl.max", patient.n=patient.n, resume=resume)
		argv			<<- unlist(strsplit(argv,' '))
		clu.tptn.max	<- hivc.prog.get.clustering.TPTN(clu.pre=clu.pre)
		#
		# compare dist.brl vs dist.max in terms of total patients in clusters
		#
		file			<- paste(indir, paste(infile,"_clust_cascvsmax_patients_",gsub('/',':',outsignat),".pdf",sep=''),sep='/')
		pdf(file=file,width=5,height=5)
		par(mar=c(4,4,1,0))
		dbwpat.cmp		<- list(	max	=lapply(seq_along(select.bs), function(i)	clu.tptn.max[["clusters.dbwpat"]][[select.brl[i]]][[select.bs[i]]]		),
									casc=lapply(seq_along(select.bs), function(i)	clu.tptn.casc[["clusters.dbwpat"]][[select.brl[i]]][[select.bs[i]]]		)	)
		ylim			<- c(0,3100)
		#sapply(dbwpat.cmp, function(x)	sapply(x, function(z) sum(as.numeric(names(z))*z)  ) )
		xlim			<- c(1,max( sapply(dbwpat.cmp, function(x)	sapply(x, function(z) max(as.numeric(names(z)))) ) ))
		
		
		plot(1,1,type='n',bty='n',xlim=xlim,ylim=ylim, xlab="cluster size", ylab="patients in clusters of size <=x")
		lapply(seq_along(dbwpat.cmp), function(i)
				{
					lapply(seq_along(dbwpat.cmp[[i]]),function(j)
							{	
								z	<- dbwpat.cmp[[i]][[j]][-1]
								z2	<- as.numeric(names(z))								
								lines(z2,cumsum(z*z2), col=cols[i], lty=j)
							})
				})
		legend("topleft",bty='n', border=NA, fill= cols, legend= c("max","casc"))
		legend("topright",bty='n', border=NA, lty= 1:2, legend= c("BS=0.8, BRL=0.1","BS=0.95, BRL=0.04"))
		dev.off()
		#
		# compare dist.brl vs dist.max in terms of %cov 
		#
		file			<- paste(indir, paste(infile,"_clust_cascvsmax_covepi_",gsub('/',':',outsignat),".pdf",sep=''),sep='/')
		pdf(file=file,width=5,height=5)
		par(mar=c(4,4,1,0))
		ylim			<- c(0,0.2)				
		plot(1,1,type='n',bty='n',xlim=xlim,ylim=ylim, xlab="cluster size", ylab="%coverage (of epi) in clusters of size <=x")
		lapply(seq_along(dbwpat.cmp), function(i)
				{
					lapply(seq_along(dbwpat.cmp[[i]]),function(j)
							{	
								z	<- dbwpat.cmp[[i]][[j]][-1]
								z2	<- as.numeric(names(z))								
								lines(z2,cumsum(z*z2)/patient.n, col=cols[i], lty=j)
							})
				})		
		legend("topleft",bty='n', border=NA, fill= cols, legend= c("max","casc"))
		legend("topright",bty='n', border=NA, lty= 1:2, legend= c("BS=0.8, BRL=0.1","BS=0.95, BRL=0.04"))
		dev.off()
		#
		# compare dist.brl vs dist.max in terms of %cov vs %fp for BS=0.8/BRL=0.1/casc  vs BS=0.95/BRL=0.05/max 
		#
		covfp.cmp.x	<- rbind( clu.tptn.casc[["fp.by.all"]]["0.8",], clu.tptn.max[["fp.by.all"]]["0.95",]	)
		covfp.cmp.y	<- rbind( clu.tptn.casc[["clusters.cov.epidemic"]]["0.8",], clu.tptn.max[["clusters.cov.epidemic"]]["0.95",]	)		
		
		file			<- paste(indir, paste(infile,"_clust_cascvsmax_covepifp_",gsub('/',':',outsignat),".pdf",sep=''),sep='/')
		pdf(file=file,width=5,height=5)
		par(mar=c(4,4,1,0))		
		plot(1,1,type='n',bty='n',xlim=range(c(0.01,covfp.cmp.x)),ylim=range(c(0.2,covfp.cmp.y)),xlab="%FP (among all)",ylab="%coverage (of epi)")
		dummy	<- sapply(seq_len(nrow(covfp.cmp.x)),function(i)
				{					
					points(covfp.cmp.x[i,],covfp.cmp.y[i,],col=cols[i],type='b')
					text(covfp.cmp.x[i,],covfp.cmp.y[i,],col=cols[i],labels=as.numeric(colnames(covfp.cmp.y)),adj=c(-0.8,0.5),cex=0.5)
				})
		legend("bottomright",border=NA,bty='n',fill=cols,legend=c("max","casc"))
		dev.off()				
	}
	if(0)
	{
		#get branch lengths between F2F transmissions, where there must be a missed intermediary
		outdir					<- indir
		outfile					<- paste(infile,"_clust_",opt.brl,"_bs",thresh.bs*100,"_brl",thresh.brl*100,'_',"hetsamegender",sep='')
		outsignat				<- insignat		
		plot.file				<- paste(outdir,'/',outfile,'_',gsub('/',':',outsignat),".pdf",sep='')		
		clu.female2female		<- hivc.clu.getplot.female2female( clu.pre$ph, clu$clustering, df.cluinfo, plot.file=paste(outdir,'/',outfile,'_',gsub('/',':',outsignat),".pdf",sep='') )				
		clu.missedintermediary	<- names(cluphy.brl.bwpat) %in% as.character( unique( clu.female2female$cluphy.df[,cluster] ) )
		brl.missedintermediary	<- unlist( lapply( which(clu.missedintermediary), function(i)	cluphy.brl.bwpat[[i]] ) )
		brl.others				<- na.omit( unlist( lapply( which(!clu.missedintermediary), function(i)	cluphy.brl.bwpat[[i]] ) ) )
		brl.breaks				<- seq(0,max(brl.others,brl.missedintermediary)*1.1,by=0.02)
		brl.cols				<- sapply( brewer.pal(3, "Set1"), function(x) my.fade.col(x,0.5) )
		
		par(mfcol=c(1,2))
		hist( brl.others, breaks=brl.breaks, col=brl.cols[1], border=NA, freq=T, main='', xlab="branch lengths, others"  )
		#legend("topright", fill= brl.cols[1:2], legend= c("others","missed intermediary"), bty='n', border=NA)
		hist( brl.missedintermediary, breaks=brl.breaks, col=brl.cols[2], border=NA, freq=T, main='', xlab="branch lengths, missed intermediary" )		
		par(mfcol=c(1,1))
		
		
		#highlight seq of same patient not in cluster
		#select clusters with TN and P51

		
		stop()
		
	}
	if(0)	#count how many unlinked pairs in clustering
	{		
		verbose	<- 1
		file	<- paste(dir.name,"derived/ATHENA_2013_03_Unlinked_SeroConv_Dead_UnlinkedAll.R", sep='/')		
		if(verbose)	cat(paste("read file",file))
		load(file)
		#str(unlinked.bytime)
		
		infile		<- "ATHENA_2013_03_FirstCurSequences_PROTRT_examlbs100"
		signat.in	<- "Fri_May_24_12/59/06_2013"
		signat.out	<- "Fri_Jun_07_09/59/23_2013"
		file		<- paste(dir.name,"tmp",paste(infile,'_',gsub('/',':',signat.in),".newick",sep=''),sep='/')
		cat(paste("read file",file))
		
		#
		#set up tree, get boostrap values and patristic distances between leaves
		#
		ph								<- ladderize( read.tree(file) )		
		ph.node.bs						<- as.numeric( ph$node.label )
		ph.node.bs[is.na(ph.node.bs)]	<- 0
		ph.node.bs						<- ph.node.bs/100
		ph$node.label					<- ph.node.bs
		#print(quantile(ph.node.bs,seq(0.1,1,by=0.1)))		
		dist.brl						<- hivc.clu.brdist.stats(ph, eval.dist.btw="leaf")		#read patristic distances -- this is the expensive step but still does not take very long
		#print(quantile(dist.brl,seq(0.1,1,by=0.1)))
		
		#
		#convert truly unlinked pairs from Patient name to ph node index
		#
		df.tips							<- data.table(Node=seq_len(Ntip(ph)), Patient=ph$tip.label )
		setkey(df.tips, Patient)																#use data.table to speed up search
		df.tips							<- df.all[df.tips]
		ph.unlinked.dead				<- lapply(seq_along(unlinked.bytime), function(j)
											{
												as.numeric( sapply(seq_along(unlinked.bytime[[j]]), function(i)	subset(df.tips, unlinked.bytime[[j]][i]==Patient, Node) ))					
											})		
		ph.unlinked.seroneg				<- as.numeric( sapply(names(unlinked.bytime), function(x)	subset(df.tips, x==Patient, Node) ))		
		tmp								<- sort( ph.unlinked.seroneg, index.return=1)			#sort ph.unlinked.seroneg and return index, so we can also sort ph.unlinked.dead
		ph.unlinked.seroneg				<- tmp$x
		names(ph.unlinked.seroneg)		<- names(unlinked.bytime)[tmp$ix]
		ph.unlinked.dead				<- lapply(tmp$ix, function(i){		ph.unlinked.dead[[i]]	})	
		names(ph.unlinked.dead)			<- names(unlinked.bytime)[tmp$ix]		
		ph.unlinked.seroneg				<- data.table(PhNode=ph.unlinked.seroneg, PhNodeUnlinked= seq_along(ph.unlinked.seroneg), Patient= names(ph.unlinked.seroneg))
		setkey(ph.unlinked.seroneg, Patient)													#set key to take right outer join with df.serocon -- temporarily destroy correspondance with 'ph.unlinked.dead'
		setkey(df.serocon, Patient)
		ph.unlinked.seroneg				<- df.serocon[ph.unlinked.seroneg]		
		setkey(ph.unlinked.seroneg, PhNode)			#use data.table to speed up search -- restore correspondance with 'ph.unlinked.dead'		
		#print(ph.unlinked.seroneg); print(ph.unlinked.dead); stop()

		#
		#determine threshold for patristic distances based on truly unlinked pairs
		#
		thresh.brl	<- NULL
		thresh.bs	<- 0.9
		ans			<- hivc.clu.clusterbytypeIerror(ph, dist.brl, ph.node.bs, ph.unlinked.seroneg, ph.unlinked.dead, thresh.brl=thresh.brl, thresh.bs=thresh.bs, level= 0.001, verbose=1)		
		#
		#determine quick estimate of expected false pos for these thresholds if clustering was random
		#				
		check		<- hivc.clu.exp.typeIerror.randomclu(ph, dist.brl, ph.node.bs, ph.unlinked.seroneg, ph.unlinked.dead, ans[["thresh.brl"]], ans[["thresh.bs"]])
		cat(paste("\n#clusters with at least one FP=",check[["fp.n"]]," , %clusters with at least one FP=",check[["fp.rate"]],sep=''))
		#
		#max cluster size distribution
		barplot(cumsum(h$counts), names.arg=seq_along(h$counts)+1, axisnames=1, xlab="maximum cluster size")
		#
		#number of sequences in cluster, and %
		cat(paste("\n#seq in cluster=",sum(ans[["clu"]][["size.tips"]]),", %seq in cluster=",sum(ans[["clu"]][["size.tips"]])/Ntip(ph),sep='')) 		
		#
		#plot final clusters, save  clusters
		#								
		file		<- paste(dir.name,"tmp",paste(infile,"_clu_",gsub('/',':',signat.out),".pdf",sep=''),sep='/')
		cat(paste("plot clusters on tree to file",file))
		hivc.clu.plot(ph, ans[["clu"]][["clu.mem"]], file=file, pdf.scaley=25)
		file		<- paste(dir.name,"tmp",paste(infile,"_clu_",gsub('/',':',signat.out),".R",sep=''),sep='/')
		cat(paste("save analysis to file",file))
		save(ans,check,ph,dist.brl,ph.node.bs,ph.unlinked.seroneg,ph.unlinked.dead, file=file)				
	}
}
######################################################################################
project.hivc.clustering.selectparticularclusters<- function()
{	
	if(1)
	{
		verbose		<- 1
		resume		<- 1
		indir		<- paste(DATA,"tmp",sep='/')		
		infile		<- "ATHENA_2013_03_NoDRAll+LANL_Sequences_examlbs100"			
		insignat	<- "Thu_Aug_01_17/05/23_2013"
		indircov	<- paste(DATA,"derived",sep='/')
		infilecov	<- "ATHENA_2013_03_AllSeqPatientCovariates"
		opt.brl		<- "dist.brl.casc" 
		thresh.brl	<- 0.096
		thresh.bs	<- 0.8
		argv		<<- hivc.cmd.preclustering(indir, infile, insignat, indircov, infilecov, resume=resume)				 
		argv		<<- unlist(strsplit(argv,' '))
		clu.pre		<- hivc.prog.get.clustering.precompute()
		argv		<<- hivc.cmd.clustering(indir, infile, insignat, opt.brl, thresh.brl, thresh.bs, resume=resume)				 
		argv		<<- unlist(strsplit(argv,' '))
		clu			<- hivc.prog.get.clustering()		
		#
		# remove singletons
		#
		if(verbose) cat(paste("\nnumber of seq in tree is n=", nrow(clu$df.cluinfo)))
		df.cluinfo	<- subset(clu$df.seqinfo, !is.na(cluster) )
		if(verbose) cat(paste("\nnumber of seq in clusters is n=", nrow(df.cluinfo)))
		if(verbose) cat(paste("\nnumber of clusters is n=", length(unique(df.cluinfo[,cluster]))))
		#
		# remove within patient clusters
		#
		tmp			<- subset(df.cluinfo[,list(clu.is.bwpat=length(unique(Patient))>1),by="cluster"], clu.is.bwpat, cluster )
		df.cluinfo	<- merge(tmp, df.cluinfo, by="cluster", all.x=1)
		#
		# plot merged clusters that have shared patients
		#
		outdir			<- indir
		outfile			<- paste(infile,"_clust_",opt.brl,"_bs",thresh.bs*100,"_brl",thresh.brl*100,'_',"sharingpatclu2",sep='')
		outsignat		<- insignat									
		tmp				<- hivc.clu.getplot.potentialsuperinfections(clu.pre$ph, clu$clustering, df.cluinfo, plot.file=paste(outdir,'/',outfile,'_',gsub('/',':',outsignat),".pdf",sep='') )		 		
		# identify potential superinfections by name
		cluphy.df		<- tmp$cluphy.
		cluphy.df		<- subset(cluphy.df, select=c(cluster, FASTASampleCode, Patient, PosSeqT))
		tmp				<- cluphy.df[, list(count= length(unique(cluster)), cluster=cluster, FASTASampleCode=FASTASampleCode, PosSeqT=PosSeqT),by="Patient"]
		tmp				<- subset(tmp,count>1)
		unique(tmp[,Patient])
	}
	if(1)
	{	
		msm<- hivc.prog.get.clustering.MSM()		
		#
		#1) plot clusters with 		small brl / npat 	--> explosive for targeted testing -- mostly acute -- serial or starlike or what ?
		#
		df.cluinfo				<- msm$df.cluinfo
		tmp						<- df.cluinfo[,list(clu.bwpat.medbrl=clu.bwpat.medbrl[1],clu.npat=clu.npat[1], select=clu.bwpat.medbrl[1]/clu.npat[1]),by="cluster"]
		cumsum( table( tmp[,clu.npat] ) / nrow(tmp) )
		tmp						<- subset(tmp, clu.npat>4)
		tmp						<- subset( tmp, select<quantile( tmp[,select], probs=0.2 ))
		cluphy.df				<- merge( subset(tmp,select=cluster), df.cluinfo, all.x=1, by="cluster" )
		if(verbose) cat(paste("\nnumber of selected sequences is n=",nrow(cluphy.df)))		
		cluphy.subtrees			<- lapply( as.character(tmp[,cluster]), function(x)  msm$cluphy.subtrees[[x]]	)
		names(cluphy.subtrees)	<- as.character(tmp[,cluster])
		outdir					<- indir
		outfile					<- paste(infile,"_clust_",opt.brl,"_bs",thresh.bs*100,"_brl",thresh.brl*100,'_',"msmexpgr_selectexplosive",sep='')
		outsignat				<- insignat									
		tmp						<- hivc.clu.polyphyletic.clusters(cluphy.df, cluphy.subtrees=cluphy.subtrees, plot.file=paste(outdir,'/',outfile,'_',gsub('/',':',outsignat),".pdf",sep=''), pdf.scaley=10, adj.tiplabel= c(-0.05,0.5), cex.tiplabel=0.3, pdf.xlim=0.36)
		cluphy					<- tmp$cluphy		
		cluphy.df				<- subset(cluphy.df, select=c(cluster, FASTASampleCode, Patient,    PosSeqT,   DateBorn, Sex,  CountryBorn, CountryInfection, RegionHospital))			
		outfile					<- paste(DATA,'/',"tmp",'/',infile,"_clust_",opt.brl,"_bs",thresh.bs*100,"_brl",thresh.brl*100,'_',"msmexpgr_selectexplosive.csv",sep='')
		write.table(cluphy.df, file=outfile, sep=",", row.names=F)
		outfile					<- paste(DATA,'/',"tmp",'/',infile,"_clust_",opt.brl,"_bs",thresh.bs*100,"_brl",thresh.brl*100,'_',"msmexpgr_selectexplosive.R",sep='')
		save(cluphy.df, file=outfile)
		#
		#2) plot clusters with 		npat>4 	very small acute	small brl / npat 	--> explosive for targeted testing -- non - acute
		#
		df.cluinfo				<- msm$df.cluinfo
		tmp						<- df.cluinfo[,list(clu.bwpat.medbrl=clu.bwpat.medbrl[1],clu.npat=clu.npat[1], clu.fPossAcute=clu.fPossAcute[1], select=clu.bwpat.medbrl[1]/clu.npat[1]/clu.npat[1]),by="cluster"]
		cumsum( table( tmp[,clu.npat] ) / nrow(tmp) )
		tmp						<- subset(tmp, clu.npat>4)
		tmp						<- subset( tmp,  clu.fPossAcute<quantile( tmp[,clu.fPossAcute], probs=0.2 ) )
		cluphy.df				<- merge( subset(tmp,select=cluster), df.cluinfo, all.x=1, by="cluster" )
		if(verbose) cat(paste("\nnumber of selected sequences is n=",nrow(cluphy.df)))		
		cluphy.subtrees			<- lapply( as.character(tmp[,cluster]), function(x)  msm$cluphy.subtrees[[x]]	)
		names(cluphy.subtrees)	<- as.character(tmp[,cluster])
		outdir					<- indir
		outfile					<- paste(infile,"_clust_",opt.brl,"_bs",thresh.bs*100,"_brl",thresh.brl*100,'_',"msmexpgr_selectlargenonacute",sep='')
		outsignat				<- insignat									
		tmp						<- hivc.clu.polyphyletic.clusters(cluphy.df, cluphy.subtrees=cluphy.subtrees, plot.file=paste(outdir,'/',outfile,'_',gsub('/',':',outsignat),".pdf",sep=''), pdf.scaley=10, adj.tiplabel= c(-0.05,0.5), cex.tiplabel=0.3, pdf.xlim=0.36)
		#
		#3) plot clusters with 		high VLI after treat			--> likely to infect -- check manually ?	large clusters are a subset of (2) so plot all
		#
		df.cluinfo				<- msm$df.cluinfo
		tmp						<- df.cluinfo[,	{
					z			<- which(PosSeqT>AnyT_T1)
					hlRNA.mx	<- ifelse(length(z), max(lRNA_aTS[z]), 0)
					hlRNA.sm	<- ifelse(length(z), sum(lRNA_aTS[z]), 0)
					list(clu.bwpat.medbrl=clu.bwpat.medbrl[1],clu.npat=clu.npat[1], clu.fPossAcute=clu.fPossAcute[1], hlRNA.mx=hlRNA.mx, hlRNA.sm=hlRNA.sm)
				},by="cluster"]
		tmp						<- subset( tmp, hlRNA.sm>=quantile(tmp[,hlRNA.sm], probs=0.8) )
		cluphy.df				<- merge( subset(tmp,select=cluster), df.cluinfo, all.x=1, by="cluster" )
		if(verbose) cat(paste("\nnumber of selected sequences is n=",nrow(cluphy.df)))		
		cluphy.subtrees			<- lapply( as.character(tmp[,cluster]), function(x)  msm$cluphy.subtrees[[x]]	)
		names(cluphy.subtrees)	<- as.character(tmp[,cluster])
		outdir					<- indir
		outfile					<- paste(infile,"_clust_",opt.brl,"_bs",thresh.bs*100,"_brl",thresh.brl*100,'_',"msmexpgr_selecthighVLduringtreat",sep='')
		outsignat				<- insignat									
		tmp						<- hivc.clu.polyphyletic.clusters(cluphy.df, cluphy.subtrees=cluphy.subtrees, plot.file=paste(outdir,'/',outfile,'_',gsub('/',':',outsignat),".pdf",sep=''), pdf.scaley=15, adj.tiplabel= c(-0.05,0.5), cex.tiplabel=0.25, pdf.xlim=0.36)
		#
		#4) plot clusters with 		long treatment interruptions			--> likely to infect -- check manually ?	large clusters are a subset of (2) so plot all
		#
		df.cluinfo				<- msm$df.cluinfo
		tmp						<- df.cluinfo[,	{
					z		<- which(!is.na(AnyT_T1))
					TrI.mx	<- ifelse(length(z),	max(TrImo_bTS[z] + TrImo_aTS[z]), 0)
					TrI.sm	<- ifelse(length(z),	sum(TrImo_bTS[z] + TrImo_aTS[z]), 0)
					list(clu.bwpat.medbrl=clu.bwpat.medbrl[1],clu.npat=clu.npat[1], clu.fPossAcute=clu.fPossAcute[1], TrI.mx=TrI.mx, TrI.sm=TrI.sm)
				},by="cluster"]
		tmp						<- subset(tmp, TrI.mx>=quantile(tmp[,TrI.mx], probs=0.8)  & clu.fPossAcute<0.7 )
		cluphy.df				<- merge( subset(tmp,select=cluster), df.cluinfo, all.x=1, by="cluster" )
		if(verbose) cat(paste("\nnumber of selected sequences is n=",nrow(cluphy.df)))		
		cluphy.subtrees			<- lapply( as.character(tmp[,cluster]), function(x)  msm$cluphy.subtrees[[x]]	)
		names(cluphy.subtrees)	<- as.character(tmp[,cluster])
		outdir					<- indir
		outfile					<- paste(infile,"_clust_",opt.brl,"_bs",thresh.bs*100,"_brl",thresh.brl*100,'_',"msmexpgr_selectLongTRI",sep='')
		outsignat				<- insignat									
		tmp						<- hivc.clu.polyphyletic.clusters(cluphy.df, cluphy.subtrees=cluphy.subtrees, plot.file=paste(outdir,'/',outfile,'_',gsub('/',':',outsignat),".pdf",sep=''), pdf.scaley=12, adj.tiplabel= c(-0.05,0.5), cex.tiplabel=0.27, pdf.xlim=0.36)
		#
		#5) plot clusters with 		NegT			--> might help to better understand what is going on ?	
		#
		df.cluinfo				<- msm$df.cluinfo
		tmp						<- df.cluinfo[,	list(clu.bwpat.medbrl=clu.bwpat.medbrl[1],clu.npat=clu.npat[1], clu.fPossAcute=clu.fPossAcute[1], fNegT=length(which(!is.na(NegT))) / clu.ntip[1]),by="cluster"]										
		tmp						<- subset(tmp, fNegT>=quantile(tmp[,fNegT], probs=0.8) )
		cluphy.df				<- merge( subset(tmp,select=cluster), df.cluinfo, all.x=1, by="cluster" )
		if(verbose) cat(paste("\nnumber of selected sequences is n=",nrow(cluphy.df)))		
		cluphy.subtrees			<- lapply( as.character(tmp[,cluster]), function(x)  msm$cluphy.subtrees[[x]]	)
		names(cluphy.subtrees)	<- as.character(tmp[,cluster])
		outdir					<- indir
		outfile					<- paste(infile,"_clust_",opt.brl,"_bs",thresh.bs*100,"_brl",thresh.brl*100,'_',"msmexpgr_selectAllWithNegT",sep='')
		outsignat				<- insignat									
		tmp						<- hivc.clu.polyphyletic.clusters(cluphy.df, cluphy.subtrees=cluphy.subtrees, plot.file=paste(outdir,'/',outfile,'_',gsub('/',':',outsignat),".pdf",sep=''), pdf.scaley=10, adj.tiplabel= c(-0.05,0.5), cex.tiplabel=0.3, pdf.xlim=0.36)
		#
		#6) plot clusters with 		long treatment interruptions			--> long TRI vs Acute ?	
		#		
		df.cluinfo				<- msm$df.cluinfo
		tmp						<- df.cluinfo[,	list(clu.bwpat.medbrl=clu.bwpat.medbrl[1],clu.npat=clu.npat[1], clu.fPossAcute=clu.fPossAcute[1], fNegT=length(which(!is.na(NegT))) / clu.ntip[1]),by="cluster"]										
		tmp						<- subset(tmp, fNegT>=quantile(tmp[,fNegT], probs=0.8) )
		tmp						<- subset(tmp, clu.npat>3)
		cluphy.df				<- merge( subset(tmp,select=cluster), df.cluinfo, all.x=1, by="cluster" )
		if(verbose) cat(paste("\nnumber of selected sequences is n=",nrow(cluphy.df)))		
		cluphy.subtrees			<- lapply( as.character(tmp[,cluster]), function(x)  msm$cluphy.subtrees[[x]]	)
		names(cluphy.subtrees)	<- as.character(tmp[,cluster])
		outdir					<- indir
		outfile					<- paste(infile,"_clust_",opt.brl,"_bs",thresh.bs*100,"_brl",thresh.brl*100,'_',"msmexpgr_selectLargeLongTRI",sep='')
		outsignat				<- insignat									
		tmp						<- hivc.clu.polyphyletic.clusters(cluphy.df, cluphy.subtrees=cluphy.subtrees, plot.file=paste(outdir,'/',outfile,'_',gsub('/',':',outsignat),".pdf",sep=''), pdf.scaley=3, adj.tiplabel= c(-0.05,0.5), cex.tiplabel=0.3, pdf.xlim=0.36)
	}
}
	
project.hivc.examl<- function(dir.name= DATA)
{
	require(ape)
	
	indir		<- paste(DATA,"derived",sep='/')
	outdir		<- paste(DATA,"tmp",sep='/')
	signat.out	<- signat.in	<- "Wed_May__1_17/08/15_2013"
	verbose		<- resume		<- 1
	infile		<- "ATHENA_2013_03_FirstAliSequences_PROTRT"
	file		<- paste(indir,'/',infile,"_",gsub('/',':',signat.in),".R",sep='')
	if(verbose) cat(paste("\nload ",file))
	load(file)
	str(seq.PROT.RT)
	
	#debug	
	
	
	cat(cmd)
	stop()
	#create ExaML binary file from phylip
	#create Parsimonator starting tree
	#run ExaML starting tree
	#delete phylip file
	
}

project.hivc.clustalo<- function(dir.name= DATA, min.seq.len=21)
{			
	if(0)
	{
		#generate fasta files per gene
		verbose<- 1
		analyze<- 1
		
		file<- paste(dir.name,"tmp/ATHENA_2013_03_Sequences.R",sep='/')
		load(file)		
		lapply(seq_along(df),function(i)
				{
					df.gene<- df[[i]]
					#print(length(df.gene[,"SampleCode"]))
					#print(length(unique(df.gene[,"SampleCode"])))
					#stop()
					df.gene[,"SampleCode"]	<- paste('>',df.gene[,"SampleCode"],sep='')
					seq.len		<- nchar(df.gene[,"Sequence"])	
					if(analyze)
					{
						cat(paste("\nprocess",names(df)[i]))
						cat(paste("\nnumber of sequences:",nrow(df.gene),"\nnumber of unique SampleCodes:", length(unique(df.gene[,"SampleCode"])),"\n") )
						hist(seq.len, breaks=100)	
						if(nrow(df.gene)!=length(unique(df.gene[,"SampleCode"])))
							o<- sapply(unique(df.gene[,"SampleCode"]),function(x)
									{
										tmp<- which(df.gene[,"SampleCode"]==x)
										if(length(tmp)>1)
										{ 
											cat("\n found duplicates\n")
											print(df.gene[tmp,])	
										}
									})
						nok.idx	<- which(seq.len<min.seq.len)
						if(verbose)	cat(paste("discarding sequences with insufficient length, n=",length(nok.idx),'\n'))
						stop("something is odd here - code acc deleted ?")
					}
					tmp			<- as.vector(t(df.gene))
					tmp			<- paste(tmp,collapse='\n')					
					file		<- paste(dir.name,"/derived/ATHENA_2013_03_Sequences_",names(df)[i],".fa",sep='')
					cat(paste("\nwrite file to",file))
					cat(tmp,file=file)	
				})			
	}
	if(0)	#edit multiple sequence alignment: remove gaps
	{
		verbose<- 1
		require(ape)
		signat.in	<- "Wed_Apr_24_09/02/03_2013"
		signat.out	<- paste(strsplit(date(),split=' ')[[1]],collapse='_',sep='')
		indir		<- paste(dir.name,"tmp",sep='/')
		pattern 	<- gsub('/',':',paste(signat.in,".clustalo$",sep=''))
		files		<- list.files(path=indir, pattern=pattern, full.names=0)
		lapply(files,function(x)
				{
					cat(paste("\nprocess",x))
					tmp		<- read.dna( paste(indir,x,sep='/'), format="fa", as.matrix=1 )
					tmp		<- hivc.seq.rmgaps(tmp, verbose=verbose)
					file	<- paste(dir.name, "tmp", gsub(pattern, paste(signat.out,"clustalo",sep='.'),x),sep='/')
					cat(paste("\nwrite to",file))
					write.dna(tmp, file, format= "fasta" )
				})		
	}	
	if(0)
	{
		#generate initial clustalo command		 
		indir	<- paste(dir.name,"derived",sep='/')
		outdir	<- paste(dir.name,"tmp",sep='/')
		infiles	<- list.files(path=indir, pattern=".fa$", full.names=0)		
		signat	<- paste(strsplit(date(),split=' ')[[1]],collapse='_',sep='')
		cmd		<- hiv.cmd.clustalo(indir, infiles, signat=signat, outdir=outdir)
		cmd		<- hiv.cmd.hpcwrapper(cmd, hpc.q="pqeph")
		
		cat(cmd[[1]])
		#lapply(cmd, cat)
		stop()
		lapply(cmd,function(x)
				{
					signat	<- paste(strsplit(date(),split=' ')[[1]],collapse='_',sep='')
					outfile	<- paste("clustalo",signat,sep='.')					
					hiv.cmd.hpccaller(outdir, outfile, x)			
				})			
	}
}
######################################################################################
hivc.prog.recombination.process.3SEQ.output<- function()
{	
	verbose		<- 1
	resume		<- 1
	indir		<- paste(DATA,"tmp",sep='/')		
	infile		<- "ATHENA_2013_03_NoDRAll+LANL_Sequences"			
	insignat	<- "Thu_Aug_01_17/05/23_2013"
	
	if(exists("argv"))
	{
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,6),
									indir= return(substr(arg,8,nchar(arg))),NA)	}))
		if(length(tmp)>0) indir<- tmp[1]		
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,7),
									infile= return(substr(arg,9,nchar(arg))),NA)	}))
		if(length(tmp)>0) infile<- tmp[1]				
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,9),
									insignat= return(substr(arg,11,nchar(arg))),NA)	}))
		if(length(tmp)>0) insignat<- tmp[1]	
		#		
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,7),
									resume= return(as.numeric(substr(arg,9,nchar(arg)))),NA)	}))
		if(length(tmp)>0) resume<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,2),
									v= return(as.numeric(substr(arg,4,nchar(arg)))),NA)	}))
		if(length(tmp)>0) verbose<- tmp[1]		
	}	
	if(verbose)
	{
		print(verbose)
		print(resume)
		print(indir)		
		print(infile)			
		print(insignat)			
	}
	file		<- paste(indir,'/',infile,'_', gsub('/',':',insignat),".R",sep='')
	if(verbose)	cat(paste("\nload file ",file))			
	load(file)
	#loaded seq.PROT.RT
	if(resume)
	{
		file		<- paste(indir,'/',infile,"_3seq_", gsub('/',':',insignat),".R",sep='')		
		options(show.error.messages = FALSE)		
		if(verbose)	cat(paste("\ntry to resume file ",file))
		readAttempt<-	try(suppressWarnings(load(file)))
		options(show.error.messages = TRUE)
		if(!inherits(readAttempt, "try-error") && verbose)	cat(paste("\nresumed file ",file))			
	}
	if(!resume || inherits(readAttempt, "try-error"))
	{
		if(verbose)	cat(paste("\ngenerate file ",file))			
		tmp			<- list.files(indir, pattern="3seq$")
		files		<- tmp[grepl(infile, tmp, fixed=1) & grepl(gsub('/',':',insignat), tmp)]
		tmp			<- sapply(strsplit(files,'_'), function(x) rev(x)[6] )		#see if filename has '_xx-xx_' string at pos 6 from end  
		files		<- files[grepl('-', tmp)]		
		if(verbose)	cat(paste("\nFound 3seq output files matching infile and insignat, n=",length(files)))
		#	figure out if any consecutive files are missing
		tmp			<- tmp[grepl('-', tmp)]
		df.seqchunks<- as.data.table( t( sapply( strsplit(tmp, '-'), as.numeric) ) )
		setnames(df.seqchunks, c("V1","V2"), c("start","end"))
		setkey(df.seqchunks, start)
		tmp			<- subset(df.seqchunks, df.seqchunks[-1,start]-1 != df.seqchunks[-nrow(df.seqchunks),end] )
		if(nrow(tmp)){		print(tmp); stop("Found missing sequence chunks")		}
		#	determine non-empty files and number of columns		
		tmp			<- sapply(files, function(x)
				{
					tmp	<- count.fields(paste(indir, '/', x, sep=''), skip=1, sep='\t')
					ifelse(length(tmp), max(tmp), 0)			
				})
		col.max		<- max(tmp)		
		files		<- files[tmp>0]
		if(verbose)	cat(paste("\nFound non-empty 3seq files, n=",length(files)))
		#	read non-empty files
		col.names		<- c("parent1","parent2","child","m","n","k","p","[p_max]","HS?","log(p)","DS(p)","DS(p)","min_rec_length")
		if(col.max>length(col.names))
			col.names	<- c(col.names, paste("bp",seq_len( col.max-length(col.names) ),sep=''))
		df.recomb		<- lapply(files, function(x)
				{	
					if(verbose)	cat(paste("\nprocess file", x))
					tmp					<- read.delim(paste(indir, '/', x, sep=''), skip=1, header=0, fill=1, col.names=col.names, stringsAsFactors=0, strip.white=T)					
					as.data.table(tmp)			
				})
		df.recomb		<- rbindlist(df.recomb)
		#	extract FASTASampleCodes from job output
		tmp				<- df.recomb[, list(m1= regexpr("-[",parent1,fixed=1), m2=regexpr("-[",parent2,fixed=1), mc=regexpr("-[",child,fixed=1)) ]
		if(any(tmp<1))	stop("Unexpected parent1 or parent2. Parsing error?")
		df.recomb		<- cbind( df.recomb, tmp )
		df.recomb[,dummy:=seq_len(nrow(df.recomb))]
		set(df.recomb, NULL, "parent1", df.recomb[,substr(parent1, 1, m1-1), by="dummy"][,V1])
		set(df.recomb, NULL, "parent2", df.recomb[,substr(parent2, 1, m2-1), by="dummy"][,V1])
		set(df.recomb, NULL, "child", df.recomb[,substr(child, 1, mc-1), by="dummy"][,V1])
		#	extract unique recombinants from job output
		setkey(df.recomb,child)
		df.recomb		<- unique(df.recomb)
		if(verbose)	cat(paste("\nFound recombinant sequences, n=",nrow(df.recomb)))
		#	order parents
		tmp				<- df.recomb[, {
					z<- sort(c(parent1, parent2))
					list(parent1=z[1], parent2=z[2])
				} ,by="dummy"]
		df.recomb		<- df.recomb[, setdiff( colnames(df.recomb), c("parent1","parent2","p","X.p_max.","HS.","DS.p.","m1","m2","mc") ), with=F]
		df.recomb		<- merge(tmp, df.recomb, by="dummy")
		#	check if triplets are unique
		tmp				<- df.recomb[, {
					z<- sort(c(parent1, parent2, child))
					list(parent1=z[1], parent2=z[2], child=z[3])
				} ,by="dummy"]
		setkeyv(tmp, c("parent1","parent2","child"))
		tmp				<- unique(tmp)
		if(nrow(tmp)<nrow(df.recomb))	warning("triplets in df.recomb not unique")
		#
		tmp				<- subset(df.recomb, select=c(parent1,parent2,child))
		setkey(tmp, parent1, parent2)
		tmp				<- unique(tmp)
		tmp[, parentpair:=seq_len(nrow(tmp))]				
		if(verbose)	cat(paste("\nNumber of unique parent pairs, n=",nrow(tmp)))
		#
		tmp				<- unique(c(df.recomb[, parent1], df.recomb[, parent2]))
		if(verbose)	cat(paste("\nNumber of unique parent sequences, n=",length(tmp)))
		#
		if(length(unique(df.recomb[,child]))!=nrow(df.recomb))	stop("Unexpected duplicate children - select parents with smallest p?")
		#
		tmp				<- subset(df.recomb, select=c(parent1,parent2))
		setkey(tmp, parent1, parent2)
		tmp				<- unique(tmp)
		tmp[, parentpair:=seq_len(nrow(tmp))]				
		if(verbose)	cat(paste("\nNumber of unique parent pairs, n=",nrow(tmp)))
		df.recomb		<- merge(df.recomb, tmp, by=c("parent1","parent2"))
		setnames(df.recomb,c("log.p.","DS.p..1"),c("logp","adjp"))
		#	candidate breakpoints bp1 etc are all overlapping, consider only bp1 for simplicity	
		#	evaluate midpoint of breapoint regions		
		tmp<- strsplit(df.recomb[,bp1], ',')
		if(any(sapply(tmp, length )!=2)) 	stop("\nunexpected bp1: missing ','")		
		df.recomb[, bp1.1:= sapply(tmp, "[", 1)]
		df.recomb[, bp1.2:= sapply(tmp, "[", 2)]
		set(df.recomb, NULL, "bp1.1", round( sapply(strsplit(df.recomb[,bp1.1], '-'), function(x)	mean(as.numeric(x))	) ) )
		set(df.recomb, NULL, "bp1.2", round( sapply(strsplit(df.recomb[,bp1.2], '-'), function(x)	min(as.numeric(x))	) ) )
		df.recomb	<- merge(df.recomb, df.recomb[, list(child.len=max(which(as.character(seq.PROT.RT[child,])!='-'))), by="dummy"], by="dummy" )
		if(any(df.recomb[, child.len-bp1.2]<0))	stop("\nunexpected breakpoint past end of child sequence")
		df.recomb[, child.start:= 1]
		tmp			<- which( df.recomb[, bp1.1<10] )
		set(df.recomb, tmp, "child.start", df.recomb[tmp, bp1.1])
		tmp			<- which( df.recomb[, child.len-bp1.2<25] )
		set(df.recomb, tmp, "child.len", as.integer(df.recomb[tmp, bp1.2]))
		#	save potential recombinants
		file		<- paste(indir,'/',infile,"_3seq_", gsub('/',':',insignat),".R",sep='')
		if(verbose)	cat(paste("\nSave candidate triplets to",file))
		save(df.recomb, file=file)		
	}
	if(0)
	{
		setnames(df.recomb, "dummy", "triplet.id")
		#
		#	get candidate recombinant sequences
		#
		df.recombseq	<- data.table( FASTASampleCode= unique( c(df.recomb[, parent1], df.recomb[, parent2], df.recomb[, child]) ) )
		df.recombseq	<- df.recombseq[ , 	{
												tmp<- subset(df.recomb, parent1==FASTASampleCode | parent2==FASTASampleCode | child==FASTASampleCode )
												list( n.triplets=nrow(tmp), min.p=min(tmp[,adjp]), med.p=median(tmp[,adjp]), triplet.id=tmp[,triplet.id] )
											} , by="FASTASampleCode"]
		setkey(df.recombseq, n.triplets, FASTASampleCode)
		#plot( df.recombseq[,n.triplets], df.recombseq[,min.p], pch=18 )
		
		#
		#determine how many other triplet sequences there are for an m candidate
		#
		df.mrecombseq	<- subset(df.recombseq, n.triplets>2)[, {
					tmp			<- triplet.id
					tmp			<- subset(df.recomb, triplet.id%in%tmp)
					triplet.seq	<- setdiff( unique( c(tmp[, parent1], tmp[, parent2], tmp[, child]) ), FASTASampleCode)
					list( triplet.seq=triplet.seq, triplet.seq.n=length(triplet.seq)  )
				}, by="FASTASampleCode"]
		df.mrecombbp	<- subset(df.recombseq, n.triplets>2)[, {
					tmp			<- triplet.id
					tmp			<- subset(df.recomb, triplet.id%in%tmp)
					list( bp1= tmp[,bp1], bp1.1= tmp[,bp1.1], bp1.2= tmp[,bp1.2] )																	
				}, by="FASTASampleCode"]														
		setkeyv(df.mrecombbp, c("FASTASampleCode", "bp1.1", "bp1.2"))												
		#	breakpoints among mrecombinants are not necessarily the same
		print(df.mrecombbp, n=250)										
	}
	df.recomb
}
######################################################################################
hivc.prog.recombination.plot.incongruence<- function()
{
	require(RColorBrewer)
	require(ape)
	
	verbose		<- 1
	resume		<- 1
	indir		<- paste(DATA,"tmp",sep='/')		
	infile		<- "ATHENA_2013_03_NoDRAll+LANL_Sequences"			
	insignat	<- "Thu_Aug_01_17/05/23_2013"
	
	id			<- NA		
	bs.n		<- 500
	select		<- ''
	#select		<- 'ng2'
	
	if(exists("argv"))
	{
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,6),
									indir= return(substr(arg,8,nchar(arg))),NA)	}))
		if(length(tmp)>0) indir<- tmp[1]		
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,7),
									infile= return(substr(arg,9,nchar(arg))),NA)	}))
		if(length(tmp)>0) infile<- tmp[1]				
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,9),
									insignat= return(substr(arg,11,nchar(arg))),NA)	}))
		if(length(tmp)>0) insignat<- tmp[1]	
		#		
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,2),
									v= return(as.numeric(substr(arg,4,nchar(arg)))),NA)	}))
		if(length(tmp)>0) verbose<- tmp[1]
		#
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,10),
									tripletid= return(as.numeric(substr(arg,12,nchar(arg)))),NA)	}))
		if(length(tmp)>0) id<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,7),
									select= return(substr(arg,9,nchar(arg))),NA)	}))
		if(length(tmp)>0) select<- tmp[1]
	}	
	if(verbose)
	{
		print(verbose)		
		print(indir)		
		print(infile)			
		print(insignat)
		print(id)		
		print(bs.n)		
		print(select)
	}
	cols				<- brewer.pal(6,"Paired")
	pattern				<- c("^in_parent1","^out_parent1","^in_parent2","^out_parent2","^in_child","^out_child")
	cex					<- 0.6
	thresh.nodesupport	<- 0.6
	edge.length.outliers<- 0.5
	
	if(is.na(id))
	{
		#	read candidate triplets
		argv				<<-	cmd.recombination.process.3SEQ.output(indir, infile, insignat, resume=1, verbose=1) 
		argv				<<- unlist(strsplit(argv,' '))
		df.recomb			<- hivc.prog.recombination.process.3SEQ.output()
		setnames(df.recomb, "dummy", "triplet.id")
		setkey(df.recomb, triplet.id)
		#	get candidate recombinant sequences
		df.recombseq	<- data.table( FASTASampleCode= unique( c(df.recomb[, parent1], df.recomb[, parent2], df.recomb[, child]) ) )
		df.recombseq	<- df.recombseq[ , 	{
					tmp<- subset(df.recomb, parent1==FASTASampleCode | parent2==FASTASampleCode | child==FASTASampleCode )
					list( n.triplets=nrow(tmp), min.p=min(tmp[,adjp]), med.p=median(tmp[,adjp]), triplet.id=tmp[,triplet.id] )
				} , by="FASTASampleCode"]
		setkey(df.recombseq, FASTASampleCode)
		#	select unclear triplets	
		if(select=="ng2")
		{
			tmp				<- unique(subset(df.recombseq, n.triplets>2)[, triplet.id])			
			df.recomb		<- df.recomb[J(setdiff(df.recomb[,triplet.id], tmp )),]
			if(verbose)	cat(paste("\nSelected unclear triplets, n=", nrow(df.recomb)))
		}
		#	select triplets that have a given FASTASampleCode
		if(select %in% unique(subset(df.recombseq, n.triplets>2)[, FASTASampleCode]))
		{
			tmp				<- subset(df.recombseq, n.triplets>2)[J(select),]
			df.recomb		<- df.recomb[J(tmp[,triplet.id]),]
		}
		#	analyze FASTASampleCodes that occur in >2 triplets
		if(select=="g2")
		{
			tmp				<- unique(subset(df.recombseq, n.triplets>2)[, FASTASampleCode])
			dummy			<- lapply(tmp, function(x)
						{
							argv			<<- cmd.recombination.plot.incongruence(indir, infile, insignat, prog= PR.RECOMB.PLOTINCONGRUENCE, opt.select=x,verbose=1)
							argv			<<- unlist(strsplit(argv,' '))
							hivc.prog.recombination.plot.incongruence()			
						})
			stop()
		}
		
		#	read available checks for triplets
		files				<- list.files(indir, pattern=".newick$")
		files				<- files[ grepl(paste('_3seqcheck_',sep=''), files) & grepl(infile, files, fixed=1) & grepl(gsub('/',':',insignat), files) ]
		tmp					<- regexpr("3seqcheck_",files)
		tmp					<- sapply(seq_along(files), function(i){		substr(files[i],tmp[i],nchar(files[i]))		})
		files.df			<- data.table(	file=files, 
											triplet.id= sapply(strsplit(tmp, '_'), function(x)		substr(x[[2]],3,nchar(x[[2]]))	),
											region= sapply(strsplit(tmp, '_'), function(x)		substr(x[[3]],2,nchar(x[[3]]))	)		)
		set(files.df, NULL, "triplet.id", as.numeric(files.df[,triplet.id]))							
		setkey(files.df, triplet.id, region)
		if(verbose)	cat(paste("\nFound files, n=", nrow(files.df)))
		files.df			<- merge(files.df, files.df[,list(region.n= length(region)), by="triplet.id"], by="triplet.id")
		files.df			<- subset(files.df, region.n>1)
		if(verbose)	cat(paste("\nFound checked triplets, n=", nrow(files.df)/2))
		if(verbose) cat(paste("\nTriplets still to check=",paste(setdiff( df.recomb[,triplet.id],unique(files.df[,triplet.id]) ), collapse=', ')))
		#	select candidate triplets for which checks available
		files.df			<- merge(files.df, df.recomb,by="triplet.id")
		if(verbose)	cat(paste("\nSelected triplets for plotting, n=", nrow(files.df)/2))
		setkey(files.df, parentpair, adjp)
		#setkey(files.df, adjp)
		
		file				<- paste(indir,'/',infile,"_3seqcheck_examlbs",bs.n,'_',select,'_',gsub('/',':',insignat),".pdf",sep='')
		if(verbose)	cat(paste("\nPlot both phylogenies to", file))		
		pdf(file, width=12, height=6)
		def.par 			<- par(no.readonly = TRUE)
		par(mar=c(2,0.5,1,0))
		layout(matrix(c(1,1,2,2), 2, 2))					
		dummy				<- lapply( unique( files.df[, triplet.id] ), function(z)
				{
					#x<- subset(files.df, triplet.id==58)
					x<- subset(files.df, triplet.id==z)
					#	read tree corresponding to 'in' recombinant region
					tmp					<- paste(indir,'/',x[1, file],sep='')
					if(verbose)	cat(paste("\nRead file", tmp))
					ph.in				<- ladderize( read.tree(tmp) )		
					tmp					<- as.numeric( ph.in$node.label )
					tmp[is.na(tmp)]		<- 0 
					ph.in$node.label	<- tmp/100
					#	read tree corresponding to 'out' recombinant region
					tmp					<- paste(indir,'/',x[2, file],sep='')		
					if(verbose)	cat(paste("\nRead file", tmp))
					ph.out				<- ladderize( read.tree(tmp) )		
					tmp					<- as.numeric( ph.out$node.label )
					tmp[is.na(tmp)]		<- 0 
					ph.out$node.label	<- tmp/100
					#	remove outlier filler sequences
					outliers			<- c()
					tmp					<- which( ph.in$edge.length>edge.length.outliers )
					if(length(tmp))
					{
						tmp				<- ph.in$edge[tmp,2]
						outliers		<- ph.in$tip.label[ tmp[tmp<Ntip(ph.in)] ]
					}
					if(length(tmp))
					{						
						tmp				<- which( ph.out$edge.length>1 )
						tmp				<- ph.out$edge[tmp,2]
						outliers		<- c(outliers, ph.out$tip.label[ tmp[tmp<Ntip(ph.out)] ])
					}
					if(length(outliers))
					{
						ph.in			<- drop.tip(ph.in, outliers)
						ph.out			<- drop.tip(ph.out, outliers)
					}
					#	show half way stable subtrees -- ideally would contain triplet sequences
					clustering.in		<- hivc.clu.clusterbythresh(ph.in, thresh.nodesupport=thresh.nodesupport, nodesupport=ph.in$node.label, retval="all")
					clustering.out		<- hivc.clu.clusterbythresh(ph.out, thresh.nodesupport=thresh.nodesupport, nodesupport=ph.out$node.label, retval="all")
					#	show filler sequences in different color
					tip.color.in		<- rep("black", Ntip(ph.in))	
					for(i in seq_along(pattern))
					{
						tmp						<- which( grepl(pattern[i],ph.in$tip) )
						if(length(tmp))	
							tip.color.in[ tmp ]	<- cols[i]
					}
					tip.color.out		<- rep("black", Ntip(ph.out))	
					for(i in seq_along(pattern))
					{
						tmp						<- which( grepl(pattern[i],ph.out$tip) )
						if(length(tmp))	
							tip.color.out[ tmp ]<- cols[i]
					}	
					#	plot both phylogenies side by side	
					dummy<- hivc.clu.plot(ph.in, clustering.in[["clu.mem"]], show.tip.label=T, cex.nodelabel=cex, cex.edge.incluster=2*cex, no.margin=F, tip.color=tip.color.in)
					mtext(paste("parent.pair=",x[,parentpair]," triplet.id=",x[1,triplet.id]," log10p=",round(log10(x[1,adjp]),d=2)," region 'in' length=",x[1,bp1.2-bp1.1]), side = 3, cex=cex)
					axisPhylo(cex=cex)
					dummy<- hivc.clu.plot(ph.out, clustering.out[["clu.mem"]], show.tip.label=T, cex.nodelabel=cex, cex.edge.incluster=2*cex, no.margin=F, tip.color=tip.color.out)		
					mtext(paste("parent.pair=",x[,parentpair]," triplet.id=",x[1,triplet.id]," log10p=",round(log10(x[1,adjp]),d=2)," region 'out' length=",x[1,child.len-child.start-bp1.2+bp1.1]), side = 3, cex=cex)
					axisPhylo(cex=cex)
				})
		par(def.par)
		dev.off()
	
		#x<- "R12-15108"
		#subset(df.recomb, parent1==x | parent2==x | child==x )
	}
	if(!is.na(id))
	{
		#	read tree corresponding to 'in' recombinant region
		file				<- paste(indir,'/',infile,"_3seqcheck_id",id,"_rIn_examlbs",bs.n,'_',gsub('/',':',insignat),".newick",sep='')
		if(verbose)	cat(paste("\nRead file", file))
		ph.in				<- ladderize( read.tree(file) )		
		tmp					<- as.numeric( ph.in$node.label )
		tmp[is.na(tmp)]		<- 0 
		ph.in$node.label	<- tmp/100
		#	read tree corresponding to 'out' recombinant region
		file				<- paste(indir,'/',infile,"_3seqcheck_id",id,"_rOut_examlbs",bs.n,'_',gsub('/',':',insignat),".newick",sep='')
		if(verbose)	cat(paste("\nRead file", file))
		ph.out				<- ladderize( read.tree(file) )		
		tmp					<- as.numeric( ph.out$node.label )
		tmp[is.na(tmp)]		<- 0 
		ph.out$node.label	<- tmp/100
		#	remove outlier filler sequences
		outliers			<- c()
		tmp					<- which( ph.in$edge.length>edge.length.outliers )
		if(length(tmp))
		{
			tmp				<- ph.in$edge[tmp,2]
			outliers		<- ph.in$tip.label[ tmp[tmp<Ntip(ph.in)] ]
		}
		if(length(tmp))
		{						
			tmp				<- which( ph.out$edge.length>1 )
			tmp				<- ph.out$edge[tmp,2]
			outliers		<- c(outliers, ph.out$tip.label[ tmp[tmp<Ntip(ph.out)] ])
		}
		if(length(outliers))
		{
			ph.in			<- drop.tip(ph.in, outliers)
			ph.out			<- drop.tip(ph.out, outliers)
		}
		#	show half way stable subtrees -- ideally would contain triplet sequences
		clustering.in		<- hivc.clu.clusterbythresh(ph.in, thresh.nodesupport=thresh.nodesupport, nodesupport=ph.in$node.label, retval="all")
		clustering.out		<- hivc.clu.clusterbythresh(ph.out, thresh.nodesupport=thresh.nodesupport, nodesupport=ph.out$node.label, retval="all")
		#	show filler sequences in different color
		tip.color.in		<- rep("black", Ntip(ph.in))	
		for(i in seq_along(pattern))
		{
			tmp						<- which( grepl(pattern[i],ph.in$tip) )
			if(length(tmp))	
				tip.color.in[ tmp ]	<- cols[i]
		}
		tip.color.out		<- rep("black", Ntip(ph.out))	
		for(i in seq_along(pattern))
		{
			tmp						<- which( grepl(pattern[i],ph.out$tip) )
			if(length(tmp))	
				tip.color.out[ tmp ]<- cols[i]
		}	
		#	plot both phylogenies side by side
		file				<- paste(indir,'/',infile,"_3seqcheck_id",id,"_examlbs",bs.n,'_',gsub('/',':',insignat),".pdf",sep='')
		if(verbose)	cat(paste("\nPlot both phylogenies to", file))
		pdf(file, width=8, height=6)
		def.par 			<- par(no.readonly = TRUE)
		layout(matrix(c(1,1,2,2), 2, 2))		
		dummy<- hivc.clu.plot(ph.in, clustering.in[["clu.mem"]], show.tip.label=T, cex.nodelabel=cex, cex.edge.incluster=2*cex, no.margin=F, tip.color=tip.color.in)
		mtext(paste("parent.pair=",x[,parentpair]," triplet.id=",x[1,triplet.id]," log10p=",round(log10(x[1,adjp]),d=2)," region 'in' length=",x[1,bp1.2-bp1.1]), side = 3, cex=cex)
		axisPhylo(cex=cex)
		dummy<- hivc.clu.plot(ph.out, clustering.out[["clu.mem"]], show.tip.label=T, cex.nodelabel=cex, cex.edge.incluster=2*cex, no.margin=F, tip.color=tip.color.out)		
		mtext(paste("parent.pair=",x[,parentpair]," triplet.id=",x[1,triplet.id]," log10p=",round(log10(x[1,adjp]),d=2)," region 'out' length=",x[1,child.len-child.start-bp1.2+bp1.1]), side = 3, cex=cex)
		axisPhylo(cex=cex)
		par(def.par)
		dev.off()
	}
}
######################################################################################
hivc.prog.recombination.check.candidates<- function()
{	
	require(ape)
	verbose		<- 1
	resume		<- 0
	indir		<- paste(DATA,"tmp",sep='/')		
	infile		<- "ATHENA_2013_03_NoDRAll+LANL_Sequences"			
	insignat	<- "Thu_Aug_01_17/05/23_2013"
	
	id			<- 51
	seq.select.n<- 10
	bs.from		<- 0
	bs.to		<- 499
	bs.n		<- 500
	
	hpc.walltime<- 36
	hpc.mem		<- "600mb"
	hpc.nproc	<- 1		
	hpc.q		<- "pqeph"
	
	if(exists("argv"))
	{
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,6),
									indir= return(substr(arg,8,nchar(arg))),NA)	}))
		if(length(tmp)>0) indir<- tmp[1]		
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,7),
									infile= return(substr(arg,9,nchar(arg))),NA)	}))
		if(length(tmp)>0) infile<- tmp[1]				
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,9),
									insignat= return(substr(arg,11,nchar(arg))),NA)	}))
		if(length(tmp)>0) insignat<- tmp[1]	
		#		
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,7),
									resume= return(as.numeric(substr(arg,9,nchar(arg)))),NA)	}))
		if(length(tmp)>0) resume<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,2),
									v= return(as.numeric(substr(arg,4,nchar(arg)))),NA)	}))
		if(length(tmp)>0) verbose<- tmp[1]
		#
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,10),
									tripletid= return(as.numeric(substr(arg,12,nchar(arg)))),NA)	}))
		if(length(tmp)>0) id<- tmp[1]
	}	
	if(verbose)
	{
		print(verbose)
		print(resume)
		print(indir)		
		print(infile)			
		print(insignat)
		print(id)
		print(seq.select.n)
		print(bs.from)
		print(bs.to)
		print(bs.n)		
	}
	if(resume)
	{
		tmp			<- 1
		file		<- paste(indir,'/',infile,"_3seqcheck_id",id,"_rIn_",gsub('/',':',insignat),".R",sep='')
		options(show.error.messages = FALSE)		
		if(verbose)	cat(paste("\ntry to load file ",file))
		readAttempt	<-	try(suppressWarnings(load(file)))
		options(show.error.messages = TRUE)
		if(!inherits(readAttempt, "try-error") && verbose)		cat(paste("\nloaded file=",file))
		if(inherits(readAttempt, "try-error"))					tmp		<- 0
		
		file		<- paste(indir,'/',infile,"_3seqcheck_id",id,"_rOut_",gsub('/',':',insignat),".R",sep='')
		options(show.error.messages = FALSE)		
		if(verbose)	cat(paste("\ntry to load file ",file))
		readAttempt	<-	try(suppressWarnings(load(file)))
		options(show.error.messages = TRUE)
		if(!inherits(readAttempt, "try-error") && verbose)		cat(paste("\nloaded file=",file))
		if(inherits(readAttempt, "try-error"))					tmp		<- 0		
	}
	if(!resume || tmp==0)
	{
		file		<- paste(indir,'/',infile,'_', gsub('/',':',insignat),".R",sep='')
		if(verbose)	cat(paste("\nload file ",file))			
		load(file)
		#	loaded seq.PROT.RT
		file		<- paste(indir,'/',infile,"_3seq_", gsub('/',':',insignat),".R",sep='')		
		options(show.error.messages = FALSE)		
		if(verbose)	cat(paste("\ntry to load file ",file))
		readAttempt	<-	try(suppressWarnings(load(file)))
		options(show.error.messages = TRUE)
		if(inherits(readAttempt, "try-error"))	stop(paste("\nCannot find 3SEQ file, run hivc.prog.recombination.process.3SEQ.output?, file=",file))			
		#	loaded df.recomb
		
		#
		#	process triplet for dummy id	
		#	
		df.recomb	<- subset(df.recomb, dummy==id)
		if(verbose)	cat(paste("\nprocess triplet number",id))
		if(verbose)	print(df.recomb)
		#	create sequence matrices corresponding to the two breakpoint regions 
		seq.in		<- seq.PROT.RT[,seq.int(df.recomb[,bp1.1],df.recomb[,bp1.2])]
		seq.out		<- if(df.recomb[,child.start]<df.recomb[,bp1.1]-1) seq.int(df.recomb[,child.start],df.recomb[,bp1.1]-1) else numeric(0) 
		seq.out		<- if(df.recomb[,bp1.2]+1<df.recomb[,child.len]) c(seq.out,seq.int(df.recomb[,bp1.2]+1, df.recomb[,child.len]))	else 	seq.out
		seq.out		<- seq.PROT.RT[,seq.out]
		seq.select.f<- ifelse(min(ncol(seq.out),ncol(seq.in))<150, 10, 10)
		if(verbose)	cat(paste("\nsetting inflation factor to",seq.select.f))
		seq.select.n<- seq.select.n * seq.select.f
		#	select background sequences for child based on sequence similarity
		if(verbose)	cat(paste("\ncompute genetic distances for parent1 parent2 child"))
		tmp				<- which( rownames(seq.PROT.RT)==df.recomb[,child] )		
		dummy			<- 0				
		seq.dist					<- 1 - sapply(seq_len(nrow(seq.in))[-tmp],function(i){		.C("hivc_dist_ambiguous_dna", seq.in[tmp,], seq.in[i,], ncol(seq.in), dummy )[[4]]			})	
		seq.dist[is.nan(seq.dist)]	<- Inf
		seq.df			<- data.table( FASTASampleCode=rownames(seq.in)[-tmp] , dist=seq.dist, group="child", region="in" ) 		
		seq.dist					<- 1 - sapply(seq_len(nrow(seq.out))[-tmp],function(i){		.C("hivc_dist_ambiguous_dna", seq.out[tmp,], seq.out[i,], ncol(seq.out), dummy )[[4]]			})	
		seq.dist[is.nan(seq.dist)]	<- Inf
		seq.df			<- rbind(seq.df, data.table( FASTASampleCode=rownames(seq.out)[-tmp] , dist=seq.dist, group="child", region="out" )) 
		#	select background sequences for parent1 based on sequence similarity
		tmp				<- which( rownames(seq.PROT.RT)==df.recomb[,parent1] )				
		seq.dist					<- 1 - sapply(seq_len(nrow(seq.in))[-tmp],function(i){		.C("hivc_dist_ambiguous_dna", seq.in[tmp,], seq.in[i,], ncol(seq.in), dummy )[[4]]			})	
		seq.dist[is.nan(seq.dist)]	<- Inf
		seq.df			<- rbind(seq.df,data.table( FASTASampleCode=rownames(seq.in)[-tmp] , dist=seq.dist, group="parent1", region="in" )) 		
		seq.dist					<- 1 - sapply(seq_len(nrow(seq.out))[-tmp],function(i){		.C("hivc_dist_ambiguous_dna", seq.out[tmp,], seq.out[i,], ncol(seq.out), dummy )[[4]]			})	
		seq.dist[is.nan(seq.dist)]	<- Inf
		seq.df			<- rbind(seq.df, data.table( FASTASampleCode=rownames(seq.out)[-tmp] , dist=seq.dist, group="parent1", region="out" )) 
		#	select background sequences for parent2 based on sequence similarity
		tmp				<- which( rownames(seq.PROT.RT)==df.recomb[,parent2] )				
		seq.dist					<- 1 - sapply(seq_len(nrow(seq.in))[-tmp],function(i){		.C("hivc_dist_ambiguous_dna", seq.in[tmp,], seq.in[i,], ncol(seq.in), dummy )[[4]]			})	
		seq.dist[is.nan(seq.dist)]	<- Inf
		seq.df			<- rbind(seq.df,data.table( FASTASampleCode=rownames(seq.in)[-tmp] , dist=seq.dist, group="parent2", region="in" )) 		
		seq.dist					<- 1 - sapply(seq_len(nrow(seq.out))[-tmp],function(i){		.C("hivc_dist_ambiguous_dna", seq.out[tmp,], seq.out[i,], ncol(seq.out), dummy )[[4]]			})	
		seq.dist[is.nan(seq.dist)]	<- Inf
		seq.df			<- rbind(seq.df, data.table( FASTASampleCode=rownames(seq.out)[-tmp] , dist=seq.dist, group="parent2", region="out" ))
		#	first pass:
		#	select closest FASTASampleCode by group and region to verify recombination breakpoint by phylogenetic incongruence
		setkeyv(seq.df, c("group","region","dist"))
		seq.df			<- subset(seq.df, dist>0)
		if(verbose)	cat(paste("\nFound related sequences with dist>0, n=",nrow(seq.df)))
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
		tmp				<- c( df.recomb[,parent1],df.recomb[,parent2],df.recomb[,child], seq.df[, FASTASampleCode] )
		tmp				<- hivc.seq.unique(seq.in[tmp,])
		seq.df			<- merge( data.table(FASTASampleCode=rownames(tmp)), seq.df, by="FASTASampleCode" )
		tmp				<- c( df.recomb[,parent1],df.recomb[,parent2],df.recomb[,child], seq.df[, FASTASampleCode] )
		tmp				<- hivc.seq.unique(seq.out[tmp,])
		seq.df			<- merge( data.table(FASTASampleCode=rownames(tmp)), seq.df, by="FASTASampleCode" )
		seq.df			<- subset(seq.df, !FASTASampleCode%in%c( df.recomb[,parent1],df.recomb[,parent2],df.recomb[,child] ) )		
		if(verbose)	cat(paste("\nfound candidates for filler sequences that are unique on both recombinant regions, n=",nrow(seq.df)))
		if(verbose)	print( seq.df[	,	list(n=length(FASTASampleCode)) ,by=c("group","region")] )
		#		
		seq.select.n	<- seq.select.n/seq.select.f 
		#	select 'seq.select.n' unique filler sequences for 'in' region, balancing by group as much as possible
		if(verbose)	cat(paste("\nSelect filler sequences for recombinant region 'in'"))
		tmp						<- subset( seq.df, region=='in' )[,FASTASampleCode]
		tmp						<- hivc.seq.unique(seq.in[tmp ,])
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
		if(overflow>0)	cat(paste("\nNot as many filler sequences as requested for recombinant region 'in', n=",nrow(seq.out.df)))
		seq.in.df				<- ans[-nrow(ans),]
		if(verbose)	cat(paste("\nSelected balancing sequences for recombinant region 'in', n=",nrow(seq.in.df)))
		if(verbose)	print( seq.in.df[	,	list(n=length(FASTASampleCode)) ,by=c("group","region")] )
		#
		#	select 'seq.select.n' unique filler sequences for 'out' region, balancing by group as much as possible
		#
		if(verbose)	cat(paste("\nSelect sequences for recombinant region 'out'"))
		tmp						<- subset( seq.df, region=='out' )[,FASTASampleCode] 
		tmp						<- hivc.seq.unique(seq.out[tmp,])		
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
		if( nrow(hivc.seq.unique( seq.in[ seq.df[, FASTASampleCode], ] ))!=nrow(seq.df) ) stop("Unexpected non-unique 'in' sequences")
		if( nrow(hivc.seq.unique( seq.out[ seq.df[, FASTASampleCode], ] ))!=nrow(seq.df) ) stop("Unexpected non-unique 'out' sequences")
		#	could be that the triplet sequences are not unique among each other 
		#	if so, make change to one of the triplet sequences
		seq.select				<- c( df.recomb[,parent1],df.recomb[,parent2],df.recomb[,child], seq.df[, FASTASampleCode] )
		tmp						<- hivc.seq.unique( seq.in[ seq.select, ] )
		if( !length(setdiff(c( df.recomb[,parent1],df.recomb[,parent2],df.recomb[,child] ),rownames(tmp))) )
			seq.in				<- tmp
		if( length(setdiff(c( df.recomb[,parent1],df.recomb[,parent2],df.recomb[,child] ),rownames(tmp))) )
		{
			tmp					<- setdiff(c( df.recomb[,parent1],df.recomb[,parent2],df.recomb[,child] ),rownames(tmp) )		#name of sequence in triplet that is identical with one other sequence in triplet
			if(verbose)	cat(paste("\nFound identical triplet sequence for region 'in'", tmp))
			seq.in				<- as.character(seq.in)
			seq.in[tmp,1]		<- ifelse(seq.in[tmp,1]=='t','c',ifelse(seq.in[tmp,1]=='c','t',ifelse(seq.in[tmp,1]=='a','g','a')))
			seq.in				<- as.DNAbin(seq.in)
			tmp					<- hivc.seq.unique( seq.in[ seq.select, ] )
			if( length(setdiff(c( df.recomb[,parent1],df.recomb[,parent2],df.recomb[,child] ),rownames(tmp))) )	stop("Unexpected duplicate for 'in'")
			seq.in				<- tmp
		}			
		seq.out					<- seq.out[seq.select,]
		tmp						<- hivc.seq.unique(seq.out)
		if( !length(setdiff(c( df.recomb[,parent1],df.recomb[,parent2],df.recomb[,child] ),rownames(tmp))) )
			seq.out				<- tmp
		if( length(setdiff(c( df.recomb[,parent1],df.recomb[,parent2],df.recomb[,child] ),rownames(tmp))) )
		{
			tmp					<- setdiff(c( df.recomb[,parent1],df.recomb[,parent2],df.recomb[,child] ),rownames(tmp) )		#name of sequence in triplet that is identical with one other sequence in triplet
			if(verbose)	cat(paste("\nFound identical triplet sequence for region 'out'", tmp))
			seq.out				<- as.character(seq.out)
			seq.out[tmp,1]		<- ifelse(seq.out[tmp,1]=='t','c',ifelse(seq.out[tmp,1]=='c','t',ifelse(seq.out[tmp,1]=='a','g','a')))
			seq.out				<- as.DNAbin(seq.out)
			tmp					<- hivc.seq.unique(seq.out[seq.select,])
			if( length(setdiff(c( df.recomb[,parent1],df.recomb[,parent2],df.recomb[,child] ),rownames(tmp))) )	stop("Unexpected duplicate for 'out'")
			seq.out				<- tmp
		}
		if(any(rownames(seq.in)!=rownames(seq.out)))	stop("Unexpected unequal sequences selected")
		#	reset rownames
		tmp						<- c( paste(c("tparent1","tparent2","tchild"),seq.select[1:3],sep='_'), seq.df[,list(label= paste(region,group,FASTASampleCode,sep='_')), by="FASTASampleCode"][,label] )
		rownames(seq.in)		<- tmp
		rownames(seq.out)		<- tmp
		#	save
		file		<- paste(indir,'/',infile,"_3seqcheck_id",id,"_rIn_",gsub('/',':',insignat),".R",sep='')
		if(verbose) cat(paste("\nsave to ",file))
		save(seq.in, file=file)		
		file		<- paste(indir,'/',infile,"_3seqcheck_id",id,"_rOut_",gsub('/',':',insignat),".R",sep='')
		if(verbose) cat(paste("\nsave to ",file))
		save(seq.out, file=file)
	}
	if(1)
	{
		#
		#	run bootstrap ExaML for region 'in', all boostraps on one processor
		#			
		cmd				<- NULL
		file			<- paste(indir,'/',infile,"_3seqcheck_id",id,"_rIn_examlbs",bs.n,'_',gsub('/',':',insignat),".newick",sep='')		
		if(!resume || !file.exists(file))
		{
			infile.exa	<- paste(infile,"_3seqcheck_id",id,"_rIn",sep='')		
			cmd			<- cmd.examl.bootstrap.on.one.machine(indir, infile.exa, gsub('/',':',insignat),gsub('/',':',insignat), bs.from=bs.from, bs.to=bs.to, bs.n=bs.n, outdir=indir, opt.bootstrap.by="nucleotide", resume=1, verbose=1)
		}
		#
		#	run bootstrap ExaML for region 'out', all boostraps on one processor
		#						
		file			<- paste(indir,'/',infile,"_3seqcheck_id",id,"_rOut_examlbs",bs.n,'_',gsub('/',':',insignat),".newick",sep='')
		if(!resume || !file.exists(file))
		{
			infile.exa	<- paste(infile,"_3seqcheck_id",id,"_rOut",sep='')		
			cmd			<- c(cmd, cmd.examl.bootstrap.on.one.machine(indir, infile.exa, gsub('/',':',insignat),gsub('/',':',insignat), bs.from=bs.from, bs.to=bs.to, bs.n=bs.n, outdir=indir, opt.bootstrap.by="nucleotide", resume=1, verbose=1))
		}
		#
		if(verbose) cat(paste("\ncreated ExaML bootstrap runs, n=",length(cmd)))
		if(!is.null(cmd))
		{
			cmd			<- paste(cmd,collapse='\n')
			#cmd		<- paste(cmd,cmd.recombination.plot.incongruence(indir, infile, gsub('/',':',insignat), triplet.id=id, verbose=1),sep='')				
			#cat(cmd)
			if(verbose) cat(paste("\nqsub ExaML bootstrap runs, hpc.walltime=",hpc.walltime," hpc.mem=",hpc.mem," hpc.nproc=",hpc.nproc," hpc.q=",hpc.q))
			cmd			<- cmd.hpcwrapper(cmd, hpc.walltime=hpc.walltime, hpc.q=hpc.q, hpc.mem=hpc.mem, hpc.nproc=hpc.nproc)
			signat		<- paste(strsplit(date(),split=' ')[[1]],collapse='_',sep='')
			outdir		<- paste(DATA,"tmp",sep='/')
			outfile		<- paste("3sc",signat,sep='.')
			#cat(cmd)			
			cmd.hpccaller(outdir, outfile, cmd)
			Sys.sleep(1)
		}
	}
}		
######################################################################################
hivc.prog.eval.clustering.bias<- function()
{
	indir				<- paste(DATA,"tmp",sep='/')		
	infile				<- "ATHENA_2013_03_NoDRAll+LANL_Sequences"		
	insignat			<- "Thu_Aug_01_17/05/23_2013"
	indircov			<- paste(DATA,"derived",sep='/')
	infilecov			<- "ATHENA_2013_03_AllSeqPatientCovariates"
	infiletree			<- paste(infile,"examlbs100",sep="_")
	
	opt.brl				<- "dist.brl.casc" 
	thresh.brl			<- 0.096
	thresh.bs			<- 0.8
	resume				<- 1
	verbose				<- 1		
	#
	#	load msm clusters
	#
	argv			<<- hivc.cmd.clustering.msm(indir, infiletree, insignat, indircov, infilecov, opt.brl, thresh.brl, thresh.bs, resume=resume)
	argv			<<- unlist(strsplit(argv,' '))		
	msm				<- hivc.prog.get.clustering.MSM()	
	#
	df.all			<- msm$df.cluinfo
	file.out.name	<- paste( infiletree,"_clust_",opt.brl,"_bs",thresh.bs*100,"_brl",thresh.brl*100,"_msmexpgr", sep='')
	#
	#	new diagnoses by CD4
	#
	df.newdiag			<- copy(subset(df.all, select=c(Patient,AnyPos_T1, CD4_T1)))
	setkey(df.newdiag,Patient)
	df.newdiag			<- unique(df.newdiag)
	msmclu.newdiagCD4 	<- hivc.db.getplot.newdiagnosesbyCD4(	df.newdiag, 
			plot.file= paste(dir.name,"/tmp/",file.out.name,"_NewDiagByCD4_",gsub('/',':',insignat),".pdf",sep=''),
			plot.file.p= paste(dir.name,"/tmp/",file.out.name,"_NewDiagByCD4_prop_",gsub('/',':',insignat),".pdf",sep=''),
			plot.ylab= "New diagnoses HIV-1 subtype B,\nin MSM cluster")
	#
	#	seem in care by risk group
	#
	df.living			<- copy(subset(df.all, select=c(Patient, AnyPos_T1, Trm, DateDied, DateLastContact)))
	setkey(df.living,Patient)
	df.living			<- unique(df.living)		
	msmclu.livExpGr 	<- hivc.db.getplot.livingbyexposure(	df.living, 
			plot.file=paste(dir.name,"/tmp/",file.out.name,"_Seen4CareByExpGroup_",gsub('/',':',insignat),".pdf",sep=''),
			plot.file.p=paste(dir.name,"/tmp/",file.out.name,"_Seen4CareByExpGroup_prop_",gsub('/',':',insignat),".pdf",sep=''),
			plot.ylab="Seen for care with HIV-1 subtype B,\nin MSM cluster", db.endtime=2013.3, db.diff.lastcontact2died=0.5, db.diff.lastcontact2now= 2.3, verbose=1)
	#
	#	seen in care by recent, untreated/CD4, treated
	#
	file.immu			<- paste(dir.name,"derived/ATHENA_2013_03_Immu.R",sep='/')
	load(file.immu)
	df.immu				<- df
	df.immu				<- subset(df.immu, select=c(Patient,PosCD4, CD4) )
	set(df.immu, NULL, "PosCD4", hivc.db.Date2numeric(df.immu[,PosCD4]))
	
	df.living			<- copy(subset(df.all, (is.na(df.all[,Trm]) | Trm%in%c("MSM","BI","HET")) & Sex!='F', select=c(Patient, AnyPos_T1, CD4_T1, isAcute, AnyT_T1, DateDied, DateLastContact)))
	setkey(df.living,Patient)
	df.living			<- unique(df.living)
	msmclu.livCD4		<- hivc.db.getplot.livingbyCD4(	df.living, df.immu, 
			plot.file=paste(dir.name,"/tmp/",file.out.name,"_MSMSeen4CareByCD4_",gsub('/',':',insignat),".pdf",sep=''),
			plot.file.p=paste(dir.name,"/tmp/",file.out.name,"_MSMSeen4CareByCD4_prop_",gsub('/',':',insignat),".pdf",sep=''), 
			plot.ylab="MSM seen for care with HIV-1 subtype B,\nin MSM cluster", db.endtime=2013.3, db.diff.lastcontact2died=0.5, db.diff.lastcontact2now= 2.3, verbose=1)
	#
	#
	#	load ATHENA_03 data, subset to MSM and BI
	#
	#
	file			<- paste(dir.name,"/derived/",infilecov,".R",sep='')
	load(file)		
	#
	#	seem in care by risk group
	#		
	df.newdiag		<- copy(subset(df.all, (is.na(df.all[,Trm]) | Trm%in%c("MSM","BI")) & Sex!='F', select=c(Patient,AnyPos_T1, CD4_T1)))
	setkey(df.newdiag,Patient)
	df.newdiag		<- unique(df.newdiag)
	msm.newdiagCD4 	<- hivc.db.getplot.newdiagnosesbyCD4(	df.newdiag, 
															plot.file= paste(dir.name,"/derived/",infilecov,"_MSMNewDiagByCD4.pdf",sep=''),
															plot.file.p= paste(dir.name,"/derived/",infilecov,"_MSMNewDiagByCD4_prop.pdf",sep=''),
															plot.ylab= "MSM new diagnoses HIV-1 subtype B,\n seq available")
	#
	df.living		<- copy(subset(df.all, (is.na(df.all[,Trm]) | Trm%in%c("MSM","BI")) & Sex!='F', select=c(Patient, AnyPos_T1, Trm, DateDied, DateLastContact)))
	setkey(df.living,Patient)
	df.living		<- unique(df.living)
	msm.livExpGr	<- hivc.db.getplot.livingbyexposure(	df.living, 
															plot.file=paste(dir.name,"/derived/",infilecov,"_MSMSeen4CareByExpGroup.pdf",sep=''),
															plot.file.p=paste(dir.name,"/derived/",infilecov,"_MSMSeen4CareByExpGroup_prop.pdf",sep=''),
															plot.ylab="MSM seen for care with HIV-1 subtype B,\n seq available", db.endtime=2013.3, db.diff.lastcontact2died=0.5, db.diff.lastcontact2now= 2.3, verbose=1)
	#
	#	seen in care by recent, untreated/CD4, treated
	#
	file.immu		<- paste(dir.name,"derived/ATHENA_2013_03_Immu.R",sep='/')
	load(file.immu)
	df.immu			<- df
	df.immu			<- subset(df.immu, select=c(Patient,PosCD4, CD4) )
	set(df.immu, NULL, "PosCD4", hivc.db.Date2numeric(df.immu[,PosCD4]))
	
	df.living		<- copy(subset(df.all, (is.na(df.all[,Trm]) | Trm%in%c("MSM","BI")) & Sex!='F', select=c(Patient, AnyPos_T1, CD4_T1, isAcute, AnyT_T1, DateDied, DateLastContact)))
	setkey(df.living,Patient)
	df.living		<- unique(df.living)	
	msm.livCD4		<- hivc.db.getplot.livingbyCD4(			df.living, df.immu, 
															plot.file=paste(dir.name,"/derived/",infilecov,"_MSMSeen4CareByCD4.pdf",sep=''),
															plot.file.p=paste(dir.name,"/derived/",infilecov,"_MSMSeen4CareByCD4_prop.pdf",sep=''),
															plot.ylab="MSM seen for care with HIV-1 subtype B,\n seq available", db.endtime=2013.3, db.diff.lastcontact2died=0.5, db.diff.lastcontact2now= 2.3, verbose=1)						
	#	new diagnoses by year	
	clu		<- msmclu.newdiagCD4$t.newdiag
	ds		<- msm.newdiagCD4$t.newdiag
	tmp		<- apply(clu, 2, sum)[as.character(1985:2012)] / apply(ds, 2, sum)[as.character(1985:2012)]
	tmp		<- round(tmp, d=2)
	#
	clu		<- msmclu.newdiagCD4$p.newdiag
	ds		<- msm.newdiagCD4$p.newdiag	
	tmp		<- round( clu[, as.character(1999:2012)] - ds[, as.character(1999:2012)], d=2)
	apply(tmp,1,mean)
	#	seen in care
	clu		<- msmclu.livCD4$p.living
	ds		<- msm.livCD4$p.living
	tmp		<- round( clu[, as.character(1996:2012)] - ds[, as.character(1996:2012)], d=2)
	
	#
}			
######################################################################################
hivc.prog.get.clustering<- function()
{
	library(ape)
	library(data.table)	
			
	opt.dist.brl	<- "dist.brl.casc"
	thresh.bs		<- 0.8
	thresh.brl		<- 0.096	#with mut rate 6e-3/yr expect 2*0.048 brl for 8 yrs apart	
	verbose			<- 1
	resume			<- 0
	indir			<- paste(DATA,"tmp",sep='/')
	infile			<- "ATHENA_2013_03_CurAll+LANL_Sequences_examlbs100"
	insignat		<- "Sat_Jun_16_17/23/46_2013"

	if(exists("argv"))
	{
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,6),
									indir= return(substr(arg,8,nchar(arg))),NA)	}))
		if(length(tmp)>0) indir<- tmp[1]		
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,7),
									infile= return(substr(arg,9,nchar(arg))),NA)	}))
		if(length(tmp)>0) infile<- tmp[1]				
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,9),
									insignat= return(substr(arg,11,nchar(arg))),NA)	}))
		if(length(tmp)>0) insignat<- tmp[1]		
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,10),
									thresh.bs= return(as.numeric(substr(arg,12,nchar(arg)))),NA)	}))
		if(length(tmp)>0) thresh.bs<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,11),
									thresh.brl= return(as.numeric(substr(arg,13,nchar(arg)))),NA)	}))
		if(length(tmp)>0) thresh.brl<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,8),
									opt.brl= return(substr(arg,10,nchar(arg))),NA)	}))
		if(length(tmp)>0) opt.dist.brl<- tmp[1]		
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,7),
									resume= return(as.numeric(substr(arg,9,nchar(arg)))),NA)	}))
		if(length(tmp)>0) resume<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,2),
									v= return(as.numeric(substr(arg,4,nchar(arg)))),NA)	}))
		if(length(tmp)>0) verbose<- tmp[1]		
	}	
	
	if(verbose)
	{
		cat(paste("\nopt.dist.brl",opt.dist.brl))
		cat(paste("\nthresh.bs",thresh.bs))
		cat(paste("\nthresh.brl",thresh.brl))
		cat(paste("\nindir",indir))
		cat(paste("\ninfile",infile))
		cat(paste("\ninsignat",insignat))				
	}
	outdir			<- indir
	outfile			<- paste(infile,"_clust_",opt.dist.brl,"_bs",thresh.bs*100,"_brl",thresh.brl*100,sep='')
	outsignat		<- insignat	
	
	if(resume)
	{
		file		<- paste(outdir,'/',outfile,'_',gsub('/',':',outsignat),".R",sep='')		
		options(show.error.messages = FALSE)		
		if(verbose)
			cat(paste("\ntry to resume file ",file))
		readAttempt<-	try(suppressWarnings(load(file)))
		options(show.error.messages = TRUE)
		if(!inherits(readAttempt, "try-error") && verbose)
			cat(paste("\nresumed file ",file))
	}
	if(!resume || inherits(readAttempt, "try-error"))
	{		
		#
		#	load preclustlink files
		#
		file		<- paste(indir,paste(infile,"_preclust_",gsub('/',':',insignat),".R",sep=''),sep='/')
		if(verbose) cat(paste("load file",file))
		options(show.error.messages = FALSE)
		readAttempt<-	try(suppressWarnings(load(file)))
		if(inherits(readAttempt, "try-error"))	stop(paste("cannot load required file", file))
		options(show.error.messages = TRUE)
		#
		#	generate clustering
		#
		dist.brl		<- switch(	opt.dist.brl, 
									"dist.brl.max"		= dist.brl.max,
									"dist.brl.med"		= dist.brl.med,
									"dist.brl.casc"		= dist.brl.casc,
									NA)
		if(any(is.na(dist.brl)))	stop("unexpected NA in dist.brl")					
		clustering		<- hivc.clu.clusterbythresh(ph, thresh.nodesupport=thresh.bs, thresh.brl=thresh.brl, dist.brl=dist.brl, nodesupport=ph.node.bs,retval="all")
		#
		#	add cluster membership to df.seqinfo
		#
		setkey(df.seqinfo, Node)
		df.seqinfo[,"cluster":= clustering[["clu.mem"]][seq_len(Ntip(ph))]]
		#	set CountryInfection for non-ATHENA seqs
		tmp				<- which( df.seqinfo[,substr(FASTASampleCode,1,2)=="TN"] )
		set(df.seqinfo, tmp, "CountryInfection", "FRGNTN")
		tmp				<- which( df.seqinfo[,substr(FASTASampleCode,1,8)=="PROT+P51"] )
		set(df.seqinfo, tmp, "CountryInfection", "FRGN")
		#
		#	compute TP clusters
		#
		clusters.tp			<- hivc.clu.truepos(clustering, ph.linked, Ntip(ph), verbose= 0)								
		missed.df			<- copy( df.seqinfo )
		plotfile			<- paste(outdir,paste(outfile,"_missedtp_",gsub('/',':',outsignat),".pdf",sep=''),sep='/')
		hivc.clu.plot.withinpatientseq.not.samecluster(clu.pre$ph, clustering, clusters.tp, missed.df, plotfile)
		#
		#	compute TN clusters
		#		
		clusters.tn		<- hivc.clu.trueneg(clustering, ph.unlinked.info, ph.unlinked, Ntip(ph), verbose=0)
		#
		#	save
		#		
		file			<- paste(outdir,'/',outfile,'_',gsub('/',':',outsignat),".R",sep='')
		cat(paste("\nwrite cluster info to file",file))
		save(df.seqinfo, clustering, clusters.tp, clusters.tn, file=file)
	}
	
	ans	<- list(thresh.bs=thresh.bs, thresh.brl=thresh.brl, opt.dist.brl=opt.dist.brl, df.seqinfo=df.seqinfo, clustering=clustering, clusters.tp=clusters.tp, clusters.tn=clusters.tn)
	ans
}
######################################################################################
hivc.prog.get.clustering.precompute<- function()
{
	library(ape)
	#library(adephylo)
	library(data.table)	
	
	verbose				<- 1
	resume				<- 0
	use.seroneg.as.is	<- 1	#use with updated "ATHENA_2013_03_AllSeqPatientCovariates"
	indir				<- paste(DATA,"tmp",sep='/')	
	infile				<- "ATHENA_2013_03_CurAll+LANL_Sequences_examlbs100"
	insignat			<- "Sat_Jun_16_17/23/46_2013"
	indircov			<- paste(DATA,"derived",sep='/')	
	infilecov			<- "ATHENA_2013_03_AllSeqPatientCovariates"	
	
	if(exists("argv"))
	{
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,6),
									indir= return(substr(arg,8,nchar(arg))),NA)	}))
		if(length(tmp)>0) indir<- tmp[1]		
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,7),
									infile= return(substr(arg,9,nchar(arg))),NA)	}))
		if(length(tmp)>0) infile<- tmp[1]				
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,9),
									insignat= return(substr(arg,11,nchar(arg))),NA)	}))
		if(length(tmp)>0) insignat<- tmp[1]		
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,9),
									indircov= return(substr(arg,11,nchar(arg))),NA)	}))
		if(length(tmp)>0) indircov<- tmp[1]		
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,10),
									infilecov= return(substr(arg,12,nchar(arg))),NA)	}))
		if(length(tmp)>0) infilecov<- tmp[1]		
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,7),
									resume= return(as.numeric(substr(arg,9,nchar(arg)))),NA)	}))
		if(length(tmp)>0) resume<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,2),
									v= return(as.numeric(substr(arg,4,nchar(arg)))),NA)	}))
		if(length(tmp)>0) verbose<- tmp[1]		
	}	
	outdir			<- indir
	outsignat		<- insignat
	outfile			<- paste(infile,"preclust",sep='_')
	outfile.dtips	<- paste(infile,"predtips",sep='_')
	if(verbose)
	{
		print(indir)
		print(infile)
		print(insignat)
		print(indircov)
		print(infilecov)
		print(resume)
	}
	if(resume)
	{
		file		<- paste(outdir,paste(outfile,'_',gsub('/',':',outsignat),".R",sep=''),sep='/')
		options(show.error.messages = FALSE)		
		if(verbose)
			cat(paste("\ntry to resume file ",file))
		readAttempt<-	try(suppressWarnings(load(file)))
		options(show.error.messages = TRUE)
		if(!inherits(readAttempt, "try-error") && verbose)
			cat(paste("\nresumed file ",file))
	}
	if(!resume || inherits(readAttempt, "try-error"))
	{	
		file							<- paste(indir,paste(infile,'_',gsub('/',':',insignat),".newick",sep=''),sep='/')
		if(verbose) cat(paste("\nread phylo from file",file))		
		ph								<- ladderize( read.tree(file) )
		gc()
		#
		#	easy: extract bootstrap support
		#
		ph$node.label[2]				<- 0								#little hack so that clustering works
		ph.node.bs						<- as.numeric( ph$node.label )
		ph.node.bs[is.na(ph.node.bs)]	<- 0
		ph.node.bs						<- ph.node.bs/100
		ph$node.label					<- ph.node.bs
		if(1)
		{
			#
			#	memory consuming: extract branch length statistic of subtree
			#
			if(verbose) cat("\ncompute dist.brl.med")
			dist.brl.med					<- hivc.clu.brdist.stats(ph, eval.dist.btw="leaf", stat.fun=median)
			gc()
			if(verbose) cat("\ncompute dist.brl.max")
			dist.brl.max					<- hivc.clu.brdist.stats(ph, eval.dist.btw="leaf", stat.fun=max)		#read patristic distances -- this is the expensive step but still does not take very long
			gc()
			if(verbose) cat("\ncompute dist.brl.casc")
			dist.brl.casc					<- hivc.clu.brdist.stats(ph, eval.dist.btw="leaf", stat.fun=hivc.clu.min.transmission.cascade)
			gc()		
		}
		#
		#	easy: extract tree specific TP and FN data sets
		#		
		file							<- paste(indircov,"/",infilecov,".R",sep='')
		load(file)
		if(verbose) cat(paste("\nfound covariates for patients, n=",nrow(df.all)))
		df.seqinfo						<- subset(df.all, !is.na(PosSeqT) )
		if(verbose) cat(paste("\nfound covariates for patients with non-NA PosSeqT, n=",nrow(df.seqinfo)))
		if(verbose) cat(paste("\nstart: compute TP and TN data tables for phylogeny"))
		tmp								<- hivc.phy.get.TP.and.TN(ph, df.seqinfo, verbose=verbose, use.seroneg.as.is= use.seroneg.as.is)		
		if(verbose) cat(paste("\nend: compute TP and TN data tables for phylogeny"))
		unlinked.byspace				<- tmp[["unlinked.byspace"]]
		unlinked.bytime					<- tmp[["unlinked.bytime"]]
		linked.bypatient				<- tmp[["linked.bypatient"]]	
		ph.linked						<- tmp[["ph.linked"]]
		ph.unlinked.info				<- tmp[["ph.unlinked.info"]]
		ph.unlinked						<- tmp[["ph.unlinked"]]
		#
		#	add node number to df.seqinfo
		#
		df.seqinfo						<- merge( df.all, data.table(Node=seq_len(Ntip(ph)), FASTASampleCode=ph$tip.label), all.y=1, by="FASTASampleCode")
		#
		#	get and plot distribution of bootstrap values for TP and TN pairs
		#
		ph.mrca				<- mrca(ph)
		#	prepare data.tables with mrca: dead/seroneg TN pairs		
		bs.unlinkedpairs	<- lapply(unlinked.bytime, function(x)
								{						
									set(x, NULL, "query.FASTASampleCode", as.character(x[,query.FASTASampleCode]))
									set(x, NULL, "FASTASampleCode", as.character(x[,FASTASampleCode]))
									ans					<- merge(x, x[, list(mrca= ph.mrca[query.FASTASampleCode,FASTASampleCode]) ,by=FASTASampleCode],by="FASTASampleCode")
									setnames(ans, c("FASTASampleCode","query.FASTASampleCode"), c("tip2","tip1"))
									subset(ans, select=c(tip1, tip2, mrca))
								})
		bs.unlinkedpairs	<- rbindlist(bs.unlinkedpairs)
		#	prepare data.tables with mrca: geographically distant seqs and any indb seq
		unlinked.byspace[,dummy:=seq_len(nrow(unlinked.byspace))]
		set(unlinked.byspace, NULL, "FASTASampleCode", as.character(unlinked.byspace[,FASTASampleCode]))	
		seq.indb			<- colnames(ph.mrca)[ which( substr(colnames(ph.mrca),1,2)!="TN" & substr(colnames(ph.mrca),1,8)!="PROT+P51" ) ]
		bs.unlinked.byspace	<- unlinked.byspace[,	list(tip1=FASTASampleCode,tip2=seq.indb, mrca= ph.mrca[seq.indb,FASTASampleCode]), by="dummy"]
		#	prepare data.tables with mrca: within patient TP pairs
		setkey(linked.bypatient,Patient)
		bs.linked.bypatient	<- linked.bypatient[, {
														tmp					<- match(FASTASampleCode, ph$tip.label)
														ans					<- t(combn(tmp, 2 ))
														ans					<- cbind(ans, apply(ans, 1, function(z)  ph.mrca[z[1],z[2]]))
														data.table(tip1=ans[,1], tip2=ans[,2], mrca=ans[,3])
													}, by="Patient"]		
		tmp								<- hivc.phy.get.TP.and.TN.bootstrapvalues(ph, bs.linked.bypatient, ph.mrca=ph.mrca, df.seqinfo=df.seqinfo, bs.unlinkedpairs=bs.unlinkedpairs, bs.unlinked.byspace=bs.unlinked.byspace, dist.brl=dist.brl.casc, thresh.brl=0.096, plot.file=paste(outdir,'/',outfile,"_distbs_",gsub('/',':',outsignat),".pdf",sep=''), verbose=verbose)		
		bs.linked.bypatient				<- tmp[["bs.linked.bypatient"]] 
		bs.unlinkedpairs				<- tmp[["bs.unlinkedpairs"]]
		bs.unlinked.byspace				<- tmp[["bs.unlinked.byspace"]]
		#
		#	memory consuming: compute branch length matrix between tips
		#
		#if(verbose) cat("\ncompute dist.root")
		#dist.root						<-  distRoot(ph, method= "patristic")
		#gc()
		#if(verbose) cat("\ncompute dist.tips.mat")
		#dist.tips.mat					<-  distTips(ph, method= "patristic")
		#gc()		
		#
		#	save output
		#
		file							<- paste(outdir,paste(outfile,'_',gsub('/',':',outsignat),".R",sep=''),sep='/')	
		if(verbose)	cat(paste("write to file",file))
		save(	ph, dist.brl.max, dist.brl.med, dist.brl.casc, ph.node.bs, ph.linked, ph.unlinked.info, ph.unlinked, 
				df.seqinfo, unlinked.byspace, unlinked.bytime, linked.bypatient,  
				bs.linked.bypatient, bs.unlinkedpairs, bs.unlinked.byspace, file=file )
		#save(ph, dist.brl.max, dist.brl.med, dist.brl.casc, ph.node.bs, ph.linked, ph.unlinked.info, ph.unlinked, df.seqinfo, unlinked.byspace, unlinked.bytime, linked.bypatient,  file=file )		
	}
	
	ans	<- list(	ph=ph, #dist.tips.mat=dist.tips.mat, #dist.root=dist.root,  
					dist.brl.max=dist.brl.max, dist.brl.med=dist.brl.med, dist.brl.casc=dist.brl.casc, 
					ph.node.bs=ph.node.bs, ph.linked=ph.linked, ph.unlinked.info=ph.unlinked.info, ph.unlinked=ph.unlinked, 
					df.seqinfo=df.seqinfo, unlinked.byspace=unlinked.byspace, unlinked.bytime=unlinked.bytime, linked.bypatient=linked.bypatient,
					bs.linked.bypatient=bs.linked.bypatient, bs.unlinkedpairs=bs.unlinkedpairs, bs.unlinked.byspace
					)
	ans				
}

hivc.prog.remove.resistancemut<- function()
{
	library(ape)
	library(data.table)
	
	#load drug resistance mutations and select unique mutants by codon
	dr		<- as.data.table( read.csv( paste( CODE.HOME,"/data/IAS_primarydrugresistance_201303.csv",sep='' ), stringsAsFactors=F ) )	
	dr[,Alignment.nuc.pos:= (Gene.codon.number-1)*3+Gene.HXB2pos ]		
	dr		<- dr[,	{	tmp<- unique(Mutant); list(Mutant=tmp, Gene.codon.number=Gene.codon.number[1], Wild.type=Wild.type[1], DR.name=DR.name[1])	}, by=Alignment.nuc.pos]
	#select nucleotide codes that are consistent with drug resistance mutants
	nt2aa	<- as.data.table( read.csv( paste( CODE.HOME,"/data/standard_nt_code.csv",sep='' ), stringsAsFactors=F ) )
	setnames(nt2aa,c("AA","NTs"),c("Mutant","Mutant.NTs"))
	nt2aa	<- subset(nt2aa, select=c(Mutant,Mutant.NTs))
	dr		<- merge(dr, nt2aa, all.x=1, by="Mutant", allow.cartesian=TRUE)
	setkey(dr, "Alignment.nuc.pos")
	#print(dr, nrows=250)
	dr		<- subset(dr, select=c(Alignment.nuc.pos, Mutant.NTs, DR.name))
	set(dr, NULL, "Mutant.NTs", tolower(dr[,Mutant.NTs]))

	indir			<- paste(DATA,"tmp",sep='/')
	infile			<- "ATHENA_2013_03_CurAll+LANL_Sequences"
	insignat		<- "Sat_Jun_16_17/23/46_2013"
	outdir			<- paste(DATA,"tmp",sep='/')
	outfile			<- "ATHENA_2013_03_NoDRAll+LANL_Sequences"
	outsignat		<- "Thu_Aug_01_17/05/23_2013"
	alignment.start	<- 2253	
	verbose			<- 1
	if(exists("argv"))
	{
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,6),
									indir= return(substr(arg,8,nchar(arg))),NA)	}))
		if(length(tmp)>0) indir<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,7),
									infile= return(substr(arg,9,nchar(arg))),NA)	}))
		if(length(tmp)>0) infile<- tmp[1]		
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,9),
									insignat= return(substr(arg,11,nchar(arg))),NA)	}))
		if(length(tmp)>0) insignat<- tmp[1]		
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,7),
									outdir= return(substr(arg,9,nchar(arg))),NA)	}))
		if(length(tmp)>0) outdir<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,8),
									outfile= return(substr(arg,10,nchar(arg))),NA)	}))
		if(length(tmp)>0) outfile<- tmp[1]		
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,10),
									outsignat= return(substr(arg,12,nchar(arg))),NA)	}))
		if(length(tmp)>0) outsignat<- tmp[1]	
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,16),
									alignment.start= return(substr(arg,18,nchar(arg))),NA)	}))
		if(length(tmp)>0) alignment.start<- tmp[1]	
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,2),
									v= return(as.numeric(substr(arg,4,nchar(arg)))),NA)	}))
		if(length(tmp)>0) verbose<- tmp[1]
	}
	if(verbose)
	{
		print(indir)
		print(infile)
		print(insignat)
		print(outdir)
		print(outfile)		
		print(outsignat)
		print(alignment.start)
		print(verbose)
	}	
	#load alignment
	file				<- paste(indir,'/',infile,'_',gsub('/',':',insignat),".R",sep='')
	load(file)	
	#modify dr table for particular alignment	
	set(dr, NULL, "Alignment.nuc.pos", dr[,Alignment.nuc.pos]-alignment.start+1)
	
	#remove	likely.nonB.outliers
	cat("\nchange infile: remove	likely.nonB.outliers")
	likely.nonB.outliers	<- c("R03-07193","2006G206","PROT+P51_B.AU.1995.C92.AF538307","2008G084")
	likely.nonB.outliers	<- which(rownames(seq.PROT.RT) %in% likely.nonB.outliers)
	seq.PROT.RT				<- seq.PROT.RT[-likely.nonB.outliers,]
	
	#if alignment equals any of the drug resistance mutants, replace with NNN	
	seq.PROT.RT			<- hivc.seq.rm.drugresistance(as.character(seq.PROT.RT), dr, verbose=verbose, rtn.DNAbin=1 )	
	
	#save No Drug resistance alignment to file
	file								<- paste(outdir,'/',outfile,'_',gsub('/',':',outsignat),".R",sep='')
	if(verbose)	cat(paste("\nwrite R file to",file))
	save(seq.PROT.RT, file=file)	
	file								<- paste(outdir,'/',outfile,'_',gsub('/',':',outsignat),".phylip",sep='')
	if(verbose)	cat(paste("\nwrite phylip file to",file))
	hivc.seq.write.dna.phylip(seq.PROT.RT, file=file)					
	file								<- paste(outdir,'/',outfile,'_',gsub('/',':',outsignat),".fasta",sep='')			
	if(verbose)	cat(paste("\nwrite fasta file to",file))
	write.dna(seq.PROT.RT, file=file, format="fasta", colsep='', colw=ncol(seq.PROT.RT), blocksep=0)
	
	seq.PROT.RT
}
######################################################################################
hivc.prog.BEAST.evalpoolrun<- function()
{	
	require(phangorn)
	indircov			<- paste(DATA,"derived",sep='/')	
	infilecov			<- "ATHENA_2013_03_AllSeqPatientCovariates"			
	file.cov			<- paste(indircov,"/",infilecov,".R",sep='')
	file.viro			<- paste(indircov,"/ATHENA_2013_03_Viro.R",sep='/')
	file.immu			<- paste(indircov,"/ATHENA_2013_03_Immu.R",sep='/')
	file.treatment		<- paste(indircov,"/ATHENA_2013_03_Regimens.R",sep='/')
		
	indir				<- paste(DATA,"beast/beast_131011",sep='/')
	#indir				<- paste(DATA,"tmp",sep='/')	
	infile				<- "ATHENA_2013_03_NoDRAll+LANL_Sequences_beast_seroneg"
	insignat			<- "Tue_Aug_26_09/13/47_2013"
	infilexml.opt		<- "txs4clu"
	infilexml.opt		<- "mph4clu"
	infilexml.opt		<- "mph4clutx4tip"
	infilexml.template	<- "um232rhU2045"
	infilexml.template	<- "um182rhU2045"
	#infilexml.opt		<- "mph4cluLdTd"
	#infilexml.template	<- "um22rhG202018"
	#infilexml.template	<- "um182rhU2045ay"	
	
	plot						<- 1
	verbose						<- 1
	resume						<- 0
	pool.n						<- 3	
	beastlabel.idx.clu			<- 1
	beastlabel.idx.hivs			<- ifelse(any(sapply(c("LdTd","LsTd"), function(x) grepl(x,infilexml.opt))),	5,	4)
	beastlabel.idx.samplecode	<- 6	
	beastlabel.idx.rate			<- 7
	
	if(exists("argv"))
	{
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,6),
									indir= return(substr(arg,8,nchar(arg))),NA)	}))
		if(length(tmp)>0) indir<- tmp[1]		
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,7),
									infile= return(substr(arg,9,nchar(arg))),NA)	}))
		if(length(tmp)>0) infile<- tmp[1]				
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,9),
									insignat= return(substr(arg,11,nchar(arg))),NA)	}))
		if(length(tmp)>0) insignat<- tmp[1]		
		#		
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,2),
									v= return(as.numeric(substr(arg,4,nchar(arg)))),NA)	}))
		if(length(tmp)>0) verbose<- tmp[1]
		#
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,14),
									infilexml.opt= return(substr(arg,16,nchar(arg))),NA)	}))
		if(length(tmp)>0) infilexml.opt<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,19),
									infilexml.template= return(substr(arg,21,nchar(arg))),NA)	}))
		if(length(tmp)>0) infilexml.template<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,7),
									pool.n= return(as.numeric(substr(arg,9,nchar(arg)))),NA)	}))
		if(length(tmp)>0) pool.n<- tmp[1]
	}	
	if(verbose)
	{
		print(indir)
		print(infile)
		print(insignat)
		print(infilexml.opt)
		print(infilexml.template)
		print(pool.n)
	}
	if(resume)
	{
		file		<- paste(indir,'/',infile,'_',infilexml.template,'_',infilexml.opt,"_mcc_",gsub('/',':',insignat),".R",sep='')		
		options(show.error.messages = FALSE)		
		if(verbose)
			cat(paste("\ntry to resume file ",file))
		readAttempt<-	try(suppressWarnings(load(file)))
		options(show.error.messages = TRUE)
		if(!inherits(readAttempt, "try-error") && verbose)
			cat(paste("\nresumed file ",file))
		return(list(cluphy=cluphy, cluphy.trmca=cluphy.trmca, cluphy.df=cluphy.df))
	}
	if(!resume || inherits(readAttempt, "try-error"))
	{
		#	
		#	read annotated mcc trees
		#
		tmp			<- list.files(indir, pattern=paste(".nex$",sep=''))
		file.nex	<- tmp[ grepl(paste('_',infilexml.template,'_',sep=''), tmp) & grepl(paste('_',infilexml.opt,'_',sep=''), tmp) & grepl(gsub('/',':',insignat), tmp) ]
		if(verbose)	cat(paste("\nFound .nex files matching input args, n=", length(file.nex)))
		tmp			<- list.files(indir, pattern=paste(".log$",sep=''))
		file.log	<- tmp[ grepl(paste('_',infilexml.template,'_',sep=''), tmp) & grepl(paste('_',infilexml.opt,'_',sep=''), tmp) & grepl(gsub('/',':',insignat), tmp) ]
		if(verbose)	cat(paste("\nFound .log files matching input args, n=", length(file.log)))		
		#		
		if(length(file.nex)==pool.n)
		{
			#	read treeannotator .nex file
			ph.beast		<- lapply(file.nex, function(x)		hivc.treeannotator.read(paste(indir,x,sep='/'), add.to.tiplabel=c("rate_median"), rate.multiplier=1e3, round.digit=2, verbose=verbose)		)
			if(verbose)	cat(paste("\nRead trees matching input args, n=", length(ph.beast)))
			#	read length of tip stems
			file.log		<- sapply(file.log, function(x)					paste(indir,x,sep='/') )
			file.xml		<- sapply(file.log, function(x)					paste(substr(x,1,nchar(x)-3), "xml", sep='') )			
			df.tstem		<- lapply(seq_along(file.log), function(i)		hivc.beast.read.log2tstem(file.log[i], file.xml[i], beastlabel.idx.samplecode=6, burn.in= 5e6, breaks.n= 30, verbose=1)		)
			df.tstem		<- rbindlist( df.tstem )
			#	load all patient covariates
			load(file.cov)								
			#	build single tree from pooled runs
			tmp				<- hivc.treeannotator.get.phy(ph.beast, beastlabel.idx.clu=beastlabel.idx.clu, beastlabel.idx.hivs=beastlabel.idx.hivs, beastlabel.idx.samplecode=beastlabel.idx.samplecode, beastlabel.idx.rate=beastlabel.idx.rate, debug=0)
			cluphy			<- tmp$cluphy 
			ph.tip.ctime	<- tmp$ph.tip.ctime 				
			ph.root.ctime	<- tmp$ph.root.ctime
			#	extract rates		
			rates.df		<- hivc.treeannotator.get.rates(cluphy, tmp$ph.tip.df, nodelabel.idx.edgewidth=5)
			#	convert tstem time into calendar time
			tmp				<- data.table( FASTASampleCode=cluphy$tip.label, tip=seq_along(cluphy$tip.label), mrca= Ancestors(cluphy, seq_along(cluphy$tip.label), type="parent")-Ntip(cluphy) )			 
			df.tstem		<- merge( df.tstem, tmp, by="FASTASampleCode" )			
			set(df.tstem, NULL, "tstem", df.tstem[, max(ph.tip.ctime)-tstem])
			df.tstem		<- subset(df.tstem, select=c(tip, mrca, tstem, density))
			#	get length of 95% TMRCAs of tip stems
			cluphy.tstem	<- df.tstem[,	{
												x<- data.table(tstem, density)
												x[,dummy:=-density]
												setkey(x,dummy)			
												x[,cdensity:= cumsum(x[,density]/sum(x[,density]))]							
												list(height_95_diff= diff(range(subset(x, cdensity<0.95, tstem))) )							
											},	by="tip"]			
			#
			cluphy.df		<- hivc.treeannotator.get.clusterprob(ph.beast, beastlabel.idx.clu=beastlabel.idx.clu, beastlabel.idx.samplecode=beastlabel.idx.samplecode)
			cluphy.df		<- merge(cluphy.df, df.all, by="FASTASampleCode")
			cluphy.tmrca	<- hivc.treeannotator.get.tmrcas(ph.beast, beastlabel.idx.hivs=beastlabel.idx.hivs) 						
			#
			#	save mcc tree in R format
			#
			file		<- paste(indir,'/',infile,'_',infilexml.template,'_',infilexml.opt,"_mcc_",gsub('/',':',insignat),".R",sep='')
			if(verbose)	cat(paste("\nSave mcc tree objects to file",file))
			save(cluphy, cluphy.tmrca, cluphy.df, cluphy.tstem, file=file)			
			#
			#	plot cluster TMRCAs
			#
			if(plot)
			{
				tmp			<- sort(cluphy.tmrca[, height_95_diff])	
				file		<- paste(indir,'/',infile,'_',infilexml.template,'_',infilexml.opt,"_mcc_clutmrca_",gsub('/',':',insignat),".pdf",sep='')
				if(verbose)	cat(paste("\nPlot cluster TMRCAs to file",file))
				pdf(file,width=7,height=7)
				plot( tmp, seq_len(nrow(cluphy.tmrca)), type="l", xlab='width of 95% interval of TMRCA of clusters', ylab='#MRCA of clusters', bty='n' )
				tmp			<- quantile(tmp, prob=c(0.25,0.5,0.8,0.9,0.95))
				legend("bottomright",bty='n',border=NA,legend=paste( apply(rbind(names(tmp), round(tmp,d=2)),2,function(x) paste(x,collapse='=',sep='')), collapse=', ',sep=''))
				legend("topleft",bty='n',border=NA,legend=paste(infilexml.template,infilexml.opt,sep='_'))
				dev.off()
			}
			#
			#	plot tip stem TMRCAs
			#
			if(plot)
			{
				tmp			<- sort(cluphy.tstem[, height_95_diff])	
				file		<- paste(indir,'/',infile,'_',infilexml.template,'_',infilexml.opt,"_mcc_tiptmrca_",gsub('/',':',insignat),".pdf",sep='')
				if(verbose)	cat(paste("\nPlot tip TMRCAs to file",file))
				pdf(file,width=7,height=7)
				plot( tmp, seq_len(nrow(cluphy.tstem)), type="l", xlab='width of 95% interval of TMRCA of tips', ylab='#MRCA of tips', bty='n' )
				tmp			<- quantile(tmp, prob=c(0.25,0.5,0.8,0.9,0.95))
				legend("bottomright",bty='n',border=NA,legend=paste( apply(rbind(names(tmp), round(tmp,d=2)),2,function(x) paste(x,collapse='=',sep='')), collapse=', ',sep=''))
				legend("topleft",bty='n',border=NA,legend=paste(infilexml.template,infilexml.opt,sep='_'))
				dev.off()
			}
			#
			#	plot clusters
			#
			if(plot)
			{
				#load patient RNA
				load(file.viro)
				df.viro				<- df				
				#load patient CD4				
				load(file.immu)
				df.immu				<- df
				#load patient regimen
				load(file.treatment)
				df.treatment		<- df
				
				youngest.tip.ctime	<- max(ph.tip.ctime)
				#youngest.tip.ctime	<- 2010.46
				file				<- paste(indir,'/',infile,'_',infilexml.template,'_',infilexml.opt,"_mcc_",gsub('/',':',insignat),".pdf",sep='')
				if(verbose)	cat(paste("\nplotting dated clusters to file", file ))
				dummy				<- hivc.treeannotator.plot(cluphy, ph.root.ctime, youngest.tip.ctime, df.all, df.viro, df.immu, df.treatment=df.treatment, df.tstem=df.tstem, df.rates=rates.df, end.ctime=2013.3, cex.nodelabel=0.5, cex.tiplabel=0.5, file=file, pdf.width=7, pdf.height=150)
			}
		}
		else if(verbose)	
			cat(paste("\nNothing evaluated - waiting for BEAST runs to complete"))
	}
}
######################################################################################
hivc.prog.BEAST2.generate.xml<- function()
{	
	require(XML)
	
	indir				<- paste(DATA,"tmp",sep='/')		
	infile				<- "ATHENA_2013_03_NoDRAll+LANL_Sequences"		
	insignat			<- "Thu_Aug_01_17/05/23_2013"
	indircov			<- paste(DATA,"derived",sep='/')
	infilecov			<- "ATHENA_2013_03_AllSeqPatientCovariates"
	infiletree			<- paste(infile,"examlbs100",sep="_")
	infilexml			<- paste(infile,'_',"seroneg",sep='')	
	outdir				<- indir
	outsignat			<- "Tue_Aug_26_09/13/47_2013"
	
	opt.brl					<- "dist.brl.casc" 
	thresh.brl				<- 0.096
	thresh.bs				<- 0.8
	pool.ntip				<- 130
	beast.mcmc.length		<- 25e6
	infilexml.opt			<- "standard"
	infilexml.template		<- "bdsky_hky"	
	
	resume				<- 1
	verbose				<- 1
	hpc.walltime		<- 71
	hpc.ncpu			<- 1
	hpc.mem				<- "1200mb"
	
	if(exists("argv"))
	{
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,6),
									indir= return(substr(arg,8,nchar(arg))),NA)	}))
		if(length(tmp)>0) indir<- tmp[1]		
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,7),
									infile= return(substr(arg,9,nchar(arg))),NA)	}))
		if(length(tmp)>0) infile<- tmp[1]				
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,9),
									insignat= return(substr(arg,11,nchar(arg))),NA)	}))
		if(length(tmp)>0) insignat<- tmp[1]		
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,9),
									indircov= return(substr(arg,11,nchar(arg))),NA)	}))
		if(length(tmp)>0) indircov<- tmp[1]		
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,10),
									infilecov= return(substr(arg,12,nchar(arg))),NA)	}))
		if(length(tmp)>0) infilecov<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,11),
									infiletree= return(substr(arg,13,nchar(arg))),NA)	}))
		if(length(tmp)>0) infiletree<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,10),
									infilexml= return(substr(arg,12,nchar(arg))),NA)	}))
		if(length(tmp)>0) infilexml<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,7),
									outdir= return(substr(arg,9,nchar(arg))),NA)	}))
		if(length(tmp)>0) outdir<- tmp[1]		
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,10),
									outsignat= return(substr(arg,12,nchar(arg))),NA)	}))
		if(length(tmp)>0) outsignat<- tmp[1]				
		#		
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,7),
									resume= return(as.numeric(substr(arg,9,nchar(arg)))),NA)	}))
		if(length(tmp)>0) resume<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,2),
									v= return(as.numeric(substr(arg,4,nchar(arg)))),NA)	}))
		if(length(tmp)>0) verbose<- tmp[1]
		#
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,8),
									opt.brl= return(substr(arg,10,nchar(arg))),NA)	}))
		if(length(tmp)>0) opt.brl<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,10),
									thresh.bs= return(as.numeric(substr(arg,12,nchar(arg)))),NA)	}))
		if(length(tmp)>0) thresh.bs<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,11),
									thresh.brl= return(as.numeric(substr(arg,13,nchar(arg)))),NA)	}))
		if(length(tmp)>0) thresh.brl<- tmp[1]
		#
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,10),
									pool.ntip= return(as.numeric(substr(arg,12,nchar(arg)))),NA)	}))
		if(length(tmp)>0) pool.ntip<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,14),
									infilexml.opt= return(substr(arg,16,nchar(arg))),NA)	}))
		if(length(tmp)>0) infilexml.opt<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,19),
									infilexml.template= return(substr(arg,21,nchar(arg))),NA)	}))
		if(length(tmp)>0) infilexml.template<- tmp[1]
	}	
	
	#	modify beast2.spec depending on infilexml.opt
	if(grepl("S4p",infilexml.opt))
	{
		beast2.spec		<- hivc.beast2.get.specifications(mcmc.length=beast.mcmc.length, bdsky.intervalNumber=4)
		beast2.spec$bdsky.sprop.changepoint.value	<- beast2.spec$bdsky.R0.changepoint.value		<- beast2.spec$bdsky.notInf.changepoint.value	<- c(9.596, 5.596, 1.596, 0.)
	}
	else if(grepl("S5p",infilexml.opt))
	{
		beast2.spec		<- hivc.beast2.get.specifications(mcmc.length=beast.mcmc.length, bdsky.intervalNumber=5)
		beast2.spec$bdsky.sprop.changepoint.value	<- beast2.spec$bdsky.R0.changepoint.value		<- beast2.spec$bdsky.notInf.changepoint.value	<- c(9.596, 7.596, 5.596, 1.596, 0.)
	}
	else if(grepl("S8p",infilexml.opt))
	{
		beast2.spec		<- hivc.beast2.get.specifications(mcmc.length=beast.mcmc.length, bdsky.intervalNumber=8)
		beast2.spec$bdsky.sprop.changepoint.value	<- beast2.spec$bdsky.R0.changepoint.value		<- beast2.spec$bdsky.notInf.changepoint.value	<- c(9.596, 8.596, 7.596, 6.596, 5.596, 1.596, 0.596, 0.)
	}
	else if(grepl("s0106",infilexml.opt))
	{
		beast2.spec		<- hivc.beast2.get.specifications(mcmc.length=beast.mcmc.length, bdsky.intervalNumber=5)
		beast2.spec$bdsky.sprop.changepoint.value	<- beast2.spec$bdsky.R0.changepoint.value		<- beast2.spec$bdsky.notInf.changepoint.value	<- c(9.596, 5.596, 1.596, 0.596, 0.)
		beast2.spec$bdsky.sprop.value				<- c(0.1, 0.2, 0.7, 0.6, 0.6)
		beast2.spec$bdsky.sprop.prior				<- c("Exponential/0.1/0","Uniform/0.1/1.0","Uniform/0.6/1.0","Uniform/0.5/1.0","Uniform/0.5/1.0")	
	}
	else if(grepl("s00106",infilexml.opt))
	{
		beast2.spec		<- hivc.beast2.get.specifications(mcmc.length=beast.mcmc.length, bdsky.intervalNumber=5)
		beast2.spec$bdsky.sprop.changepoint.value	<- beast2.spec$bdsky.R0.changepoint.value		<- beast2.spec$bdsky.notInf.changepoint.value	<- c(9.596, 5.596, 1.596, 0.596, 0.)
		beast2.spec$bdsky.sprop.value				<- c(0.01, 0.2, 0.7, 0.6, 0.6)
		beast2.spec$bdsky.sprop.prior				<- c("Exponential/0.01/0","Uniform/0.1/1.0","Uniform/0.6/1.0","Uniform/0.5/1.0","Uniform/0.5/1.0")	
	}	
	else if(grepl("s0108",infilexml.opt))
	{
		beast2.spec		<- hivc.beast2.get.specifications(mcmc.length=beast.mcmc.length, bdsky.intervalNumber=5)
		beast2.spec$bdsky.sprop.changepoint.value	<- beast2.spec$bdsky.R0.changepoint.value		<- beast2.spec$bdsky.notInf.changepoint.value	<- c(9.596, 5.596, 1.596, 0.596, 0.)
		beast2.spec$bdsky.sprop.value				<- c(0.1, 0.5, 0.9, 0.8, 0.8)
		beast2.spec$bdsky.sprop.prior				<- c("Exponential/0.01/0","Uniform/0.4/1.0","Uniform/0.8/1.0","Uniform/0.7/1.0","Uniform/0.7/1.0")	
	}
	else if(grepl("s124",infilexml.opt))
	{
		beast2.spec		<- hivc.beast2.get.specifications(mcmc.length=beast.mcmc.length, bdsky.intervalNumber=5)
		beast2.spec$bdsky.sprop.changepoint.value	<- beast2.spec$bdsky.R0.changepoint.value		<- beast2.spec$bdsky.notInf.changepoint.value	<- c(9.596, 5.596, 1.596, 0.596, 0.)
		beast2.spec$bdsky.sprop.value				<- c(0.1, 0.5, 0.9, 0.8, 0.8)
		beast2.spec$bdsky.sprop.prior				<- c("Exponential/0.01/0","Exponential/0.1/0","Uniform/0.8/1.0","Uniform/0.2/1.0","Beta/2.5/4.0/0")	
	}	
	else if(grepl("s424",infilexml.opt))
	{
		beast2.spec		<- hivc.beast2.get.specifications(mcmc.length=beast.mcmc.length, bdsky.intervalNumber=5)
		beast2.spec$bdsky.sprop.changepoint.value	<- beast2.spec$bdsky.R0.changepoint.value		<- beast2.spec$bdsky.notInf.changepoint.value	<- c(9.596, 5.596, 1.596, 0.596, 0.)
		beast2.spec$bdsky.sprop.value				<- c(0.1, 0.5, 0.9, 0.8, 0.8)
		beast2.spec$bdsky.sprop.prior				<- c("Exponential/0.01/0","Exponential/0.4/0","Uniform/0.8/1.0","Uniform/0.2/1.0","Beta/2.5/4.0/0")	
	}
	else if(grepl("s024",infilexml.opt))
	{
		beast2.spec		<- hivc.beast2.get.specifications(mcmc.length=beast.mcmc.length, bdsky.intervalNumber=5)
		beast2.spec$bdsky.sprop.changepoint.value	<- beast2.spec$bdsky.R0.changepoint.value		<- beast2.spec$bdsky.notInf.changepoint.value	<- c(9.596, 5.596, 1.596, 0.596, 0.)
		beast2.spec$bdsky.sprop.value				<- c(0.1, 0.5, 0.9, 0.8, 0.8)
		beast2.spec$bdsky.sprop.prior				<- c("Exponential/0.01/0","Exponential/0.01/0","Uniform/0.8/1.0","Uniform/0.2/1.0","Beta/2.5/4.0/0")	
	}
	else if(grepl("s184",infilexml.opt))
	{
		beast2.spec		<- hivc.beast2.get.specifications(mcmc.length=beast.mcmc.length, bdsky.intervalNumber=5)
		beast2.spec$bdsky.sprop.changepoint.value	<- beast2.spec$bdsky.R0.changepoint.value		<- beast2.spec$bdsky.notInf.changepoint.value	<- c(9.596, 5.596, 1.596, 0.596, 0.)
		beast2.spec$bdsky.sprop.value				<- c(0.1, 0.5, 0.9, 0.6, 0.3)
		beast2.spec$bdsky.sprop.prior				<- c("Exponential/0.01/0","Exponential/0.1/0","Uniform/0.2/1.0","Beta/4.0/3.0/0","Beta/2.5/4.0/0")	
	}
	else if(grepl("sartest",infilexml.opt))
	{
		beast2.spec		<- hivc.beast2.get.specifications(mcmc.length=beast.mcmc.length, bdsky.intervalNumber=5)
		beast2.spec$bdsky.sprop.changepoint.value	<- beast2.spec$bdsky.R0.changepoint.value		<- beast2.spec$bdsky.notInf.changepoint.value	<- c(9.596, 5.596, 1.596, 0.596, 0.)
		beast2.spec$bdsky.sprop.value				<- c(0.1, 0.5, 0.9, 0.6, 0.3)
		beast2.spec$bdsky.sprop.prior				<- c("Exponential/0.01/0","Exponential/0.1/0","Uniform/0.2/1.0","Beta/4.0/3.0/0","Beta/2.5/4.0/0")
		beast2.spec$sasky.r.value					<- c(0.1, 0.2, 0.5, 0.7, 0.7)
		beast2.spec$sasky.r.prior					<- c("Uniform/0.0/0.5","Uniform/0.0/0.5","Uniform/0.0/1.0","Uniform/0.5/0.8","Uniform/0.5/1.0")
	}
	else if(grepl("d999",infilexml.opt))
	{
		beast2.spec		<- hivc.beast2.get.specifications(mcmc.length=beast.mcmc.length, bdsky.intervalNumber=5)
		beast2.spec$bdsky.sprop.changepoint.value	<- beast2.spec$bdsky.R0.changepoint.value		<- beast2.spec$bdsky.notInf.changepoint.value	<- c(9.596, 5.596, 1.596, 0.596, 0.)
		beast2.spec$bdsky.sprop.value				<- c(0.1, 0.5, 0.9, 0.6, 0.3)
		beast2.spec$bdsky.sprop.prior				<- c("Exponential/0.01/0","Exponential/0.1/0","Uniform/0.2/1.0","Beta/4.0/3.0/0","Beta/2.5/4.0/0")
		beast2.spec$bdsky.notInf.value				<- 1/c(9, 9, 9, 9, 9)
		beast2.spec$bdsky.notInf.prior				<- c("Exponential/0.11/0","Exponential/0.11/0","Exponential/0.11/0","Exponential/0.11/0","Exponential/0.11/0")
	}
	else if(grepl("d774",infilexml.opt))
	{
		beast2.spec		<- hivc.beast2.get.specifications(mcmc.length=beast.mcmc.length, bdsky.intervalNumber=5)
		beast2.spec$bdsky.sprop.changepoint.value	<- beast2.spec$bdsky.R0.changepoint.value		<- beast2.spec$bdsky.notInf.changepoint.value	<- c(9.596, 5.596, 1.596, 0.596, 0.)
		beast2.spec$bdsky.sprop.value				<- c(0.1, 0.5, 0.9, 0.6, 0.3)
		beast2.spec$bdsky.sprop.prior				<- c("Exponential/0.01/0","Exponential/0.1/0","Uniform/0.2/1.0","Beta/4.0/3.0/0","Beta/2.5/4.0/0")
		beast2.spec$bdsky.notInf.value				<- 1/c(7, 7, 7, 4, 4)
		beast2.spec$bdsky.notInf.prior				<- c("Exponential/0.14/0","Exponential/0.14/0","Exponential/0.14/0","Exponential/0.25/0","Exponential/0.25/0")
	}
	else if(grepl("d543",infilexml.opt))
	{
		beast2.spec		<- hivc.beast2.get.specifications(mcmc.length=beast.mcmc.length, bdsky.intervalNumber=5)
		beast2.spec$bdsky.sprop.changepoint.value	<- beast2.spec$bdsky.R0.changepoint.value		<- beast2.spec$bdsky.notInf.changepoint.value	<- c(9.596, 5.596, 1.596, 0.596, 0.)
		beast2.spec$bdsky.sprop.value				<- c(0.1, 0.5, 0.9, 0.6, 0.3)
		beast2.spec$bdsky.sprop.prior				<- c("Exponential/0.01/0","Exponential/0.1/0","Uniform/0.2/1.0","Beta/4.0/3.0/0","Beta/2.5/4.0/0")
		beast2.spec$bdsky.notInf.value				<- 1/c(5, 4, 4, 3, 3)
		beast2.spec$bdsky.notInf.prior				<- c("Exponential/0.2/0","Exponential/0.25/0","Exponential/0.25/0","Exponential/0.33/0","Exponential/0.33/0")
	}
	else if(grepl("dg543",infilexml.opt))
	{
		beast2.spec		<- hivc.beast2.get.specifications(mcmc.length=beast.mcmc.length, bdsky.intervalNumber=5)
		beast2.spec$bdsky.sprop.changepoint.value	<- beast2.spec$bdsky.R0.changepoint.value		<- beast2.spec$bdsky.notInf.changepoint.value	<- c(9.596, 5.596, 1.596, 0.596, 0.)
		beast2.spec$bdsky.sprop.value				<- c(0.1, 0.5, 0.9, 0.6, 0.3)
		beast2.spec$bdsky.sprop.prior				<- c("Exponential/0.01/0","Exponential/0.1/0","Uniform/0.2/1.0","Beta/4.0/3.0/0","Beta/2.5/4.0/0")
		beast2.spec$bdsky.notInf.value				<- 1/c(5, 4, 4, 3, 3)
		beast2.spec$bdsky.notInf.prior				<- c("Gamma/5/0.03/0.1","Gamma/5/0.03/0.1","Exponential/0.25/0","Gamma/5/0.05/0.1","Exponential/0.33/0")
		beast2.spec$sasky.r.value					<- c(0.9, 0.5, 0.5, 0.5, 0.5)
		beast2.spec$sasky.r.prior					<- c("Beta/2/1/0","Uniform/0.0/1.0","Uniform/0.0/1.0","Uniform/0.0/1.0","Uniform/0.0/1.0")
	}
	else if(grepl("r7543",infilexml.opt))
	{
		beast2.spec		<- hivc.beast2.get.specifications(mcmc.length=beast.mcmc.length, bdsky.intervalNumber=5)
		beast2.spec$bdsky.sprop.changepoint.value	<- beast2.spec$bdsky.R0.changepoint.value		<- beast2.spec$bdsky.notInf.changepoint.value	<- c(9.596, 5.596, 1.596, 0.596, 0.)
		beast2.spec$bdsky.sprop.value				<- c(0.1, 0.5, 0.9, 0.6, 0.3)
		beast2.spec$bdsky.sprop.prior				<- c("Exponential/0.01/0","Exponential/0.1/0","Uniform/0.2/1.0","Beta/4.0/3.0/0","Beta/2.5/4.0/0")
		beast2.spec$bdsky.notInf.value				<- 1/c(5, 4, 4, 3, 3)
		beast2.spec$bdsky.notInf.prior				<- c("Gamma/5/0.03/0.1","Gamma/5/0.03/0.1","Exponential/0.25/0","Gamma/5/0.05/0.1","Exponential/0.33/0")
		beast2.spec$sasky.r.value					<- c(0.5, 0.5, 0.5, 0.5, 0.5)
		beast2.spec$sasky.r.prior					<- c("Uniform/0.0/0.7","Uniform/0.0/0.7","Uniform/0.0/0.7","Uniform/0.0/0.7","Uniform/0.0/0.7")
	}
	else if(grepl("r5543",infilexml.opt))
	{
		beast2.spec		<- hivc.beast2.get.specifications(mcmc.length=beast.mcmc.length, bdsky.intervalNumber=5)
		beast2.spec$bdsky.sprop.changepoint.value	<- beast2.spec$bdsky.R0.changepoint.value		<- beast2.spec$bdsky.notInf.changepoint.value	<- c(9.596, 5.596, 1.596, 0.596, 0.)
		beast2.spec$bdsky.sprop.value				<- c(0.1, 0.5, 0.9, 0.6, 0.3)
		beast2.spec$bdsky.sprop.prior				<- c("Exponential/0.01/0","Exponential/0.1/0","Uniform/0.2/1.0","Beta/4.0/3.0/0","Beta/2.5/4.0/0")
		beast2.spec$bdsky.notInf.value				<- 1/c(5, 4, 4, 3, 3)
		beast2.spec$bdsky.notInf.prior				<- c("Gamma/5/0.03/0.1","Gamma/5/0.03/0.1","Exponential/0.25/0","Gamma/5/0.05/0.1","Exponential/0.33/0")
		beast2.spec$sasky.r.value					<- c(0.4, 0.4, 0.4, 0.4, 0.4)
		beast2.spec$sasky.r.prior					<- c("Uniform/0.0/0.5","Uniform/0.0/0.5","Uniform/0.0/0.5","Uniform/0.0/0.5","Uniform/0.0/0.5")
	}
	else if(grepl("r1543",infilexml.opt))
	{
		beast2.spec		<- hivc.beast2.get.specifications(mcmc.length=beast.mcmc.length, bdsky.intervalNumber=5)
		beast2.spec$bdsky.sprop.changepoint.value	<- beast2.spec$bdsky.R0.changepoint.value		<- beast2.spec$bdsky.notInf.changepoint.value	<- c(9.596, 5.596, 1.596, 0.596, 0.)
		beast2.spec$bdsky.sprop.value				<- c(0.1, 0.5, 0.9, 0.6, 0.3)
		beast2.spec$bdsky.sprop.prior				<- c("Exponential/0.01/0","Exponential/0.1/0","Uniform/0.2/1.0","Beta/4.0/3.0/0","Beta/2.5/4.0/0")
		beast2.spec$bdsky.notInf.value				<- 1/c(5, 4, 4, 3, 3)
		beast2.spec$bdsky.notInf.prior				<- c("Gamma/5/0.03/0.1","Gamma/5/0.03/0.1","Exponential/0.25/0","Gamma/5/0.05/0.1","Exponential/0.33/0")
		beast2.spec$sasky.r.value					<- c(0.05, 0.05, 0.05, 0.05, 0.05)
		beast2.spec$sasky.r.prior					<- c("Uniform/0.0/0.1","Uniform/0.0/0.1","Uniform/0.0/0.1","Uniform/0.0/0.1","Uniform/0.0/0.1")
	}
	else if(grepl("ori",infilexml.opt))
	{
		beast2.spec		<- hivc.beast2.get.specifications(mcmc.length=beast.mcmc.length, bdsky.intervalNumber=5)
		beast2.spec$bdsky.sprop.changepoint.value	<- beast2.spec$bdsky.R0.changepoint.value		<- beast2.spec$bdsky.notInf.changepoint.value	<- c(9.596, 5.596, 1.596, 0.596, 0.)
		beast2.spec$bdsky.sprop.value				<- c(0.1, 0.1, 0.6, 0.6, 0.3)
		beast2.spec$bdsky.sprop.prior				<- c("Exponential/0.01/0","Exponential/0.1/0","Beta/4.0/3.0/0","Beta/4.0/3.0/0","Beta/2.5/4.0/0")
		beast2.spec$bdsky.notInf.value				<- 1/c(5, 4, 4, 3, 3)
		beast2.spec$bdsky.notInf.prior				<- c("LogNormal/0.2/0.6/0.1/true","LogNormal/0.25/0.8/0.1/true","LogNormal/0.3/0.8/0.1/true","LogNormal/0.5/1.2/0.1/true","LogNormal/0.5/1.2/0.1/true")
		beast2.spec$sasky.r.value					<- rep(0.05, 5)
		beast2.spec$sasky.r.prior					<- rep("Uniform/0.0/0.1",5)		
		beast2.spec$bdsky.origin.prior				<- as.numeric(substr(infilexml.opt,4,nchar(infilexml.opt)))		
		beast2.spec$bdsky.origin.prior				<- paste("Uniform/",min(20,beast2.spec$bdsky.origin.prior-20),"/",beast2.spec$bdsky.origin.prior,sep='')
print(beast2.spec$bdsky.origin.prior	)
		stop()
		
	}	
	else stop("unknown infilexml.opt")
	#		 
	#
	#
	#
	if(grepl("sasky",infilexml.template))
	{
		beast2.spec$treemodel			<- "SampledAncestorSkylineModel"
		prog.beast						<- PR.BEAST2SA
	}
	if(grepl("bdsky",infilexml.template))
	{
		beast2.spec$treemodel			<- "BirthDeathSkylineModel"
		prog.beast						<- PR.BEAST2
	}
	#	
	if(verbose)
	{
		print(indir)
		print(infile)
		print(insignat)
		print(indircov)
		print(infilecov)
		print(infiletree)
		print(infilexml)
		print(outdir)		
		print(outsignat)
		print(resume)
		print(opt.brl)
		print(thresh.brl)
		print(thresh.bs)
		print(pool.ntip)		 
		print(infilexml.opt)
		print(infilexml.template)
	}	
	#
	#	load complete tree to generate starting tree
	#
	#file					<- paste(indir,'/',infiletree,'_',gsub('/',':',insignat),".R",sep='')
	#if(verbose)	cat(paste("\nload complete tree to generate starting tree from file",file))
	#load(file)	#load object 'ph'	
	#
	#	load sequences
	#
	file		<- paste(indir,'/',infile,'_',gsub('/',':',insignat),".R",sep='')
	if(verbose)	cat(paste("\nload sequences from file",file))
	load( file )
	#print( seq.PROT.RT )	
	#
	#	load msm clusters
	#
	argv		<<- hivc.cmd.clustering.msm(indir, infiletree, insignat, indircov, infilecov, opt.brl, thresh.brl, thresh.bs, resume=resume)
	argv		<<- unlist(strsplit(argv,' '))		
	msm			<- hivc.prog.get.clustering.MSM()	
	#
	#	select seroconverters
	#
	df.cluinfo				<- msm$df.cluinfo
	tmp						<- df.cluinfo[,	list(clu.bwpat.medbrl=clu.bwpat.medbrl[1],clu.npat=clu.npat[1], clu.fPossAcute=clu.fPossAcute[1], fNegT=length(which(!is.na(NegT))) / clu.ntip[1]),by="cluster"]										
	tmp						<- subset(tmp, fNegT>=quantile(tmp[,fNegT], probs=0.8) )
	cluphy.df				<- merge( subset(tmp,select=cluster), df.cluinfo, all.x=1, by="cluster" )
	if(verbose) cat(paste("\nnumber of selected sequences is n=",nrow(cluphy.df)))
	cluphy.df				<- hivc.beast.addBEASTLabel( cluphy.df )
	#
	#	create sets of cluster pools for BEAST
	#
	df.clupool				<- hivc.beast.poolclusters(cluphy.df, pool.ntip= pool.ntip, verbose=1)
	#
	#	load xml template file
	#	
	file			<- paste(CODE.HOME,"/data/BEAST2_template_",infilexml.template,".xml",sep='')
	bxml.template	<- xmlTreeParse(file, useInternalNodes=TRUE, addFinalizer = TRUE)			
	#
	#	create BEAST2 XML file	
	#
	bfile			<- lapply(seq_len(length(df.clupool$pool.df)), function(pool.id)
						{
							df							<- df.clupool$pool.df[[pool.id]]
							setkey(df, cluster)
							#outfile				<- paste(indir,'/',outfile,".nex",sep='')
							#tmp					<- seq.PROT.RT[df[, FASTASampleCode],]		
							#rownames(tmp)		<- df[, BEASTlabel]
							#dummy				<- hivc.seq.write.dna.nexus(tmp, file=outfile )
							beast2.spec$xml.dir			<- indir
							beast2.spec$xml.filename	<- paste(infilexml,'-',pool.ntip,'-',pool.id,'_',infilexml.template,'_',infilexml.opt,'_',gsub('/',':',outsignat),sep='')
							beast2.xml					<- hivc.beast2.get.xml( bxml.template, seq.PROT.RT, df, beast2.spec, ph=NULL, verbose=1)			
							file						<- paste(beast2.spec$xml.dir,'/',beast2.spec$xml.filename,".xml", sep='')
							if(verbose)	cat(paste("\nwrite xml file to",file))
							saveXML(beast2.xml, file=file)
							paste(infilexml,'-',pool.ntip,'-',pool.id,'_',infilexml.template,'_',infilexml.opt,sep='')
						})
	#
	#	generate BEAST commands and run
	#
	sapply(bfile, function(x)
			{
				cmd			<- hivc.cmd.beast2.runxml(indir, x, outsignat, hpc.ncpu=hpc.ncpu, prog.beast=prog.beast, prog.opt.Xmx="1200m", hpc.tmpdir.prefix="beast2")
				#cmd		<- paste(cmd,hivc.cmd.beast.evalrun(outdir, infilexml, outsignat, infilexml.opt, infilexml.template, length(bfile), verbose=1),sep='')				
				cmd			<- cmd.hpcwrapper(cmd, hpc.walltime=hpc.walltime, hpc.q="pqeph", hpc.mem=hpc.mem,  hpc.nproc=hpc.ncpu)					
				cat(cmd)				
				outfile		<- paste("b2",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='.')
				cmd.hpccaller(outdir, outfile, cmd)
				stop()
			})	
}
######################################################################################
hivc.prog.BEAST.generate.xml<- function()
{	
	require(XML)
	
	indir				<- paste(DATA,"tmp",sep='/')		
	infile				<- "ATHENA_2013_03_NoDRAll+LANL_Sequences"		
	insignat			<- "Thu_Aug_01_17/05/23_2013"
	indircov			<- paste(DATA,"derived",sep='/')
	infilecov			<- "ATHENA_2013_03_AllSeqPatientCovariates"
	infiletree			<- paste(infile,"examlbs100",sep="_")
	infilexml			<- paste(infile,'_',"beast",'_',"seroneg",sep='')
	
	outdir				<- indir
	outsignat			<- "Tue_Aug_26_09/13/47_2013"
		
	opt.brl					<- "dist.brl.casc" 
	thresh.brl				<- 0.096
	thresh.bs				<- 0.8
	pool.ntip				<- 130
	beast.mcmc.chainLength	<- 100000000
	infilexml.opt			<- "mph4cluLdTd"
	infilexml.template		<- "standard"	
	
	resume				<- 1
	verbose				<- 1
	hpc.walltime		<- 171
	hpc.ncpu			<- 4
	hpc.mem				<- "1800mb"

	
	if(exists("argv"))
	{
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,6),
									indir= return(substr(arg,8,nchar(arg))),NA)	}))
		if(length(tmp)>0) indir<- tmp[1]		
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,7),
									infile= return(substr(arg,9,nchar(arg))),NA)	}))
		if(length(tmp)>0) infile<- tmp[1]				
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,9),
									insignat= return(substr(arg,11,nchar(arg))),NA)	}))
		if(length(tmp)>0) insignat<- tmp[1]		
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,9),
									indircov= return(substr(arg,11,nchar(arg))),NA)	}))
		if(length(tmp)>0) indircov<- tmp[1]		
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,10),
									infilecov= return(substr(arg,12,nchar(arg))),NA)	}))
		if(length(tmp)>0) infilecov<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,11),
									infiletree= return(substr(arg,13,nchar(arg))),NA)	}))
		if(length(tmp)>0) infiletree<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,10),
									infilexml= return(substr(arg,12,nchar(arg))),NA)	}))
		if(length(tmp)>0) infilexml<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,7),
									outdir= return(substr(arg,9,nchar(arg))),NA)	}))
		if(length(tmp)>0) outdir<- tmp[1]		
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,10),
									outsignat= return(substr(arg,12,nchar(arg))),NA)	}))
		if(length(tmp)>0) outsignat<- tmp[1]				
		#		
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,7),
									resume= return(as.numeric(substr(arg,9,nchar(arg)))),NA)	}))
		if(length(tmp)>0) resume<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,2),
									v= return(as.numeric(substr(arg,4,nchar(arg)))),NA)	}))
		if(length(tmp)>0) verbose<- tmp[1]
		#
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,8),
									opt.brl= return(substr(arg,10,nchar(arg))),NA)	}))
		if(length(tmp)>0) opt.brl<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,10),
									thresh.bs= return(as.numeric(substr(arg,12,nchar(arg)))),NA)	}))
		if(length(tmp)>0) thresh.bs<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,11),
									thresh.brl= return(as.numeric(substr(arg,13,nchar(arg)))),NA)	}))
		if(length(tmp)>0) thresh.brl<- tmp[1]
		#
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,10),
									pool.ntip= return(as.numeric(substr(arg,12,nchar(arg)))),NA)	}))
		if(length(tmp)>0) pool.ntip<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,14),
									infilexml.opt= return(substr(arg,16,nchar(arg))),NA)	}))
		if(length(tmp)>0) infilexml.opt<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,19),
									infilexml.template= return(substr(arg,21,nchar(arg))),NA)	}))
		if(length(tmp)>0) infilexml.template<- tmp[1]
	}	
	
	#	load complete tree to generate starting tree
	file					<- paste(indir,'/',infiletree,'_',gsub('/',':',insignat),".R",sep='')
	if(verbose)	cat(paste("\nload complete tree to generate starting tree from file",file))
	load(file)	#load object 'ph'
	
	xml.monophyly4clusters			<- ifelse(grepl("mph4clu",infilexml.opt),1,0)		
	pool.includealwaysbeforeyear	<- ifelse(grepl("fx03",infilexml.opt),2003, NA)
	xml.prior4tipstem				<- ifelse(grepl("u4tip",infilexml.opt),"uniform",NA)
	xml.taxon4tipstem				<- ifelse(	any(sapply(c("u4tip","tx4tip"), function(x) grepl(x,infilexml.opt))),	1, 0)
	df.resetTipDate					<- NA
	if(grepl("LdTd",infilexml.opt))				df.resetTipDate<- "LdTd"
	else if(grepl("LsTd",infilexml.opt))		df.resetTipDate<- "LsTd"
	else if(grepl("UmTd",infilexml.opt))		df.resetTipDate<- "UmTd"
	#case NoTd is dealt with below
	#
	#
	if(verbose)
	{
		print(indir)
		print(infile)
		print(insignat)
		print(indircov)
		print(infilecov)
		print(infiletree)
		print(infilexml)
		print(outdir)		
		print(outsignat)
		print(resume)
		print(opt.brl)
		print(thresh.brl)
		print(thresh.bs)
		print(pool.ntip)
		print(pool.includealwaysbeforeyear)
		print(infilexml.opt)
		print(infilexml.template)
		print(xml.monophyly4clusters)
		print(xml.taxon4tipstem)
		print(xml.prior4tipstem)
		print(df.resetTipDate)
	}	
	#
	#	load sequences
	#
	file		<- paste(indir,'/',infile,'_',gsub('/',':',insignat),".R",sep='')
	if(verbose)	cat(paste("\nload sequences from file",file))
	load( file )
	#print( seq.PROT.RT )	
	#
	#	load msm clusters
	#
	argv		<<- hivc.cmd.clustering.msm(indir, infiletree, insignat, indircov, infilecov, opt.brl, thresh.brl, thresh.bs, resume=resume)
	argv		<<- unlist(strsplit(argv,' '))		
	msm			<- hivc.prog.get.clustering.MSM()	
	#
	#	select seroconverters
	#
	df.cluinfo				<- msm$df.cluinfo
	tmp						<- df.cluinfo[,	list(clu.bwpat.medbrl=clu.bwpat.medbrl[1],clu.npat=clu.npat[1], clu.fPossAcute=clu.fPossAcute[1], fNegT=length(which(!is.na(NegT))) / clu.ntip[1]),by="cluster"]										
	tmp						<- subset(tmp, fNegT>=quantile(tmp[,fNegT], probs=0.8) )
	cluphy.df				<- merge( subset(tmp,select=cluster), df.cluinfo, all.x=1, by="cluster" )
	if(verbose) cat(paste("\nnumber of selected sequences is n=",nrow(cluphy.df)))
	cluphy.df				<- hivc.beast.addBEASTLabel( cluphy.df, df.resetTipDate=df.resetTipDate )
	#
	#	create sets of cluster pools for BEAST
	#
	df.clupool				<- hivc.beast.poolclusters(cluphy.df, pool.ntip= pool.ntip, pool.includealwaysbeforeyear=pool.includealwaysbeforeyear, verbose=1)			
	#
	if(0)		#used to generate standard xml file
	{
		outfile		<- paste(infile,"beast","seroneg",sep='_')
		outsignat	<- "Thu_Aug_01_17/05/23_2013"
		hivc.beast.writeNexus4Beauti(seq.PROT.RT, cluphy.df, file=paste(outdir,'/',outfile,'_',gsub('/',':',outsignat),".nex",sep=''))
	}
	#
	#	load xml template file
	#	
	file		<- paste(indir,'/',infilexml,'_',infilexml.template,'_',gsub('/',':',insignat),".xml",sep='')
	btemplate	<- xmlTreeParse(file, useInternalNodes=TRUE, addFinalizer = TRUE)		
	#
	#	produce xml files for each cluster pool from template
	#	
	bfile		<- lapply(seq_len(length(df.clupool$pool.df)), function(pool.id)
			{
				df					<- df.clupool$pool.df[[pool.id]]
				setkey(df, cluster)
				#print( unique(df[,cluster]) )
				#	get xml file 
				outfile				<- paste(infilexml,'-',pool.ntip,'-',pool.id,'_',infilexml.template,'_',infilexml.opt,'_',gsub('/',':',outsignat),sep='')
				beast.label.datepos	<- ifelse(!is.na(df.resetTipDate), 5, 4)
				beast.usingDates	<- ifelse(grepl("NoTd",infilexml.opt), "false", "true")
				bxml				<- hivc.beast.get.xml(btemplate, seq.PROT.RT, df, outfile, ph=ph, xml.monophyly4clusters=xml.monophyly4clusters, xml.taxon4tipstem=xml.taxon4tipstem, xml.prior4tipstem=xml.prior4tipstem, beast.label.datepos=beast.label.datepos, beast.label.sep= '_', beast.date.direction= "forwards", beast.date.units= "years", beast.usingDates=beast.usingDates, beast.mcmc.chainLength=beast.mcmc.chainLength, verbose=1)
				getNodeSet(bxml, "//*[@id='tmrca(c1)']")		
				#	write xml file
				file				<- paste(outdir,'/',outfile,".xml",sep='')
				if(verbose)	cat(paste("\nwrite xml file to",file))
				saveXML(bxml, file=file)		
				#	freed through R garbage collection (I hope!) since addFinalizer=TRUE throughout
				paste(infilexml,'-',pool.ntip,'-',pool.id,'_',infilexml.template,'_',infilexml.opt,sep='')
			})
	#
	#	generate BEAST commands and run
	#
	sapply(bfile, function(x)
			{
				cmd			<- hivc.cmd.beast.runxml(outdir, x, outsignat, hpc.tmpdir.prefix="beast", hpc.ncpu=hpc.ncpu)				
				cmd			<- paste(cmd,hivc.cmd.beast.evalrun(outdir, infilexml, outsignat, infilexml.opt, infilexml.template, length(bfile), verbose=1),sep='')				
				cmd			<- cmd.hpcwrapper(cmd, hpc.walltime=hpc.walltime, hpc.q="pqeph", hpc.mem=hpc.mem,  hpc.nproc=hpc.ncpu)					
				cat(cmd)
				stop()
				outfile		<- paste("bea",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='.')
				cmd.hpccaller(outdir, outfile, cmd)				
			})		
}
######################################################################################
hivc.prog.get.geneticdist<- function()
{
	library(bigmemory)
	library(ape)
	
	indir		<- outdir<- paste(DATA,"tmp",sep='/')
	infile		<- "ATHENA_2013_03_FirstAliSequences_PROTRT"
	resume		<- verbose <- 1	
	signat 		<- "Wed_May__1_17/08/15_2013"
	gd.max		<- NA
	out.phylip	<- 0
	if(exists("argv"))
	{
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,6),
									indir= return(substr(arg,8,nchar(arg))),NA)	}))
		if(length(tmp)>0) indir<- tmp[1]		
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,7),
									infile= return(substr(arg,9,nchar(arg))),NA)	}))
		if(length(tmp)>0) infile<- tmp[1]				
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,7),
									outdir= return(substr(arg,9,nchar(arg))),NA)	}))		
		if(length(tmp)>0) outdir<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,7),
									signat= return(substr(arg,9,nchar(arg))),NA)	}))
		if(length(tmp)>0) signat<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,7),
									resume= return(as.numeric(substr(arg,9,nchar(arg)))),NA)	}))
		if(length(tmp)>0) resume<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,2),
									v= return(as.numeric(substr(arg,4,nchar(arg)))),NA)	}))
		if(length(tmp)>0) verbose<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,6),
									maxgd= return(as.numeric(substr(arg,8,nchar(arg)))),NA)	}))
		if(length(tmp)>0) gd.max<- tmp[1]
	}
	if(1)
	{
		print(indir)
		print(infile)
		print(outdir)
		print(signat)
		print(resume)
		print(gd.max)
	}	
	
	tmp			<- list.files(outdir, pattern=paste(".gdm$",sep=''))		
	file		<- tmp[ grepl(paste(infile,'_',sep=''), tmp) & grepl(gsub('/',':',signat), tmp) ]
	if(verbose)	cat(paste("\nFound .gdm files matching input args, n=", length(file)))	
	if(resume)
	{
		if(verbose)	cat(paste("\nloading",file))
		gd.bigmat			<- read.big.matrix(file, has.row.names=1, sep=',', type="char")		
	}
	if(!resume || !length(file))
	{		
		if(verbose)	cat(paste("\ncreate gdm file"))
		file				<- paste(indir,"/",infile,"_",gsub('/',':',signat),".R",sep='')
		if(verbose)	cat(paste("\nload",file))
		load(file)
		str(seq.PROT.RT)		
			#tmp				<- tmp[1:10,]
		gd.bigmat			<- hivc.seq.dist(  seq.PROT.RT )
		file				<- paste(outdir,"/",infile,"_",gsub('/',':',signat),".gdm",sep='')
		if(verbose) cat(paste("\nwrite to",file))
		write.big.matrix(gd.bigmat, file, row.names= 1, col.names=0, sep=',')		
		#gd.bigmat.d 		<- describe( gd.bigmat )
		#file				<- paste(outdir,"/ATHENA_2013_03_FirstAliSequences_PROTRT_",gsub('/',':',signat),".gdm",sep='')		
		#dput(gd.bigmat.d, file=file)
	}
	if(!is.na(gd.max))
	{
		if(verbose)	cat(paste("\nselect sequences with at least one other sequence within gd.max",gd.max))
		#	now have pairwise genetic distances in gd.bigmat
		gd.bigmat.min	<- sapply(seq.int(1,nrow(gd.bigmat)-1),function(i)
								{
									tmp<- which.min(gd.bigmat[i,])
									c(i,tmp, gd.bigmat[i,tmp])
								})
		rownames(gd.bigmat.min)<- c("seq.idx1","seq.idx2","gd")
		gd.bigmat.seq.i	<- which( gd.bigmat.min["gd",] < (gd.max * 1e3) )
		gd.seqs			<- gd.bigmat.min[c("seq.idx1","seq.idx2"),gd.bigmat.seq.i]
		gd.seqs			<- sort(unique(as.vector(gd.seqs)))
		if(verbose) cat(paste("\ncomputed sequences sth at least one pair with less than ",gd.max*100,"% genetic distance, n=",length(gd.seqs)))		
		file			<- paste(indir,"/",infile,"_",gsub('/',':',signat),".R",sep='')
		if(verbose)	cat(paste("\nload",file))
		load(file)
		seq.PROT.RT.gd	<- seq.PROT.RT[ rownames(gd.bigmat)[gd.seqs], ]
		#	save selected sequences in R
		file 			<- paste(outdir,"/",infile,"_Gd",gd.max*1000,"_",gsub('/',':',signat),".R",sep='')		
		save(seq.PROT.RT.gd, file=file)		
		#	save selected sequences in phylip
		if(out.phylip)
		{
			file 			<- paste(outdir,"/",infile,"_Gd",gd.max*1000,"_",gsub('/',':',signat),".phylip",sep='')
			if(verbose) cat(paste("\nwrite to ",file))
			hivc.seq.write.dna.phylip(seq.PROT.RT.gd, file=file)
		}
	}
	else
	{
		seq.PROT.RT.gd<- NULL		
	}
	list(gd.bigmat=gd.bigmat, seq.PROT.RT.gd=seq.PROT.RT.gd)						
}

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

hivc.pipeline.ExaML<- function()
{
	dir.name<- DATA
	if(0)	#compute one ExaML tree, no bootstrapping
	{		
		indir	<- paste(dir.name,"tmp",sep='/')
		infile	<- "ATHENA_2013_03_FirstAliSequences_PROTRT"
		outdir	<- paste(dir.name,"tmp",sep='/')
		cmd		<- paste(cmd,cmd.examl(indir,infile,gsub('/',':',signat.out),gsub('/',':',signat.out),outdir=outdir,resume=1,verbose=1),sep='')
		cmd		<- paste(cmd,cmd.examl.cleanup(outdir),sep='')
	}
	if(0)	#compute ExaML trees with bootstrap values. Bootstrap is over initial starting trees to start ML search.
	{
		bs.from	<- 0
		bs.to	<- 0
		bs.n	<- 100
		signat.in	<- "Sat_May_11_14/23/46_2013"
		signat.out	<- "Sat_May_11_14/23/46_2013"				
		indir	<- paste(dir.name,"tmp",sep='/')
		infile	<- "ATHENA_2013_03_FirstCurSequences_PROTRT"
		infile	<- "ATHENA_2013_03_FirstCurSequences_PROTRTCD3"
		outdir	<- paste(dir.name,"tmp",sep='/')
		cmd		<- cmd.examl.bsstarttree(indir,infile,gsub('/',':',signat.out),gsub('/',':',signat.out),bs.from=bs.from,bs.to=bs.to,bs.n=bs.n,outdir=outdir, resume=1, verbose=1)
		#check if we have all bs.n files. if yes, combine and cleanup
		
		outdir	<- paste(dir.name,"tmp",sep='/')							
		lapply(cmd, function(x)
				{				
					x		<- cmd.hpcwrapper(x, hpc.walltime=36, hpc.q="pqeph")
					signat	<- paste(strsplit(date(),split=' ')[[1]],collapse='_',sep='')
					outfile	<- paste("pipeline",signat,sep='.')
					cat(x)
					#cmd.hpccaller(outdir, outfile, x)
					#Sys.sleep(1)
				})
		stop()
	}
	if(1)	#compute ExaML trees with bootstrap values. Bootstrap is over codon in alignment and over initial starting trees to start ML search.
	{
		bs.from		<- 0
		bs.to		<- 100
		bs.n		<- 500
		
		indir		<- paste(dir.name,"tmp",sep='/')
		#infile		<- "ATHENA_2013_03_CurAll+LANL_Sequences"
		#signat.in	<- "Sat_Jun_16_17:23:46_2013"								
		#infile		<- "ATHENA_2013_03_NoDRAll+LANL_Sequences"
		#signat.in	<- "Thu_Aug_01_17/05/23_2013"	
		#infile		<- "ATHENA_2013_03_NoRCDRAll+LANL_Sequences"	
		#signat.in	<- "Fri_Nov_01_16/07/23_2013"
		infile		<- "ATHENA_2013_03_-DR-RC-SH+LANL_Sequences"
		signat.in	<- "Wed_Dec_18_11/37/00_2013"
		
		#infile		<- "UKCA_2013_07_TNTPHIVnGTR"
		#signat.in	<- "Mon_Sep_22_17/23/46_2013"
		
		outdir		<- paste(dir.name,"tmp",sep='/')
		cmd			<- cmd.examl.bootstrap(indir,infile,gsub('/',':',signat.in),gsub('/',':',signat.in),bs.from=bs.from,bs.to=bs.to,bs.n=bs.n,outdir=outdir, resume=1, verbose=1)				
		outdir		<- paste(dir.name,"tmp",sep='/')							
		lapply(cmd, function(x)
				{				
					x		<- cmd.hpcwrapper(x, hpc.walltime=24, hpc.q=NA, hpc.mem="3850mb", hpc.nproc=8)
					#x		<- cmd.hpcwrapper(x, hpc.walltime=24, hpc.q="pqeph", hpc.mem="3850mb", hpc.nproc=8)
					signat	<- paste(strsplit(date(),split=' ')[[1]],collapse='_',sep='')
					outfile	<- paste("exa",signat,sep='.')
					#cat(x)
					cmd.hpccaller(outdir, outfile, x)
					Sys.sleep(1)
				})
		stop()
	}	
}

hivc.pipeline.clustering<- function()
{	
	if(1)	#clustering: precompute clustering objects, evaluate TPTN, get default clustering, refine to capture MSM transmission
	{	
		resume		<- 1
		verbose		<- 1
		#seq project
		indir		<- paste(DATA,"tmp",sep='/')		
		infile		<- "ATHENA_2013_03_CurAll+LANL_Sequences_examlbs100"
		insignat	<- "Sat_Jun_16_17/23/46_2013"
		infile		<- "ATHENA_2013_03_NoDRAll+LANL_Sequences_examlbs100"
		#infile		<- "ATHENA_2013_03_NoDRAll+LANL_Sequences_examlbs500"
		insignat	<- "Thu_Aug_01_17/05/23_2013"		
		infile		<- "ATHENA_2013_03_NoRCDRAll+LANL_Sequences_examlbs500"	
		insignat	<- "Fri_Nov_01_16/07/23_2013"
		
		
		#seq covariates
		indircov	<- paste(DATA,"derived",sep='/')
		infilecov	<- "ATHENA_2013_03_AllSeqPatientCovariates"
		#default clustering tptn		
		patient.n	<- 15700
		#default clustering	
		opt.brl		<- "dist.brl.casc" 
		thresh.brl	<- 0.096
		thresh.bs	<- 0.8		
		
		cmd			<- hivc.cmd.preclustering(indir, infile, insignat, indircov, infilecov)		
		cmd			<- paste(cmd, hivc.cmd.clustering.tptn(indir, infile, insignat, indircov, infilecov, opt.brl="dist.brl.casc", patient.n=patient.n, resume=resume),sep='')
		cmd			<- paste(cmd, hivc.cmd.clustering.tptn(indir, infile, insignat, indircov, infilecov, opt.brl="dist.brl.max", patient.n=patient.n, resume=resume),sep='')
		cmd			<- paste(cmd, hivc.cmd.clustering(indir, infile, insignat, opt.brl, thresh.brl, thresh.bs, resume=resume),sep='')		
		cmd			<- paste(cmd, hivc.cmd.clustering.msm(indir, infile, insignat, indircov, infilecov, opt.brl, thresh.brl, thresh.bs, resume=resume),sep='')			
		cmd			<- cmd.hpcwrapper(cmd, hpc.walltime=71, hpc.q="pqeph", hpc.mem="15850mb",  hpc.nproc=1)
		
		cat(cmd)
		#stop()
		outdir		<- paste(DATA,"tmp",sep='/')
		outfile		<- paste("clust",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='.')
		cmd.hpccaller(outdir, outfile, cmd)
		quit("no")
	}
}

hivc.pipeline.BEAST<- function()
{
	if(0)	#run BEAST 1.7.5 GMRF skyline
	{
		indir				<- paste(DATA,"tmp",sep='/')		
		infile				<- "ATHENA_2013_03_NoDRAll+LANL_Sequences"		
		insignat			<- "Thu_Aug_01_17/05/23_2013"
		indircov			<- paste(DATA,"derived",sep='/')
		infilecov			<- "ATHENA_2013_03_AllSeqPatientCovariates"
		infiletree			<- paste(infile,"examlbs100",sep="_")
		infilexml			<- paste(infile,'_',"beast",'_',"seroneg",sep='')
		#infilexml.template	<- "um22rhU2050"
		#infilexml.template	<- "um22rhG202018"
		#infilexml.template	<- "rhU65rho753"
		#infilexml.template	<- "rhU65rho903"
		#infilexml.template	<- "rhU65rho906"
		#infilexml.template	<- "rhU65rho909"	
		#infilexml.template	<- "um181rhU2045"
		#infilexml.template	<- "um182rhU2045"
		#infilexml.template	<- "um183rhU2045"
		#infilexml.template	<- "um182us45"
		#infilexml.template	<- "um182us60"
		#infilexml.template	<- "um182rhU2045ay"
		infilexml.template	<- "um232rhU2045"
		#infilexml.template	<- "um232rhU2045ay"
		#infilexml.opt		<- "txs4clu"
		#infilexml.opt		<- "txs4clufx03"
		#infilexml.opt		<- "mph4clu"
		#infilexml.opt		<- "mph4clutx4tip"
		#infilexml.opt		<- "mph4clufx03"
		#infilexml.opt		<- "mph4cluLdTd"
		#infilexml.opt		<- "mph4cluUmTd"
		#infilexml.opt		<- "mph4cluLsTd"
		#infilexml.opt		<- "mph4cluNoTd"
		#infilexml.opt		<- "mph4cluu4tipLdTd"
		infilexml.opt		<- "mph4clutx4tipLdTd"
		infilexml.opt		<- "mph4clutx4tipLsTd"
		#infilexml.opt		<- "mph4clutx4tip"
		
		outdir				<- indir
		outsignat			<- "Tue_Aug_26_09/13/47_2013"
		
		opt.brl				<- "dist.brl.casc" 
		thresh.brl			<- 0.096
		thresh.bs			<- 0.8
		pool.ntip			<- 130		
		#pool.ntip			<- 150
		#pool.ntip			<- 190
		#pool.ntip			<- 400
		resume				<- 1
		verbose				<- 1
		
		argv				<<- hivc.cmd.beast.poolrunxml(indir, infile, insignat, indircov, infilecov, infiletree, infilexml, outsignat, pool.ntip, infilexml.opt=infilexml.opt, infilexml.template=infilexml.template, opt.brl=opt.brl, thresh.brl=thresh.brl, thresh.bs=thresh.bs, resume=resume, verbose=1)
		argv				<<- unlist(strsplit(argv,' '))
		hivc.prog.BEAST.generate.xml()		
		quit("no")
	}
	if(1)		#generate BEAST2 BDSKYline xml file
	{
		indir				<- paste(DATA,"tmp",sep='/')		
		infile				<- "ATHENA_2013_03_NoDRAll+LANL_Sequences"		
		insignat			<- "Thu_Aug_01_17/05/23_2013"
		indircov			<- paste(DATA,"derived",sep='/')
		infilecov			<- "ATHENA_2013_03_AllSeqPatientCovariates"
		infiletree			<- paste(infile,"examlbs100",sep="_")
		infilexml			<- paste(infile,'_',"seroneg",sep='')
		outdir				<- indir
		outsignat			<- "Tue_Aug_26_09/13/47_2013"		
		opt.brl				<- "dist.brl.casc" 
		thresh.brl			<- 0.096
		thresh.bs			<- 0.8
		pool.ntip			<- 130		
		resume				<- 1
		verbose				<- 1
		
		infilexml.template	<- "bdsky_hky" 
		infilexml.template	<- "sasky_hky"
		infilexml.opt		<- "S4p"
		#infilexml.opt		<- "S5p"
		#infilexml.opt		<- "S8p"
		infilexml.opt		<- "s0106"
		#infilexml.opt		<- "s0108"
		#infilexml.opt		<- "s00106"
		infilexml.opt		<- "s124"
		infilexml.opt		<- "s024"
		infilexml.opt		<- "s424"
		infilexml.opt		<- "s184"
		infilexml.opt		<- "sartest"
		infilexml.opt		<- "d999"
		infilexml.opt		<- "d774"
		infilexml.opt		<- "d543"
		infilexml.opt		<- "dg543"
		infilexml.opt		<- "r7543"
		infilexml.opt		<- "r5543"
		infilexml.opt		<- "r1543"
		#infilexml.opt		<- "ori40"
		#infilexml.opt		<- "ori50"
		#infilexml.opt		<- "ori60"
		infilexml.opt		<- "ori70"
		argv				<<- hivc.cmd.beast.poolrunxml(indir, infile, insignat, indircov, infilecov, infiletree, infilexml, outsignat, pool.ntip, infilexml.opt=infilexml.opt, infilexml.template=infilexml.template, opt.brl=opt.brl, thresh.brl=thresh.brl, thresh.bs=thresh.bs, resume=resume, verbose=1)
		argv				<<- unlist(strsplit(argv,' '))		
		hivc.prog.BEAST2.generate.xml()
	}
	if(0)		#run BEAST2 BDSKYline
	{
		indir				<- paste(DATA,"tmp",sep='/')
		infile				<- "ATHENA_2013_03_NoDRAll+LANL_Sequences_bdsky_seroneg-130-1_standard_standard"
		infile				<- "ATHENA_2013_03_NoDRAll+LANL_Sequences_bdsky_seroneg-130-1_standard_fxGpInvS4"
		infile				<- "ATHENA_2013_03_NoDRAll+LANL_Sequences_bdsky_seroneg-130-1_standard_fxGpInvS6"
		infile				<- "ATHENA_2013_03_NoDRAll+LANL_Sequences_bdsky_seroneg-130-1_standard_fxGpInvS8"
		infile				<- "ATHENA_2013_03_NoDRAll+LANL_Sequences_bdsky_seroneg-130-1_standard_fxGpInvS10"
		infile				<- "ATHENA_2013_03_NoDRAll+LANL_Sequences_bdsky_seroneg-130-1_standard_fxS4"
		infile				<- "ATHENA_2013_03_NoDRAll+LANL_Sequences_bdsky_seroneg-130-1_standard_fxS5"
		infile				<- "ATHENA_2013_03_NoDRAll+LANL_Sequences_bdsky_seroneg-130-1_standard_fxS8"
		insignat			<- "Tue_Aug_26_09/13/47_2013"
		hpc.ncpu			<- 1

		cmd					<- hivc.cmd.beast2.runxml(indir, infile, insignat, hpc.ncpu=hpc.ncpu, prog.opt.Xmx="1200m")
		cmd					<- cmd.hpcwrapper(cmd, hpc.q="pqeph", hpc.nproc=hpc.ncpu, hpc.walltime=72, hpc.mem="1200mb")
		
		cat(cmd)
		stop()
		signat		<- paste(strsplit(date(),split=' ')[[1]],collapse='_',sep='')
		outdir		<- paste(DATA,"tmp",sep='/')
		outfile		<- paste("b2",signat,sep='.')					
		cmd.hpccaller(outdir, outfile, cmd)
	}
	
}

hivc.proj.pipeline<- function()
{
	stop()
	dir.name<- DATA		 
	signat.in	<- "Wed_May__1_17/08/15_2013"
	signat.out	<- "Wed_May__1_17/08/15_2013"		
	cmd			<- ''

	if(0)	#align sequences in fasta file with Clustalo
	{
		indir	<- paste(dir.name,"tmp",sep='/')
		outdir	<- paste(dir.name,"tmp",sep='/')
		infile	<- "ATHENA_2013_03_FirstAliSequences_HXB2PROTRT_Wed_May__1_17:08:15_2013.fasta"
		infile	<- "ATHENA_2013_03_All+LANL_Sequences_Sat_Jun_16_17:23:46_2013.fasta"
		cmd		<- hivc.cmd.clustalo(indir, infile, signat='', outdir=outdir)
		cmd		<- cmd.hpcwrapper(cmd, hpc.q="pqeph", hpc.nproc=1)
	}	
	if(0)	#extract first sequences for each patient as available from ATHENA data set
	{
		indir		<- paste(dir.name,"tmp",sep='/')
		infile		<- "ATHENA_2013_03_SeqMaster.R"		
		outdir		<- paste(dir.name,"tmp",sep='/')
		cmd			<- paste(cmd,hivc.cmd.get.firstseq(indir, infile, signat.in, signat.out, outdir=outdir),sep='')
	}
	if(0)	#compute raw pairwise genetic distances accounting correctly for ambiguous IUPAC nucleotides 
	{				
		gd.max		<- 0.045
		indir		<- paste(DATA,"tmp",sep='/')
		#infile		<- "ATHENA_2013_03_FirstCurSequences_PROTRT"
		#insignat	<- "Sat_May_11_14/23/46_2013"				
		infile		<- "ATHENA_2013_03_NoDRAll+LANL_Sequences"
		insignat	<- "Thu_Aug_01_17/05/23_2013"
		cmd			<- hivc.cmd.get.geneticdist(indir, infile, insignat, gd.max, outdir=indir)
		
		outfile		<- paste("pipeline",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='.')					
		cmd			<- cmd.hpcwrapper(cmd, hpc.nproc= 1, hpc.q="pqeph")
		cat(cmd)
		cmd.hpccaller(outdir, outfile, cmd)
		quit("no")
	}	
	
	signat	<- paste(strsplit(date(),split=' ')[[1]],collapse='_',sep='')
	outdir	<- paste(dir.name,"tmp",sep='/')
	outfile	<- paste("pipeline",signat,sep='.')					
	lapply(cmd, function(x)
			{				
				#x<- cmd.hpcwrapper(x, hpc.q="pqeph")
				cat(x)
				cmd.hpccaller(outdir, outfile, x)
			})							
}


project.hivc.test<- function()
{
	require(ape)
	if(1)
	{
		x<- as.DNAbin( c("a","c","g","t","m","r","w","s","y","k","v","h","d","b","n") )
		x<- as.DNAbin( as.matrix(c("a","c","g","t","m","r","w","s","y","k","v","h","d","b","n"),1,15) )
		print(x)
		.C("hivc_printdna", x, length(x) ) 
	}
	if(0)
	{
		x<- as.DNAbin( c("a","c","g","t","m","r","w","s","y","k","v","h","d","b","n") )
		y<- as.DNAbin( c("s","c","g","t","m","r","w","s","y","k","v","h","d","b","n") )
		print(as.character(x))
		print(as.character(y))
		tmp<- 0
		gd<- .C("hivc_dist_ambiguous_dna", x, y, length(x), tmp )[[4]]
		print(gd)
	}
	quit("no")
}

