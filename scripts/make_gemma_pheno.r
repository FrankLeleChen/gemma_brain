#Read family file data to establish subject ordering  
#CHANGE grac to studyname

args <- commandArgs(trailingOnly = TRUE)
subjects <- args[1]
inclusions <- args[2]
phenotypes <- args[3]
phenolist <- args[4]
covars <- args[5]
covlist <- args[6] 
studycov <- args[7]
outfilename <- args[8]

#String splitting functions
 unlist_split <- function(x,  ...)
 {
     toret <- unlist(strsplit(x, ...) )
     return(t(toret))
 }
 unlist_split2 <- function(x, colnum ,...)
 {
     toret <- unlist(strsplit(x, ...) )[colnum]
     return(toret)
 }

#disabled: Median imputation function (missing covariates will be imputed to median value)
 # imputemed=function(x){
   # x<-as.numeric(as.character(x)) #first convert each column into numeric if it is from factor
   # x[is.na(x)] =median(x, na.rm=TRUE) #convert the item with NA to median value from the column
   # x #display the column
 # }
#Make list of phenotypes R readable
 phenolist_r <- c(t(unlist_split(phenolist,","))) 

#Make the list of covariates R readable
 covlist_r <- c(t(unlist_split(covlist,","))) 

#Load subject list from .fam file
 subjs <- read.table(subjects,header=F,stringsAsFactors=F,na.strings=c("NA","-9"))
 names(subjs)[1:2] <- c("FID","IID")
 subjs$order <- 1:dim(subjs)[1] #Note .fam file order for data ordering purposes, very important, matches GRM. No IDs in exported file!!


#Load list of people to include in analysis
 if (inclusions != "xxxx")
 {
  include <- read.table(inclusions,header=T,stringsAsFactors=F,na.strings=c("NA","-9"))
  names(include) <- c("FID","IID")
  include$keep_sub <- 1 

  subjs <- merge(subjs,include,by=c("FID","IID"),all.x=T)
 } else { subjs$keep_sub <- 1 }
 
 
#Load in phenotype data. Blank space delimited. Must have header. Header must have FID, IID as first two columns. 
 phenodata <- read.table(phenotypes,header=T,stringsAsFactors=F,na.strings=c("NA","-9"))

#Join phenos to subject list
 dat0  <- merge(subjs,phenodata,all.x=TRUE,by=c("FID","IID"),suffixes=c("","_phen"))
 
#Load in other covariates (blank space delimited file with header,s must have FID,IID)
 if(covars != "xxxx")
 {
  covdata <- read.table(covars,header=T,stringsAsFactors=F)
  
  dat <- merge(dat0,covdata,all.x=T,by=c("FID","IID"),suffixes=c("","_cov"))
  
 } else { dat <- dat0 }

#Sort to original order
 dat <- dat[order(dat$order),]
 
#For all subjects who should be out of the analysis, set all of their phenotypes to NA
if( length(which(dat$keep_sub == 1)) != dim(dat)[1])
{
 dat[-which(dat$keep_sub == 1),phenolist_r] <- NA
}

#Need intercept for GEMMA
 dat$intercept <- 1 

#Make a study covariate
 if (studycov == TRUE)
 {
  
  dat$studycov <- sapply(dat$FID,unlist_split2,colnum=3,split="_") #In Ricopili format, study site is the third item in the FID
  removetemp <- which(is.na(dat[,phenolist_r[1]]))

 
  dattemp <- dat[-removetemp,] #[,] #Will only make a study covariate for subjects with a phenotype value for the first phenotype in the data! 
  print(dim(dattemp))
  recoded_site <-  model.matrix(~studycov, data=dattemp, na.action='na.fail')[,-1] #Make the model matrix, fail if any subject is missing a covariate

  colnames(recoded_site) <- paste(colnames(recoded_site),sep="_") #Rename if you want to here..
  indicator_cols <- colnames(recoded_site) #List the names of these indicator columns, used to subset covariate sheet
  recoded_site <- data.frame(recoded_site)
  recoded_site$FID <- dattemp$FID
  recoded_site$IID <- dattemp$IID

  dat2 <- merge(dat,recoded_site,by=c("FID","IID"),all.x=TRUE)
 } else {dat2 <- dat}
 
#Get a dataframe of covariates
 print(names(dat2))
#Make sure the order is correct
 dat2 <- dat2[order(dat2$order),]
 
 #If there are covariates, keep them
 if (covlist_r[1] != "xxxx")
 {
  datexp2 <- subset(dat2, select=c("intercept",covlist_r,indicator_cols))
 } else { datexp2 <- subset(dat2, select=c("intercept",indicator_cols)) }
 
 #datexp2 <- data.frame(apply(datexp2a,2,imputemed))
 
 # # #Look at phenotype available subjects. If in this set a study covariate is all uniform values, remove it
 # # phenosubs <- which(!is.na(dat[,phenolist_r[1]]))
 
 # # for (scov in colnames(recoded_site))
 # # {
  # # if( all(dat[phenosubs,scov] == 0))
  # # {
   # # print(scov)
   # # print("must be removed!")
   # # datexp2 <- datexp2[,names(datexp2a) != scov] #remove it
  # # }
 # # }
 
#Check for remaining NA covariates - they are not allowed!!

 write.table(dat2[,phenolist_r], paste(outfilename,".pheno",sep=""),row.names=F,col.names=F,quote=F)
 # print(dim(datexp2))
 print(head(as.matrix(datexp2)))
 write.table(datexp2,paste(outfilename,".covar",sep=""),row.names=F,col.names=F,quote=F)

  #NOTE: SEGFAULT WILL HAPPEN IF DIMENSION OF DATA DOES NOT MATCH DIMENSION OF GRM!
  #NOTE: SEGFAULT WILL HAPPEN IF THERE ARE NA VALUES IN COVARIATE SHEET!
  #NOTE: POTENTIAL SEGFAULT DUE TO COLLINEARITY OF DUMMY CODES, E.G. Covariate for sex, covariate for study, but all subjects in study are female
 #Export individual phenotypes to table
  for (i in phenolist_r)
  {
    write.table(dat2[,i], paste(outfilename,".",i,".pheno",sep=""),row.names=F,col.names=F,quote=F)
  }