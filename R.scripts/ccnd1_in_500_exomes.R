##Load libraries
library(RMySQL)
library(deepSNV)

## Get list of 463 samples to test (both NORMAL and TUMOUR)

## Function to fetch data
sqldo = function(this.sql){
  ## Load MySQL driver and connect to db
  drv = dbDriver("MySQL")
  con = dbConnect(drv,dbname="XYREMSEQ",host="home",username="homeweb",password="opensesame")
  res=dbSendQuery(con,this.sql) # Execute the query
  resdf=fetch(res,-1) # Get query results as a dataframe
  dbClearResult(res) # Clear query
  dbDisconnect(con) # Disconnect
  return(resdf)
}

## Grab project info/translocation/mlpa tables to subdivide data
sql = "SELECT * FROM PROJECT_INFO
LEFT JOIN TRANSLOCATIONS ON
PROJECT_INFO.Patient_sample = TRANSLOCATIONS.patient_sample
LEFT JOIN MLPA ON
PROJECT_INFO.Patient_sample = MLPA.patient_sample"
pinfo = sqldo(sql)

## Get the 926 normal/tumour pairs in the 500 exome project
the926=pinfo[pinfo[,"500_exomes"]=="Yes",]
the926normal=pinfo[pinfo[,"500_exomes"]=="Yes" & pinfo[,"SampleType"]=="PB",]
the926tumour=pinfo[pinfo[,"500_exomes"]=="Yes" & pinfo[,"SampleType"]=="BM",]

## Loop through the926 to fetch values and store in an object
normalresults=as.data.frame(matrix(NA,nrow=nrow(the926tumour),ncol=4))
tumourresults=as.data.frame(matrix(NA,nrow=nrow(the926tumour),ncol=4))

for(i in 1:nrow(the926tumour)){
  theCounts = function(samplename){
    message(paste("PROCESSING:",samplename,i))
    filepath=paste0("/mnt/exomes/",samplename,"/realignedBAM/dedup.realigned.bam")
    counts = bam2R(file=filepath,chr="11",start=69462910,stop=69462910)
    gcount = counts[,"G"]+counts[,"g"]
    acount = counts[,"A"]+counts[,"a"]
    depth = sum(counts[,colnames(counts)%in%c("A","C","G","T","a","c","g","t")])
    message(paste("SUCCESS:",samplename,i))
    return(c(acount,gcount,depth))
  }
  samplename = the926normal[i,"Patient_sample"]
  normalresults[i,1] = samplename
  normalresults[i,2:4] = theCounts(samplename)

  samplename = the926tumour[i,"Patient_sample"]
  tumourresults[i,1] = samplename
  tumourresults[i,2:4] = theCounts(samplename)
}

## Now merge tumour/normal sesults into one object 
mergedResults = cbind(normalresults,tumourresults)

## Write output to file
write.table(mergedResults,file="CCDN1_snp_463exomes.txt",quote=FALSE,sep="\t",row.names=FALSE,col.names=TRUE)



