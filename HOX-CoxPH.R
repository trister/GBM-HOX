### Andrew Trister
### Sage Bionetworks
### Seattle, WA
### 20121203

### Build a model of survival based on HOX genes


#first load the necessary libraries
library(synapseClient)
library(survival)
library(ggplot2)
library(affy)
library(org.Hs.eg.db)


#library(probemapper)
#setwd("./trister/")


#here we will load the TCGA data from Synapse
synapseLogin()



#'syn313583', 'syn372761' - these are the SNM normalized data
dataReturn <- loadEntity('syn313583')
eset <- exprs(dataReturn$objects$eset)

metadataLoad <- loadEntity('syn673127') #all of the coherent metadata
metadata <- metadataLoad$objects$metadata #extract the R object
metadataAll <- metadata


# get the intersection of the HGU133A files that also have metadata
rows133a <- unlist(lapply(colnames(eset),function(x){ 
  return(grep(x,metadata[,1]))
}))

metadata <- metadata[rows133a,]
rownames(metadata) <- metadata[,1]

inCommon <- intersect(colnames(eset), rownames(metadata))
eset <- eset[,inCommon]
metadata <- metadata[inCommon,]





#Now let's get an object to have geneIDs from Entrez
x <- org.Hs.egSYMBOL
# Get the gene symbol that are mapped to an entrez gene identifiers
mapped_genes <- mappedkeys(x)
# Convert to a list
xx <- as.list(x[mapped_genes])
if(length(xx) > 0) {
  # Get the SYMBOL for the first five genes
  xx[1:5]
  # Get the first one
  xx[[1]]
}

HOXA3 <- paste(names(xx[grep("^HOXA3",xx)]),"_eg",sep="")
MEOX1 <- paste(names(xx[grep("^MEOX1",xx)]),"_eg",sep="")
LRP6 <- paste(names(xx[grep("^LRP6$",xx)]),"_eg",sep="")




## let's make a survival object
metadata$vital_status[metadata$vital_status=="[Not Available]"] <- NA
gbmPat <- c()
gbmPat$survTime <- as.numeric(metadata$days_to_death)
gbmPat$surv <- ifelse(metadata$vital_status=="DECEASED", 1,0)
gbmPat$survTime[ which(metadata$vital_status=="LIVING")] <- as.numeric(metadata$days_to_last_followup)[ which(metadata$vital_status=="LIVING")]

tmpSurv <- Surv(gbmPat$survTime,gbmPat$surv)






#let's look at all HOX genes
temp <- grep("^HOX",xx)
HOX <- unlist(lapply(temp,function(y){
  return(paste(names(xx[y]),"_eg",sep=""))
}))
HOX <- c(HOX,LRP6,MEOX1)
HOX <- intersect(HOX,rownames(eset))

HOXnames <- unlist(lapply(HOX,function(y){
  return(xx[grep(paste("^",strsplit(y,"_eg")[[1]],"$",sep=""),names(xx))])
}))


cox.output <- lapply(HOX,function(y){
  return(summary(coxph(tmpSurv~eset[y,])))
})

pvals <- unlist(lapply(cox.output,function(y){
  return(y$waldtest[3])
}))


pvals.adjust <- p.adjust(pvals,method="fdr")

print(HOXnames)
print(pvals)
print(pvals.adjust)



HOX2 <- c(HOXA3,MEOX1,LRP6)

cox.output <- lapply(HOX2,function(y){
  return(summary(coxph(tmpSurv~eset[y,])))
})

pvals <- unlist(lapply(cox.output,function(y){
  return(y$waldtest[3])
}))


pvals.adjust <- p.adjust(pvals,method="fdr")

print(pvals)
print(pvals.adjust)





#now look only at the patients that are coded as having had radiation therapy
metadata.culled <- metadata[which(metadata[,"radiation_therapy"]=="YES"),]
#metadata.culled <- metadata[which(!is.na(metadata[,"radiation_type"])),]

eset.culled <- eset[,metadata.culled[,1]]



## let's make a survival object

gbmPat.culled <- c()
gbmPat.culled$survTime <- as.numeric(metadata.culled$days_to_death)
gbmPat.culled$surv <- ifelse(metadata.culled$vital_status=="DECEASED", 1,0)
gbmPat.culled$survTime[ which(metadata.culled$vital_status=="LIVING")] <- as.numeric(metadata.culled$days_to_last_followup)[ which(metadata.culled$vital_status=="LIVING")]

tmpSurv.culled <- Surv(gbmPat.culled$survTime,gbmPat.culled$surv)

cox.output.rad <- lapply(HOX2,function(y){
  return(summary(coxph(tmpSurv.culled~eset.culled[y,])))
})

pvals.rad <- unlist(lapply(cox.output.rad,function(y){
  return(y$waldtest[3])
}))


pvals.adjust.rad <- p.adjust(pvals.rad,method="fdr")

print(pvals.rad)
print(pvals.adjust.rad)

plot(survfit(coxph(tmpSurv.culled~eset.culled[LRP6,]), ))
