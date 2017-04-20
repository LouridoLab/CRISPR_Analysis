
#################### USER-DEFINED VARIABLES ####################

## Define the local directory. If the local directory is on the desktop, replace "CRISPR Analysis" with the name of your directory. 
setwd("~/Desktop/CRISPR Analysis/") # Modify the text in quotations

## Define the directory that contains the FASTQ files
dir.name <- "FASTQ" # Modify the text in quotations

## Define the file containing the sgRNA sequences
lib.name <-"sgRNAs_nr_empty" # Modify the text in quotations

## Indicate the file names for the different samples exclude the *.fastq extension
L <- "B4_L" # Modify the text in quotations

P3 <- "B4_P3" # Modify the text in quotations

################################################################

## Retrieve the contents of the local directory.
file.contents <- list.files()

## Check for the presence of the directory containing the FASTQ files
try(if(!(dir.name %in% file.contents)) stop("ERROR: Directory not found!"))

## Check for the presence of the file containing the sgRNA sequences
try(if(!(lib.name %in% file.contents)) stop("ERROR: sgRNA file not found!"))

## Retrieve the contents of the directory containing the FASTQ files
file.contents <- list.files(dir.name)

## Check for the presence of the first sample
try(if(!(paste0(L,".fastq") %in% file.contents)) stop("ERROR: Sample file not found!")) 

## Check for the presence of the second sample
try(if(!(paste0(P3,".fastq") %in% file.contents)) stop("ERROR: Sample file not found!")) 


## Ensure that Countess is executable
system(command = "chmod 770 Countess")

## Run Countess on the FASTQ files found in the defined directory

com <- paste("./Countess -l ",lib.name, " -d ", dir.name)
system(command = com)

## Aggregate samples into data frame "DF.sg"
DF.sg <- data.frame()
count.dir.name <- paste(dir.name, "_Counted", sep = "")
file_list <- list.files(path = count.dir.name)

for (file in file_list){
  file.path <- paste(count.dir.name,"/",file, sep = "")
  # if the merged dataset doesn't exist, create it
  if (nrow(DF.sg)==0){
    DF.sg <- read.table(file.path, row.names=1)
    colnames(DF.sg) <- c(file)
  } else {
    temp_dataset <-read.table(file.path, row.names=1)
    colnames(temp_dataset) <- c(file)
    DF.sg<-merge(DF.sg, temp_dataset, by="row.names", all=T)
    rownames(DF.sg)<-DF.sg[,1]
    DF.sg[,1]<-NULL
    rm(temp_dataset)
  }
} 

## Match normalize
for (i in 1:ncol(DF.sg))
{
  x.sum <- sum(DF.sg[,i])
  DF.sg[,i] <- c(DF.sg[,i]/x.sum)
}

## Replace zeros
for (i in 1:ncol(DF.sg))
{
  x <- DF.sg[,i]
  x[x == 0] <- NA
  x[is.na(x)] <- .9*min(x, na.rm=T)
  DF.sg[,i] <- x
  rm(x)
}

## Log2 transform
DF.sg <- log2(DF.sg)

## Define the quantile to exclude
Q <- 0.05 #quantile to be used in pruning

## Remove the sgRNAs in the lowest quantile of the library
DF.sg <- DF.sg[DF.sg[,L] > quantile(DF.sg[,L], Q, na.rm = T),]

##Calculate fold change
DF.sg$FC<-c(DF.sg[,P3]-DF.sg[,L])

## Extract the gene IDs from the sgRNA IDs
DF.sg$ID <- sub(".*(TGGT1_[0-9]+[A-Z]?).*","\\1",rownames(DF.sg))

## Create a new data frame called "DF.genes" with the gene IDs
DF.genes <- data.frame(ID=unique(DF.sg$ID))

## Define the function to calculate standard error of the mean (SEM)
sem <- function(x) sd(x, na.rm=T)/sqrt(sum(!is.na(x)))

## Calculate the phenotype score and its SEM.
## n equals the number of sgRNAs used in the calculation.
for (i in 1:length(DF.genes$ID))
{
  gene <- as.character(DF.genes$ID[i])
  scores <- subset(DF.sg, DF.sg$ID == gene)[,"FC"]
  scores <- scores[!is.na(scores)]
  if (length(scores) > 5)
  {
    scores <- scores[order(-scores)]
    DF.genes$phenotype[i] <- mean(scores[1:5])
    DF.genes$sem[i] <- sem(scores[1:5])
    DF.genes$n[i] <- 5
  }
  else if (length(scores) >= 2)
  {
    DF.genes$phenotype[i] <- mean(scores)
    DF.genes$sem[i] <- sem(scores)
    DF.genes$n[i] <- length(scores)
  }
  else
  {
    DF.genes$phenotype[i] <- NA
    DF.genes$sem[i] <- NA
    DF.genes$n[i] <- length(scores)
  }
}

## Rank-order the genes and plot by phenotype
DF.genes <- DF.genes[order(DF.genes$phenotype),]
DF.genes$pos <- c(1:length(DF.genes$ID))

## Define parameters for graphing and plot
SE.up = DF.genes$phenotype + DF.genes$sem
SE.dn = DF.genes$phenotype - DF.genes$sem
plot(DF.genes$phenotype~DF.genes$pos,pch=20, cex=0.1, las=1, ylab="phenotype", xlab="rank-ordered genes")
abline(0,0)
arrows(DF.genes$pos,SE.dn,DF.genes$pos,SE.up,code=3,length=0,angle=90,col='gray85')
points(DF.genes$phenotype~DF.genes$pos, pch=20, cex=.5)

Neg.Ctrl.Genes <- c("TGGT1_292055","TGGT1_295760","TGGT1_217600","TGGT1_240390","TGGT1_232130","TGGT1_235380","TGGT1_225490","TGGT1_205380","TGGT1_297880","TGGT1_255190","TGGT1_283510","TGGT1_246940","TGGT1_228400","TGGT1_243400","TGGT1_309870","TGGT1_253770","TGGT1_306620","TGGT1_305860","TGGT1_205480","TGGT1_310160","TGGT1_286590","TGGT1_201840","TGGT1_205250","TGGT1_233030","TGGT1_260820","TGGT1_262730","TGGT1_204130","TGGT1_277080","TGGT1_263180","TGGT1_250880","TGGT1_227620","TGGT1_312480","TGGT1_309590","TGGT1_233460","TGGT1_315760","TGGT1_294690","TGGT1_254555","TGGT1_282055","TGGT1_254470","TGGT1_319560")

Neg.Ctrl.Genes <- subset(DF.genes, DF.genes$ID %in% Neg.Ctrl.Genes)
Ref.mean <- mean(Neg.Ctrl.Genes$phenotype)

points(Neg.Ctrl.Genes$pos, Neg.Ctrl.Genes$phenotype, col="orange", pch=19)

abline(h=Ref.mean, lty=2, lwd=2, col="orange")

legend("bottomright", inset=.08, title="", c("mean phenotype", "known dispensable"), horiz=F, cex=1, bty="n", lty=2, lwd=2, col=c("orange", "transparent"), text.col="transparent")
legend("bottomright", inset=.08, title="", c("mean phenotype", "known dispensable"), horiz=F, cex=1, bty="n", pch=19, col=c("transparent", "orange"))

## Extract the sgRNA values for the control genes
control <- subset(DF.sg, DF.sg$ID %in% Neg.Ctrl.Genes$ID)

## Calculate the p value for each gene
for (i in 1:length(DF.genes$ID))
{
  scores <- DF.sg$FC[DF.sg$ID==as.character(DF.genes$ID[i])]
  scores <- scores[!is.na(scores)]
  if (length(scores) > 2)
  {
    DF.genes$p.value[i] <- wilcox.test(scores,control$FC, alternative = "less")$p.value
  } else {
    DF.genes$p.value[i] <- NA
  }
}

## Correct the p values for multiple hypothesis testing
DF.genes$FDR = p.adjust(DF.genes$p.value, method="fdr")

## Define fitness conferring genes as those with FDR less than or equal to 0.05
for (i in 1:length(DF.genes$ID))
{
  if (!is.na(DF.genes$FDR[i]) & DF.genes$FDR[i] <= 0.05)
  {
    DF.genes$fitness.conferring[i] <- "yes"  
  } else {
    DF.genes$fitness.conferring[i] <- "no"  
  }
}

## Extract the fitness conferring genes
fitness.conferring <- DF.genes[DF.genes$fitness.conferring=="yes",]

## Plot highlighting the fitness-conferring genes
plot(DF.genes$phenotype~DF.genes$pos,pch=20, cex=0.1, las=1, ylab="phenotype", xlab="rank-ordered genes")
abline(0,0)
arrows(DF.genes$pos,SE.dn,DF.genes$pos,SE.up,code=3,length=0,angle=90,col='gray85')
points(DF.genes$phenotype~DF.genes$pos, pch=20, cex=.5)

points(fitness.conferring$pos, fitness.conferring$phenotype, col="red", pch=19)
points(Neg.Ctrl.Genes$pos, Neg.Ctrl.Genes$phenotype, col="orange", pch=19)

abline(h=Ref.mean, lty=2, lwd=2, col="orange")

legend("bottomright", inset=.08, title="", c("mean phenotype", "known dispensable", "fitness-conferring"), horiz=F, cex=1, bty="n", lty=2, lwd=2, col=c("orange", "transparent", "transparent"), text.col="transparent")
legend("bottomright", inset=.08, title="", c("mean phenotype", "known dispensable", "fitness-conferring"), horiz=F, cex=1, bty="n", pch=19, col=c("transparent", "orange", "red"))

## Save data frame
write.csv(DF.genes, file="GenePhenotypes.csv", row.names=F)

