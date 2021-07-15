library(tidyverse)
library(dplyr)
library(cowplot)

###########Information###########
diamond_id <-"207254"
dates <- c("20190703",
           "20190705",
           "20190708",
           "20190710",
           "20190712",
           "20190715",
           "20190719",
           "20190724",
           "20190726")
samflags <- "16"
specimen_type <- "plasma"
prefix<-"SLX-19850.DIAMOND.CAMRAD."
postfix<-".s_3.r_1_val_1_GRCh38_bismark_bt2_pe.deduplicated.bam" 
path <- "../bams/"

##for test purposes to limit reads, take out below during loading as well if not needed!!
#number_of_reads<-"100000"

current_run<-paste(gsub(":","",gsub(" ","_",Sys.time())),"read_length_analysis",sep="_")
system(paste("mkdir ",current_run,sep=""))

###########Functions###########
read_sam <- function(sam_file){
  #read sam file
  sam <- read_tsv(sam_file, col_names = F, col_types = cols(
    X1 = col_character(),
    X2 = col_double(),
    X3 = col_character(),
    X4 = col_double(),
    X5 = col_double(),
    X6 = col_character(),
    X7 = col_character(),
    X8 = col_double(),
    X9 = col_double(),
    X10 = col_character(),
    X11 = col_character(),
    X12 = col_character(),
    X13 = col_character(),
    X14 = col_character(),
    X15 = col_character(),
    X16 = col_character()
  ))
  
  if(nrow(sam) == 0){
    stop("Sam file is empty.")
  }
  
  colnames(sam) <- c("QNAME", "FLAG", "RNAME", "POS", "MAPQ", "CIGAR", "RNEXT", "PNEXT",
                     "TLEN", "SEQ", "QUAL", "EDIT", "MISMAT", "METHYL", "RCON", "GCON")
  return(sam)
}

process_sam <- function(processed_sam, diamond_id, date, specimen_type){
  print(paste("Starting to process",diamond_id,date,specimen_type))
  
  pdf(paste("./",current_run,"/",paste("CAMRAD",diamond_id,date,"processesed.pdf",sep="_"),sep=""))
  print(paste(paste("CAMRAD",diamond_id,date,"processesed.pdf",sep="_"),"created"))
  
  #creates insert length distribution 
  distribution <- as.data.frame(table(filter(processed_sam, TLEN > 0)$TLEN))
  names(distribution)<-c("TLEN","Freq")
  #make TLEN numbers rather than factors
  distribution$TLEN <- as.numeric(as.character(distribution$TLEN))
  total_reads<-sum(distribution$Freq)
  
  #censor the bin at 126 -- this is a hack and needs to be taken out
  distribution$Freq[which(distribution$TLEN == 126)] <- distribution$Freq[which(distribution$TLEN == 126)-1]
  
  #collect all reads with TLEN 126 into CSV file
  write.csv(processed_sam[processed_sam$TLEN == 126,],paste("./",current_run,"/",paste("CAMRAD",diamond_id,date,"_126_reads.csv",sep="_"),sep=""))
  
  hist(processed_sam$TLEN, xlim = c(0,max(distribution$TLEN)),breaks=max(distribution$TLEN),
       main=paste("Read Length Histogram for",diamond_id,date,specimen_type),
       xlab="Read length",
       ylab="Frequency")
  
  hist(processed_sam$TLEN, xlim = c(100,200),breaks=max(distribution$TLEN),
       main=paste("Focussed Read Length Histogram for",diamond_id,date,specimen_type),
       xlab="Read length",
       ylab="Frequency")
  
  print(paste("Max Peak at: ",distribution$TLEN[which(distribution$Freq == max(distribution$Freq))]))
  
  #adds column with difference in insert sizes to find largest gap
  #distribution$gap <- append(diff(as.numeric(as.vector(distribution$TLEN))), 1)
  distribution$gap <- append(diff(distribution$TLEN), 1)
  
  #block 1 to only contain the first valley not present as true low
  block1 <- subset(distribution, TLEN > 150 & TLEN < 300)
  #find row with the biggest gap, get TLEN from this but add half the read length gap to land in the middle
  cut1<-block1[which.max(block1$gap),]$TLEN+floor(max(block1$gap)/2)
  
  #block 2 to only contain the second valley and find any read length not present as true low
  block2 <- subset(distribution, TLEN > 300)
  #find row with the biggest gap, get TLEN from this but add half the read length gap to land in the middle
  cut2<-block2[which.max(block2$gap),]$TLEN+floor(max(block2$gap)/2)
  
  #if no read length not present then select the read length with the fewest reads instead
  #note that these cutoffs depend on the arbitrary search intervals in the rows above!
  if(cut1==151){
    cut1<-block1[which.min(block1$Freq),]$TLEN
  }
  
  if(cut2==301){
    cut2<-block2[which.min(block2$Freq),]$TLEN
  }
  
  #define peak regions based on cut points
  
  peak1<-subset(distribution, TLEN < cut1)
  peak2<-subset(distribution, TLEN > cut1 & TLEN < cut2)
  peak3<-subset(distribution, TLEN > cut2)
  
  results <- data.frame(c(diamond_id),
                        c(date),
                        c(specimen_type),
                        c(cut1),
                        c(cut2),
                        c(median(peak1$TLEN)),
                        c(round(sum(peak1$Freq)/total_reads,digits=3)),
                        c(median(peak2$TLEN)),
                        c(round(sum(peak2$Freq)/total_reads,digits=3)),
                        c(median(peak3$TLEN)),
                        c(round(sum(peak3$Freq)/total_reads,digits=3))
  )
  
  names(results)<-c("Diamond ID",
                    "Date",
                    "Specimen Type",
                    "Valley 1",
                    "Valley 2", 
                    "Median Size Peak 1", 
                    "Proportion of Reads in Peak 1", 
                    "Median Size Peak 2", 
                    "Proportion of Reads in Peak 2", 
                    "Median Size Peak 3", 
                    "Proportion of Reads in Peak 3")
  
  dev.off()
  return(results)
}


###########Main Programme###########

#generate SAM from BAM
print("Starting...")
system(paste("samtools view -F ",samflags,"  ",path, prefix, diamond_id,".",dates[1],".",specimen_type,postfix," > tmp_",current_run,".sam",sep=""))
#load SAM
sam <- read_sam(paste("tmp_",current_run,".sam",sep=""))
#process SAM
results<-process_sam(sam,diamond_id,dates[1],specimen_type)
print("Sam 1 finished")

for(i in 2:length(dates)){
  #generate SAM from BAM
  
  #if limiting reads use the version with "head"
  #system(paste("samtools view -F ",samflags,"  ",path, prefix, diamond_id,".",dates[i],".",specimen_type,postfix," | head -n ",number_of_reads," > tmp.sam",sep=""))
  system(paste("samtools view -F ",samflags,"  ",path, prefix, diamond_id,".",dates[i],".",specimen_type,postfix," > tmp_",current_run,".sam",sep=""))
  
  #load SAM
  sam <- read_sam(paste("tmp_",current_run,".sam",sep=""))
  #process SAM
  results <- rbind(results,process_sam(sam,diamond_id,dates[i],specimen_type))
  print(paste("Sam",i,"finished"))
}

results$Date<-as.character(results$Date)
results$Date<-as.Date(results$Date,format="%Y%m%d")
results$daysontreat <- difftime(results$Date,results$Date[1], units = "days")


plot1<-ggplot(data=results, aes(x=daysontreat,y=`Median Size Peak 1`,group=1)) +
  labs(title="Median Fragment Length Peak 1", x ="Days on Treatment", y = "Median Fragment Length")+
  geom_line(col="green")+
  geom_point()


plot2 <-ggplot(data=results, aes(x=daysontreat,y=`Median Size Peak 2`,group=1)) +
  labs(title="Median Fragment Length Peak 2", x ="Days on Treatment", y = "Median Fragment Length")+
  geom_line(col="green")+
  geom_point()

plot3 <-ggplot(data=results, aes(x=daysontreat,y=`Median Size Peak 3`,group=1)) +
  labs(title="Median Fragment Length Peak 3", x ="Days on Treatment", y = "Median Fragment Length")+
  geom_line(col="green")+
  geom_point()

plot4 <- ggplot(data=results, aes(x=daysontreat,y=`Proportion of Reads in Peak 1`,group=1)) +
  labs(title="Proportion of reads in Peak 1", x ="Days on Treatment", y = "Median Fragment Length")+
  geom_line(col="red")+
  geom_point()

plot5 <- ggplot(data=results, aes(x=daysontreat,y=`Proportion of Reads in Peak 2`,group=1)) +
  labs(title="Proportion of reads in Peak 1", x ="Days on Treatment", y = "Median Fragment Length")+
  geom_line(col="red")+
  geom_point()

plot6 <- ggplot(data=results, aes(x=daysontreat,y=`Proportion of Reads in Peak 3`,group=1)) +
  labs(title="Proportion of reads in Peak 1", x ="Days on Treatment", y = "Median Fragment Length")+
  geom_line(col="red")+
  geom_point()

pdf(paste("./",current_run,"/",paste("CAMRAD",diamond_id,"plots.pdf",sep="_"),sep=""))
plot_grid(plot1, plot2, plot3, plot4, plot5,plot6, labels = "AUTO")
dev.off()

write.csv(results,paste("./",current_run,"/",paste("CAMRAD",diamond_id,"read.lengths.csv",sep="_"),sep=""), row.names = FALSE)
system(paste("rm tmp_",current_run,".sam",sep=""))
print("All done and tidy")
