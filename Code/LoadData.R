#############################################################################
################################# LOAD DATA #################################
#############################################################################
LoadPeaks<-function(){
  # Peaks data first!
  peaks<-read.delim(paste0(dir,"/Data/peak_data.txt"),header = F)
  # It's two lines, first describes properties of , 
  # the second the DNA sequence of this specific
  peaks<-data.frame(desc=peaks[seq(1,nrow(peaks),by=2),1],
                    dna=peaks[seq(2,nrow(peaks),by=2),1])
  # House-cleaning
  peaks$desc<-str_replace_all(str_replace_all(peaks[,1],">>",""),"_\\+","")
  # 'Random' is a sequence that was in an unfinished state, but whose chromosome is known
  peaks<-peaks[!grepl(x = peaks$desc,"random"),]; print("Removing random sequences")
  # Chromosome 'Un' are contigs that cannot be confidently placed on a specific hormone
  peaks<-peaks[!grepl(x = peaks$desc,"chrUn"),]; print("Removing chrUn sequences")
  # Break down the description into its constituent characteristics
  tmp<-sapply(peaks$desc,FUN = function(p) str_split(p,"_")[[1]],simplify = F)
  # Check for anomalous entries
  if(length(tmp[unlist(sapply(tmp,length))>6])>0 | length(tmp[unlist(sapply(tmp,length))<6])>0) print("Some incorrect sequence descriptions were found")
  # Get this variable into data.frame format
  tmp%<>%as.data.frame()%>%as.matrix()%>%t()%>%as.data.frame()
  # Change the column names
  colnames(tmp)<-c("GenomeBuild","Chromosome","StartP","EndP","Other","FoldEnrich")
  # Add this to the peaks data frame
  peaks%<>%cbind(tmp)
  # Neaten up the data
  peaks$Other<-str_replace_all(peaks$Other,"reg","")
  peaks$Chromosome<-str_replace_all(peaks$Chromosome,"chr","")
  # Split the 'Other' column and extract the two sections
  tmp<-str_split(peaks$Other,"\\.")%>%as.data.frame()%>%as.matrix()%>%t()%>%as.data.frame()
  colnames(tmp)<-c("Namer","pScore")
  # Bind to the peaks data frame
  peaks%<>%cbind(tmp); rm(tmp)
  # Housecleaning
  peaks$FoldEnrich%<>%as.numeric(); peaks$StartP%<>%as.integer(); peaks$EndP%<>%as.integer(); peaks$Namer%<>%as.integer()
  # Add a variable for the length of the DNA section
  peaks$lengthDNA<-abs(peaks$EndP-peaks$StartP)
  # Check for NA values
  apply(peaks,2,function(i) any(is.na(i)))
  # Cleaning!
  peaks[,c("Other", "desc")]<-NULL
  # Modify the p-value of the signal to control
  peaks$pScore<-str_replace_all(peaks$pScore,"p","")%>%as.integer()
  
  return(peaks)
}

LoadShuffle<-function(){
  # Raw data
  shuffle<-read.delim(paste0(dir,"/Data/shuffled_data.txt"),header = F)
  # It's two lines, first describes properties of , 
  # the second the DNA sequence of this specific
  shuffle<-data.frame(desc=shuffle[seq(1,nrow(shuffle),by=2),1],
                    dna=shuffle[seq(2,nrow(shuffle),by=2),1])
  # House-cleaning
  shuffle$desc<-str_replace_all(str_replace_all(shuffle[,1],">>",""),"_\\+","")
  # Split up the description to extract the important properties
  tmp<-str_split(shuffle$desc,"_")%>%as.data.frame()%>%as.matrix()%>%t()%>%as.data.frame()
  # Rename
  colnames(tmp)<-c("GenomeBuild","StartP","EndP")
  # Add to shuffle data frame
  shuffle%<>%cbind(tmp)
  # Class modifications
  shuffle$StartP%<>%as.integer(); shuffle$EndP%<>%as.integer()
  # House-cleaning!
  rm(tmp); shuffle$desc<-NULL
  # Add length of the DNA section:
  shuffle$lengthDNA<-abs(shuffle$StartP-shuffle$EndP)
  
  return(shuffle)
}






