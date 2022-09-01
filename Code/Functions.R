# Convert a string to a one-hot matrix
# Input: 
#     s - character string
#     n - the encodings to be used
# Output: 
#     One-hot array
oneHot <- function(s, n = c("A", "C", "G", "T"), maxL=300) {
  # Construct matrix
  s <- toupper(s)
  # Make sure the output length is that what we need
  if(nchar(s)>maxL) stop("Length of one-hot array must be more than the largest sequence length")
  # Split the sequence
  seq_split <- unlist(strsplit(x = s, split = ""))
  # Create the array skeleton
  seq_mat <- matrix(data = rep(0, maxL * length(n)), nrow = length(n))
  # Dimensionality
  rownames(seq_mat) <- n
  colnames(seq_mat) <- c(seq_split,rep("",maxL-length(seq_split)))
  # Encode
  for (i in n) seq_mat[rownames(seq_mat) == i, colnames(seq_mat) == i] <- 1
  
  return(seq_mat)
}

# Function that converts from a DNA sequence to a one-hot array indicating for the bases
# Input:
#     dna is a character vector that contains the sequence, i.e. ACGCTGCGAAAA...
#     maxL is the largest base to be used (determines the array dimension for CNN)
# Output:
#     a list of length length(dna) of one-hot arrays, each of dimension c(maxL,4) 
#     where 4 represents A,C,G,T.
ConvOneHot<-function(dna,maxL=NULL,lennyD=NULL,twoD=F,aug=F){
  # Extract the maximum DNA sequence length
  maxxer<-max(str_length(dna))
  # Check for NULL maxL values and default to the max string length of dna
  if(is.null(maxL)) {maxL<-maxxer} else if(maxL<maxxer) maxL<-maxxer
  # Create the one-hot arrays, one-by-one
  tmp<-sapply(dna,function(inp) oneHot(inp, maxL=maxL),simplify = F)
  # Condense this into one single array
  # This dimensionality of this single array depends on the model:
  if(twoD){
    # For 2D CNNs, add the length data as an extra layer
    if(is.null(lennyD)) stop("2D CNN model requires the length data input to one-hot encoder")
    # Skeleton
    out<-array(NA,dim=c(length(tmp),4,maxL,2))
    for(i in 1:length(tmp)) {
      # One-hot
      out[i,,,1]<-tmp[[i]]
      # length of DNA sequence
      out[i,,,2]<-tmp[[i]]*lennyD[i]
    }
  } else{
    # For 1D CNNs, either use only one-hot values or the length data
    # Skeleton
    out<-array(NA,dim=c(length(tmp),4,maxL))
    # One-hot
    if(is.null(lennyD)) {for(i in 1:length(tmp)) out[i,,]<-tmp[[i]]
    # length of DNA sequence
    }else for(i in 1:length(tmp)) out[i,,]<-tmp[[i]]*lennyD[i]
  }
  # Do Transcription Factors read both backwards and forwards and yet produce the same protein?
  if(aug)  stop("not ready yet")
  return(out)
}

# Permute DNA bases
PermuteBase<-function(OneHot,ord=1:4){
  # Skeleton
  tmp<-OneHot
  # Re-order
  # For 1D CNN
  if(length(dim(OneHot))==4) {for(i in 1:4) OneHot[,i,,]<-tmp[,ord[i],,]
  # For 2D CNN
  } else if(length(dim(OneHot))==5) {for(i in 1:4) OneHot[,i,,,]<-tmp[,ord[i],,,]
  } else stop("Incorrect dimensions of one-hot array given to permutation function")
  
  return(OneHot)
}

PermutateData<-function(Data,ord=1:4){
  Data$X_Train%<>%PermuteBase(ord)
  Data$X_Test%<>%PermuteBase(ord)
  return(Data)
}

SortSCV<-function(indies,OneH,cv,maxL,SCV=5){
  OtherP<-OtherS<-c()
  for(ii in (1:SCV)[1:SCV!=cv]) {
    OtherP<-c(OtherP,indies$P[[ii]])
    OtherS<-c(OtherS,indies$S[[ii]])  
  }
  # Reduce the shuffle training data to be compatible with the length of the peak training data
  iStrain<-sample(OtherS,length(OtherP))
  iStest<-c(OtherS[!OtherS%in%iStrain],indies$S[[cv]])
  # Calculate the lengths
  lTrainP<-length(OtherP)
  lTrainS<-length(iStrain)
  lTestP<-length(indies$P[[cv]])
  lTestS<-length(iStest)
  # Templates
  X_Train<-array(NA,c(lTrainP+lTrainS,4,maxL))
  X_Test<-array(NA,c(lTestP+lTestS,4,maxL))
  Y_Train<-rep(NA,(lTrainP+lTrainS))
  Y_Test<-rep(NA,(lTestP+lTestS))
  ################ TRAINING DATA ################
  # Peak data:
  X_Train[1:lTrainP,,]<-OneH$peaks[OtherP,,]
  Y_Train<-rep(1,lTrainP)
  # Shuffle data:
  X_Train[(lTrainP+1):(lTrainP+lTrainS),,]<-OneH$shuffle[iStrain,,]
  Y_Train<-c(Y_Train,rep(0,lTrainS))
  ################## TEST DATA ##################
  # Peak data:
  X_Test[1:lTestP,,]<-OneH$peaks[indies$P[[cv]],,]
  Y_Test<-rep(1,lTestP)
  # Shuffle data:
  X_Test[(lTestP+1):(lTestP+lTestS),,]<-OneH$shuffle[iStest,,]
  Y_Test<-c(Y_Test,rep(0,lTestS))
  
  return(list(X_Train=array_reshape(X_Train, c(dim(X_Train)[1], 4, maxL, 1)),
              X_Test=array_reshape(X_Test, c(dim(X_Test)[1], 4, maxL, 1)),
              Y_Train=to_categorical(Y_Train, 2),
              Y_Test=to_categorical(Y_Test, 2)))
  
}
