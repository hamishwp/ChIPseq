dir<-getwd()
# Install & Load all the necessary packages & required functions
source(paste0(dir,"/Code/Packages.R"))
source(paste0(dir,"/Code/LoadData.R"))
source(paste0(dir,"/Code/Functions.R"))
# Load in the ChIP-seq datasets
peaks<-LoadPeaks()
# Load also the shuffle data
shuffle<-LoadShuffle()

################################################################################
######################### SUPPORT VECTOR MACHINE - SVM #########################
################################################################################
# Data wrangling:
peaks$Bound<-"Unbound"; shuffle$Bound="Bound"; shuffle$Chromosome=100
ALL<-rbind(dplyr::select(peaks,c(lengthDNA,StartP,EndP,Chromosome,Bound)),
           dplyr::select(shuffle,c(lengthDNA,StartP,EndP,Chromosome,Bound)))
ALL%<>%group_by(Chromosome)%>%mutate(modStart=StartP/min(StartP),modEnd=EndP/min(EndP))%>%ungroup()
ALL%<>%dplyr::select(-Chromosome)
# Peaks:
process<-ALL%>%preProcess(method=c("range"))
ALL<-predict(process, ALL)
ALL$Bound%<>%factor()
# SVM
# Repeated Cross-Validation
train_control <- trainControl(method="repeatedcv", number=5, repeats=10,
                              classProbs = TRUE, summaryFunction = twoClassSummary)
SVMout<-data.frame()
# Run the SVM!
for(meth in c("svmLinear","svmRadial","svmPoly")){
  # Using the chromosome-based normalisation of position
  svmDNA <- caret::train(Bound ~., 
                         data = dplyr::select(ALL,Bound,lengthDNA,modStart,modEnd), 
                         method =meth, trControl = train_control,metric = "ROC")
  # tuneGrid = expand.grid(C = seq(0.3, 3, length = 5)))#,  preProcess = c("center","scale"),)
  SVMout%<>%rbind(cbind(data.frame(model=meth,data="Mod"),
                        dplyr::select(as.data.frame(svmDNA$results),
                                      c(C,ROC,Sens,Spec))))
  # Using the original position
  svmDNA <- caret::train(Bound ~., 
                         data = dplyr::select(ALL,Bound,lengthDNA,StartP,EndP), 
                         method =meth, trControl = train_control,metric = "ROC")
  # tuneGrid = expand.grid(C = seq(0.3, 3, length = 5)))#,  preProcess = c("center","scale"),)
  SVMout%<>%rbind(cbind(data.frame(model=meth,data="Orig"),
                        dplyr::select(as.data.frame(svmDNA$results),
                                      c(C,ROC,Sens,Spec))))
}
saveRDS(SVMout,paste0(dir,"/Results/SVM.RData"))

################################################################################
######################## CONVOLUTIONAL NEURAL NET - CNN ########################
################################################################################
# Seeing as length of the DNA seems to be predictive of separating the bound and unbound data,
# let's use it later on!
# Model 1: one-hot 1D CNN model using a (4,x1) CNN layer filter
#          & Try out some different CNN structures
# Model 2: one-hot 1D CNN model using (2,x1) CNN layer filter then (1,x2)
# Model 3: one-hot*lengthDNA 1D CNN model
# Model 4: {one-hot,lengthDNA} 2D CNN model
# Other models: compare different #filters and other hyperparameters

################# DATA WRANGLING #################
# Min-max normalise the data
# Peaks:
process<-peaks%>%preProcess(method=c("range"))
peaks<-predict(process, peaks)
# Shuffle:
process<-shuffle%>%preProcess(method=c("range"))
shuffle<-predict(process, shuffle)
# Find the minimum width of arrays for CNN
maxL<-max(max(str_length(peaks$dna)),max(str_length(shuffle$dna))); print(maxL)
# One-hot encoding
OneH<-list(peaks=ConvOneHot(peaks$dna,maxL),
           shuffle=ConvOneHot(shuffle$dna,maxL))
# Calculate all the possible permutations of the one-hot array:
# permy<-gtools::permutations(n = 4, r = 4, v = 1:4)
# CNN hyperparameters 
Hyperparams<-list(cnnfilters=30,
                  kerneldim=c(4,30),
                  poolsize=2,
                  denselayers=30,
                  epocher=40,
                  SCV=5)

RunCNN<-function(OneH, Hyperparams){
  # How many bound & unbound sequences do we have to play with?
  lenP<-dim(OneH$peaks)[1]
  lenS<-dim(OneH$shuffle)[1]
  maxL<-dim(OneH$peaks)[3]
  # Setup the Stratified Cross-Validation method
  TrainTestValid<-floor(c(0.8,0.2)*lenP)
  # Output file
  performance<-data.frame()
  cnamers<-c("Loss","AUC","Recall","Precision","CVfold","j")
  for(j in 1:10){
    #@@@@@@@@@@@@@@@@@@@@@ STRATIFIED CROSS-VALIDATION @@@@@@@@@@@@@@@@@@@@@#
    for(cv in 1:Hyperparams$SCV){
      # Get all peaks and some shuffle data together and split into test, training and validation
      # Get a nSCV-column indices array that was generated from sampling from 1:lenP
      indies<-list(P=createFolds(1:lenP, k = Hyperparams$SCV, list = T, returnTrain = FALSE),
                   S=createFolds(1:lenS, k = Hyperparams$SCV, list = T, returnTrain = FALSE))
      # First work out the indices for both the peak and shuffle training datasets:
      Orig<-Data<-SortSCV(indies,OneH,cv,maxL,SCV=Hyperparams$SCV)
      # Run the model for each one-hot permutation:
      # for(ppp in 1:nrow(permy)){
      #   permSCV<-paste0("A=1,C=2,G=3,T=4 ordered by ",paste0(permy[ppp,],collapse = ","))
      #   Data<-Orig%<>%PermutateData(permy[ppp,])
      ################# CNN SECTION #################
      cnn_model <- keras_model_sequential() %>%
        # kernel_size taken from DeFine article - c(4,24)
        layer_conv_2d(filters = Hyperparams$cnnfilters, kernel_size = Hyperparams$kerneldim,
                      activation = 'relu', input_shape = c(4,maxL,1)) %>%
        layer_max_pooling_2d(pool_size = c(1, Hyperparams$poolsize)) %>%
        # layer_conv_2d(filters = Hyperparams$cnnfilters, kernel_size = c(2,10),
        #               activation = 'relu') %>%
        layer_dropout(rate = 0.2) %>%
        layer_batch_normalization() %>%
        layer_flatten() %>%
        layer_dense(units = Hyperparams$denselayers, activation = 'relu') %>%
        layer_dropout(rate = 0.5) %>%
        # layer_dense(units = 10, activation = 'relu') %>%
        # softmax to create probabilities which can build the AUC
        layer_dense(units = 2, activation = 'softmax') 
      # summary(cnn_model)
      # Compile it
      cnn_model %>% compile(
        loss = loss_categorical_crossentropy, # for some reason this was preferred over binary
        optimizer = optimizer_adam(), # could use optimizer_adadelta() or optimizer_adam() or optimizer_sgd()
        metrics = c(tf$keras$metrics$AUC(),tf$keras$metrics$Recall(),tf$keras$metrics$Precision())
      )
      # Fit the model!
      cnn_history <- cnn_model %>% fit(
        Data$X_Train, Data$Y_Train,
        # Batch size taken from DeFine
        batch_size = 139,
        # epochs taken from DeFine
        epochs = Hyperparams$epocher,
        validation_split = 0.0,
        verbose=0,
        callbacks = list(callback_early_stopping(monitor = "loss", patience = 5, restore_best_weights = TRUE))
      )
      # Check the performance on the test data
      tmp<-cnn_model%>%evaluate(Data$X_Test,Data$Y_Test)
      # Bind to the output file
      performance%<>%rbind(cbind(as.data.frame(as.list(t(tmp)),col.names=c(cnamers[1:4])),data.frame(CVfold=cv,j=j)))
      # }
    }
  }
  colnames(performance)<-cnamers
  
  return(performance)
}

output<-RunCNN(OneH, Hyperparams)

output%>%summarise(Loss=mean(Loss),AUC=mean(AUC),Recall=mean(Recall),Precision=mean(Precision))%>%
  saveRDS("./ChIPseq/Results/TwoCNNLayers.RData")