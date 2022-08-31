################################################################################
######################### SUPPORT VECTOR MACHINE - SVM #########################
################################################################################


################################################################################
######################## CONVOLUTIONAL NEURAL NET - CNN ########################
################################################################################
# Seeing as length of the DNA seems to be predictive of separating the bound and unbound data,
# let's use it later on!
# Model 1: one-hot 1D CNN model using a (4,x1) CNN layer filter
#          & Try out some different CNN structures
# Model 2: one-hot 1D CNN model using (2,x1) CNN layer filter then (1,x2)
# Model 3: augmented (flipped) one-hot 1D CNN model
# Model 4: one-hot*lengthDNA 1D CNN model
# Model 5: {one-hot,lengthDNA} 2D CNN model
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
# CNN hyperparameters 
epocher<-50
SCV<-5
# Setup the Stratified Cross-Validation method
TrainTestValid<-floor(c(0.8,0.2)*lenP)
# Calculate all the possible permutations of the one-hot array:
permy<-gtools::permutations(n = 4, r = 4, v = 1:4)
# Run the model for each one:
for(ppp in 1:nrow(permy)){
  permSCV<-paste0("A=1,C=2,G=3,T=4 ordered by ",paste0(permy[ppp,],collapse = ","))
  # How many bound sequences do we have to play with?
  lenP<-dim(OneH$peaks)[1]
  lenS<-dim(OneH$shuffle)[1]
  # Get all peaks and some shuffle data together and split into test, training and validation
  # Get a nSCV-column indices array that was generated from sampling from 1:lenP
  indies<-list(P=createFolds(1:lenP, k = SCV, list = T, returnTrain = FALSE),
               S=createFolds(1:lenS, k = SCV, list = T, returnTrain = FALSE))
  # Output file
  performance<-data.frame()
  cnamers<-c("Loss","AUC","Recall","Precision","permSCV","CVfold")
  #@@@@@@@@@@@@@@@@@@@@@ STRATIFIED CROSS-VALIDATION @@@@@@@@@@@@@@@@@@@@@#
  for(cv in 1:SCV){
    # First work out the indices for both the peak and shuffle training datasets:
    Data<-SortSCV(indies,OneH,cv,maxL,SCV=5)
    
    ################# CNN SECTION #################
    cnn_model <- keras_model_sequential() %>%
      # kernel_size taken from DeFine article - c(4,24)
      layer_conv_2d(filters = 32, kernel_size = c(4,24), activation = 'relu', input_shape = c(4,maxL,1)) %>% 
      layer_max_pooling_2d(pool_size = c(1, 4)) %>%
      layer_conv_2d(filters = 32, kernel_size = c(1,5), activation = 'relu') %>%
      layer_max_pooling_2d(pool_size = c(1, 2)) %>%
      layer_dropout(rate = 0.2) %>% 
      layer_batch_normalization() %>%
      layer_flatten() %>% 
      layer_dense(units = 128, activation = 'relu') %>% 
      layer_dropout(rate = 0.5) %>% 
      layer_dense(units = 128, activation = 'relu') %>% 
      layer_dense(units = 2, activation = 'softmax')
    
    summary(cnn_model)
    
    cnn_model %>% compile(
      loss = loss_categorical_crossentropy, # for some reason this was preferred over binary
      optimizer = optimizer_adam(), # could use optimizer_adadelta()
      metrics = c(tf$keras$metrics$AUC(),tf$keras$metrics$Recall(),tf$keras$metrics$Precision())
    )
    
    cnn_history <- cnn_model %>% fit(
      Data$X_Train, Data$Y_Train,
      # Batch size taken from DeFine
      batch_size = 139,
      # epochs taken from DeFine
      epochs = epocher,
      validation_split = 0.0
    )
    
    tmp<-cnn_model%>%evaluate(Data$X_Test,Data$Y_Test)
    
    performance%<>%rbind(data.frame(tmp,permSCV=permSCV,CVfold=cv))
    
    
  }
  colnames(performance)<-cnamers
}









