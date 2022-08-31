################################################################################
######################### SUPPORT VECTOR MACHINE - SVM #########################
################################################################################


################################################################################
######################## CONVOLUTIONAL NEURAL NET - CNN ########################
################################################################################
# Get the data in:
# Find the minimum width of arrays for CNN
maxL<-max(max(str_length(peaks$dna)),max(str_length(shuffle$dna))); print(maxL)
# One-hot encoding
OneHpeaks<-ConvOneHot(peaks$dna,maxL)
OneHshuffle<-ConvOneHot(shuffle$dna,maxL)

# STRATIFIED CROSS-VALIDATION

# Seeing as length of the DNA seems to be predictive of separating the bound and unbound data,
# let's use it later on!
# Model 1: one-hot 1D CNN model using a (4,x1) CNN layer filter
# Model 2: one-hot 1D CNN model using (2,x1) CNN layer filter then (1,x2)
# Model 3: augmented (flipped) one-hot 1D CNN model
# Model 4: one-hot*lengthDNA 1D CNN model
# Model 5: {one-hot,lengthDNA} 2D CNN model

# Other models: compare different #filters and other hyperparameters

# Hyperparameters 
epocher<-50
optimmy<-"optimizer_adam"
SCV<-4
# Normalise LengthDNA column


# Get rid of validation bit
# Use test data & predict function to build a TPR & FPR, including all the remaining shuffle data
# Make sure to use all of the shuffle data across the SCV splitting

# Try out some different CNN structures
# Then move on to the different models that we want to try


performance<-data.frame()
cnamers<-c("loss_Test","AUC_Test","loss_Train","AUC_Train","Algorithm","permSCV","repSCV")
# For SCV, we can simply permutate the arrays wrt the bases



# How many bound sequences do we have to play with?
lenP<-dim(OneHpeaks)[1]
lenS<-dim(OneHshuffle)[1]
# Get all peaks and some shuffle data together and split into test, training and validation
# Setup the Stratified validation data
TrainTestValid<-floor(c(0.7,0.15,0.15)*lenP)
# Templates
X_Validate<-array(dim=c(TrainTestValid[3]*2L,4,maxL))
Y_Validate<-rep(NA,TrainTestValid[3]*2L)
X_Train<-array(dim=c(TrainTestValid[1]*2L,4,maxL))
Y_Train<-rep(NA,TrainTestValid[1]*2L)
X_Test<-array(dim=c(TrainTestValid[2]*2L,4,maxL))
Y_Test<-rep(NA,TrainTestValid[2]*2L)
##################### VALIDATION #####################
# Randomly sample some values from this array
indies<-sample(1:lenP,TrainTestValid[3]); arrI<-rep(F,lenP); arrI[indies]<-T
# Save into validation data and reduce the peaks array
X_Validate[1:TrainTestValid[3],,]<-OneHpeaks[arrI,,]
Y_Validate[1:TrainTestValid[3]]<-1
OneHpeaks<-OneHpeaks[!arrI,,]
# Randomly sample some values from this array
indies<-sample(1:lenS,TrainTestValid[3]); arrI<-rep(F,lenS); arrI[indies]<-T
# Save into validation data and reduce the peaks array
X_Validate[(TrainTestValid[3]+1):(2L*TrainTestValid[3]),,]<-OneHshuffle[arrI,,]
Y_Validate[(TrainTestValid[3]+1):(2L*TrainTestValid[3])]<-0
OneHshuffle<-OneHshuffle[!arrI,,]
################### TEST & TRAINING ###################
lenP<-dim(OneHpeaks)[1]
lenS<-dim(OneHshuffle)[1]
# Randomly sample some values from this array - TRAINING
indies<-sample(1:lenP,TrainTestValid[1]); arrI<-rep(F,lenP); arrI[indies]<-T
# Save into data frames
X_Train[1:TrainTestValid[1],,]<-OneHpeaks[arrI,,]
if(sum(!arrI)>TrainTestValid[2]) arrI[sample((1:lenP)[!arrI],1)]<-T
X_Test[1:TrainTestValid[2],,]<-OneHpeaks[!arrI,,]
Y_Train[1:TrainTestValid[1]]<-1
Y_Test[1:TrainTestValid[2]]<-1
# Randomly sample some values from this array
indies<-sample(1:lenS,TrainTestValid[1]); arrI<-rep(F,lenS); arrI[indies]<-T
# Save into data frames
X_Train[(1+TrainTestValid[1]):(2L*TrainTestValid[1]),,]<-OneHshuffle[arrI,,]
if(sum(!arrI)>TrainTestValid[2]) arrI[sample((1:lenS)[!arrI],(sum(!arrI)-TrainTestValid[2]))]<-T
X_Test[(1+TrainTestValid[2]):(2L*TrainTestValid[2]),,]<-OneHshuffle[!arrI,,]
Y_Train[(1+TrainTestValid[1]):(2L*TrainTestValid[1])]<-0
Y_Test[(1+TrainTestValid[2]):(2L*TrainTestValid[2])]<-0
# Clean-up
rm(indies,arrI,OneHpeaks,OneHshuffle)
# Convert to categorical
Y_Train<-to_categorical(Y_Train, 2)
Y_Test<-to_categorical(Y_Test, 2)
Y_Validate<-to_categorical(Y_Validate, 2)
# Reshape for TensorFlow
X_Train <- array_reshape(X_Train, c(dim(X_Train)[1], 4, maxL, 1))
X_Test <- array_reshape(X_Test, c(dim(X_Test)[1], 4, maxL, 1))

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
  optimizer = match.fun(optimmy) , # could use optimizer_adadelta()
  metrics = tf$keras$metrics$AUC()
)

cnn_history <- cnn_model %>% fit(
  X_Train, Y_Train,
  # Batch size taken from DeFine
  batch_size = 139,
  # epochs taken from DeFine
  epochs = epocher,
  validation_split = 0.0
)

tmp<-cnn_model %>% evaluate(X_Test, Y_Test)
tmp<-c(unlist(unname(tmp)),cnn_history$metrics$loss[epocher],cnn_history$metrics$auc_1[epocher],optimmy)

performance%<>%rbind(tmp)



colnames(performance)<-cnamers










