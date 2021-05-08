##This is where we set up and train the neural network based on the dataset we generated
rm(list=ls())
library(spatstat)

library(hypergeo)
library(keras)

numBSs = 10
numBSi = 20
oSINR = 1:numBSs

#idices for the values in the dataset
mP = (numBSs+1):(2*numBSs)
Load = (2*numBSs+1):(3*numBSs)
Int =  (3*numBSs+1):(23*numBSs)
Network = (23*numBSs+1):(24*numBSs)
Distances= (24*numBSs+1):(25*numBSs)

##Load the dataset
load(file="trainingDataset.RData")
    
    
#model generation and training
    model <- keras_model_sequential()
    model %>%
      layer_dense(units=221,input_shape = c(221)) %>%
     layer_dense(units = 500, activation = 'relu') %>%
      layer_dense(units = 500, activation = 'relu') %>%
      layer_dense(units = 10, activation = 'softmax')
    
    
    
    model %>% compile(
      optimizer = 'adamax', 
      loss = 'categorical_crossentropy',
      metrics = c('accuracy')
    )
    
    #Take the training entries from the dataset
    train_images = results[1:480000,c(1:(numBSs*25+7))]
    
    
    ##Before beginning the training process, I (try to) normalise the data as much as possible, so it's as close to ranging between 0 and 1 as possible
    ##Maybe this isn't necessary, I'm not sure
    maxSINR = max(train_images[,(numBSs*25+7)])
    maxP = max(train_images[,c(2,mP)])
    minP = min(train_images[,c(2,mP)])
    maxInt = max(train_images[,c(4,Int)])
    minInt = min(train_images[,c(4,Int)])
    maxDist = max(train_images[,c(6,Distances)])
    minDist = min(train_images[,c(6,Distances)])
    for(i in 1:480000){
      
      ##normalise the observed SINR values
      mnfoo=min(train_images[i,oSINR])
      mxfoo=max(train_images[i,oSINR])
      if(mxfoo>0){
        train_images[i,oSINR]=(train_images[i,oSINR]-mnfoo)/(mxfoo-mnfoo)
      }  
      
    ##normalise the measured power  
    foo=min(train_images[i,mP])
    train_images[i,mP]= log10(train_images[i,mP]/foo)  
      
    mnfoo=min(train_images[i,mP])
    mxfoo=max(train_images[i,mP])
    if(mxfoo>0){
    train_images[i,mP]=(train_images[i,mP]-mnfoo)/(mxfoo-mnfoo)
    }
    
    ##normalise the interference score values
    mnfoo=min(train_images[i,Int])
    mxfoo=max(train_images[i,Int])
    if(mxfoo>0){
   # train_images[i,4]=(train_images[i,4]-mnfoo)/(mxfoo-mnfoo)
    train_images[i,Int]=(train_images[i,Int]-mnfoo)/(mxfoo-mnfoo)
    }
    
    #normalise the BS distances
    mnfoo=min(train_images[i,Distances])
    mxfoo=max(train_images[i,Distances])
    if(mxfoo>0){
  #    train_images[i,6]=(train_images[i,6]-mnfoo)/(mxfoo-mnfoo)
      train_images[i,Distances]=(train_images[i,Distances]-mnfoo)/(mxfoo-mnfoo)
    }
    
    }
    
    train_images[,Load]=train_images[,Load]/100

    ##normalise UAV height
    train_images[,(numBSs*25+3)]=train_images[,(numBSs*25+3)]/300
    
    #normalise BS height and building scale parameter
    train_images[,(numBSs*25+4)]=train_images[,(numBSs*25+4)]/50
    train_images[,(numBSs*25+5)]=train_images[,(numBSs*25+5)]/50
    
    #normalise beamwitdth
    train_images[,(numBSs*25+6)]=train_images[,(numBSs*25+6)]/pi
    
    #normalise current directional SINR
    train_images[,(numBSs*25+7)]=max(0,log10(train_images[,(numBSs*25+7)])/log10(maxSINR))
   
    
    #since we're not using the RB load parameters, get rid of them for the training
    train_images=train_images[,c(mP,Int,Distances,(numBSs*25+3))]
   
   #convert the labels to categorical format (ie ) 
    train_labels = to_categorical(results[1:480000,(numBSs*25+8)]-1,10)
    
    
    #now do all that with the test data
    test_images = results[480001:500000,c(1:(numBSs*25+7))]
    
    for(i in 1:20000){
      
      ##normalise the observed SINR values
      mnfoo=min(test_images[i,oSINR])
      mxfoo=max(test_images[i,oSINR])
      if(mxfoo>0){
        test_images[i,oSINR]=(test_images[i,oSINR]-mnfoo)/(mxfoo-mnfoo)
      }  
      
      ##normalise the measured power  
      foo=min(test_images[i,mP])
      test_images[i,mP]= log10(test_images[i,mP]/foo)  
      
      mnfoo=min(test_images[i,mP])
      mxfoo=max(test_images[i,mP])
      if(mxfoo>0){
        test_images[i,mP]=(test_images[i,mP]-mnfoo)/(mxfoo-mnfoo)
      }
      
      ##normalise the interference score values
      mnfoo=min(test_images[i,Int])
      mxfoo=max(test_images[i,Int])
      if(mxfoo>0){
        test_images[i,Int]=(test_images[i,Int]-mnfoo)/(mxfoo-mnfoo)
      }
      
      #normalise the BS distances
      mnfoo=min(test_images[i,Distances])
      mxfoo=max(test_images[i,Distances])
      if(mxfoo>0){
        test_images[i,Distances]=(test_images[i,Distances]-mnfoo)/(mxfoo-mnfoo)
      }
      
    }
    
    #normalise the load
    test_images[,Load]=test_images[,Load]/100
    
    ##normalise UAV height
    test_images[,(numBSs*25+3)]=test_images[,(numBSs*25+3)]/300
    
    #normalise BS height and building scale parameter
    test_images[,(numBSs*25+4)]=test_images[,(numBSs*25+4)]/50
    test_images[,(numBSs*25+5)]=test_images[,(numBSs*25+5)]/50
    
    #normalise beamwitdth
    test_images[,(numBSs*25+6)]=test_images[,(numBSs*25+6)]/pi
    
    #normalise current directional SINR
    test_images[,(numBSs*25+7)]=max(0,log10(test_images[,(numBSs*25+7)])/log10(maxSINR))
    
    
    #since we're not using the RB load parameters, get rid of them for the testing
    test_images=test_images[,c(mP,Int,Distances,(numBSs*25+3))]
    
    #convert the labels to categorical format (ie ) 
    test_labels = to_categorical(results[480001:500000,(numBSs*25+8)]-1,10)
    
    
    
    model %>% fit(train_images, train_labels, epochs = 50, validation_data=list(test_images, test_labels),batch_size=200)
    
    score <- model %>% evaluate(test_images, test_labels)
    
    cat('Test loss:', score$loss, "\n")
    cat('Test accuracy:', score$acc, "\n")
    
    
   model %>% save_model_weights_hdf5("NNweights.h5")
    