packsneed <- c('rcdk','rcdklibs','randomForest','leaps','caret','corrplot','tidyverse',
               'mlr','dplyr','Metrics','ggpubr','ggplot2','miceadds','rio','openxlsx','ggrepel')
lapply(packsneed, require, character.only = TRUE)

########develop a function for model retrain with retenton time data from local machine
#function to correct the smiles to canonical smiles
canosmiles <- function(input){
  mols <- parse.smiles(input)
  targets <- unlist(lapply(mols,get.smiles,smiles.flavors(c('Canonical')))) #return smiles in list element
  res_list = list()
  for (i in targets) {
    res_list = c(res_list, i) #this will add i as list into the new list
  }
  return(unlist(res_list)) #unlist the result to get smiles as the elements in list
}

#define functions for molecular descriptors without preprocess #SMILES as the input
getdesc_nopp <- function(input){
  mols <- parse.smiles(input)
  descNames <- unique(unlist(sapply(get.desc.categories(),get.desc.names)))
  descs_tot <- data.frame()
  for (mol in mols) {
    descs_temp <- eval.desc(mol,descNames,verbose = FALSE)
    if (dim(descs_tot)[1] == 0) {
      descs_tot <- descs_temp
    }
    descs_tot <- rbind(descs_tot,descs_temp)
  }
  descs_start <- descs_tot[-c(1),]
  return(descs_start)
}

#define functions for molecular descriptors with preprocess
getdesc <- function(smilist) {
  mols <- parse.smiles(smilist)
  descNames <- unique(unlist(sapply(get.desc.categories(), get.desc.names))) #get 51 types of mds
  descs_tot <- data.frame()
  for (mol in mols) {
    descs_temp <- eval.desc(mol,descNames,verbose = FALSE)
    if (dim(descs_tot)[1] == 0) {
      descs_tot <- descs_temp
    }
    descs_tot <- rbind(descs_tot,descs_temp)
  }
  descs_start <- descs_tot[-c(1),]
  descs_start2<- descs_start[,!apply(descs_start, 2,function(x) any(is.na(x)))]  #220 left
  descs_clean2 <- descs_start2[,!apply(descs_start2,2,function(x) length(unique(x))==1)] #166 left
  r2 <- which(cor(descs_clean2)^2 >.9, arr.ind = TRUE)
  r2 <- r2[r2[,1]>r2[,2],]
  descs_output <- descs_clean2[,-unique(r2[,2])]
  return(descs_output)
}

#normalization by scaling with mean for each column
getnorm <- function(a,b){
  input_norm <- predict(preProcess(a), a)
  input_norm$RT <- b
  return(input_norm)
}

#get normalization settings
getnormsettings <- function(input){
  normsettings <- preProcess(input)
  return(normsettings)
}

#function for features selection after normalized data
getseldesc <- function(input){
  set.seed(123)
  smp_size <- round(nrow(input)*0.75)
  index <- sample(seq_len(nrow(input)), size = smp_size)
  train <- input[index, ]
  test <- input[-index, ]
  #rfe feature selection
  subsets <- c(5,10, 15, 20, 25,30,50)
  ctrl <- rfeControl(functions = rfFuncs,
                     method = 'repeatedcv',
                     repeats = 10,
                     verbose = FALSE)
  y<-train$RT
  x<-subset(train,select = -c(RT))
  lmprofile <- rfe(x,y,
                   sizes = subsets,
                   rfeControl = ctrl)
  seldesc <- predictors(lmprofile)
  return(seldesc)
}

#function for getting model
getmodel <- function(input,picname){
  #data splitting
  set.seed(123)
  smp_size <- round(nrow(input)*0.75)
  index <- sample(seq_len(nrow(input)), size = smp_size)
  train <- input[index, ]
  test <- input[-index, ]
  #rfe feature selection
  subsets <- c(10, 15, 20, 25,30,50)
  ctrl <- rfeControl(functions = rfFuncs,
                     method = 'repeatedcv',
                     repeats = 10,
                     verbose = FALSE)
  y<-train$RT
  x<-subset(train,select = -c(RT))
  lmprofile <- rfe(x,y,
                   sizes = subsets,
                   rfeControl = ctrl)
  #random forest modeling
  selmds <- predictors(lmprofile)
  x <- subset(x,select = c(selmds))
  rf.model = randomForest(y~., data=x,
                          mtry=9,
                          ntree =500,
                          importance = TRUE)
  #output model performance with prediction error
  test_y<-test$RT
  test_x<-subset(test,select = -c(RT))
  prd_rt = predict(rf.model, newdata = test_x)
  #data prep for plot
  prediction_table <- as.data.frame(cbind(prd_rt, test_y),row.names = NULL)
  #plot
  plot(lmprofile, type=c('g','o'), main='RMSE of 5cv', lwd=1.0, cex = 1.0)
  plot(rf.model)
  p <- ggscatter(prediction_table, x = "prd_rt", y = "test_y",
                 add = "reg.line", conf.int = TRUE,
                 cor.coef = TRUE, cor.method = "pearson",
                 xlab = "RT prediction (mins)", ylab = "RT experiment (mins)")
  p+ geom_text(x=8.5,y=14, label=paste('mae =',round(mae(prd_rt,test_y),2)),size = 5) #need to change x,y positions
  ggsave(paste (picname,"prediction vs experiment.png"))
  return(rf.model)
}

##model retrain function
RTmodelretrain = function(input){
  ## input candidate data
  input <- choose.files()
  candidate <- import(input)
  smilst = canosmiles(candidate$SMILES)
  RT = candidate$RT
  datadescs = getdesc(smilst)
  norm_data = getnorm(datadescs,RT)
  RTmodel = getmodel(norm_data, picname='RT_model')
  #save model settings
  RTmodel_descs <- colnames(datadescs) #descripotrs before norm
  RTmodel_descnorsettings <- getnormsettings(datadescs)
  save(RTmodel,file = file.path(getwd(),'ntaRTmodel2.RData'))
  save(RTmodel_descs,file = file.path(getwd(),'ntaRTdescs2.RData'))
  save(RTmodel_descnorsettings,file = file.path(getwd(),'ntaRTnormsettings2.RData'))
}
