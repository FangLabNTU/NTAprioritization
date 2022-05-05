packsneed <- c('rcdk','rcdklibs','randomForest','leaps','caret','mlr')
lapply(packsneed, require, character.only = TRUE)
#####RT prediction model functions###########
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

#define functions for molecular descriptors without preprocess
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

#get prediction from new data, #input as smiles lists
predictionRT <- function(input){
  descs <- getdesc_nopp(input) #get mds without any data imputation
  #find a way to load colnames data, settings,and model
  load.Rdata('ntaRTmodel.RData','predmodel')
  load.Rdata('ntaRTdescs.RData','descname')
  load.Rdata('ntaRTnormsettings.RData','settings')
  descsdata <- descs[,descname] #select mds to prepare same dimensions for normalzation settings
  descsdata <- predict(settings,descsdata) #get normalized mds for modeling
  predictionrt <- predict(predmodel, newdata = descsdata) #model select the mds as it trained
  output <- cbind('smiles'=input,'predictionRT'= predictionrt)
  return(output)
}
