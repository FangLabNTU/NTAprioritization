packsneed <- c('rcdk','rcdklibs','randomForest','leaps','caret','corrplot','tidyverse',
               'mlr','dplyr','Metrics','ggpubr','ggplot2','miceadds','rio','openxlsx','ggrepel')
lapply(packsneed, require, character.only = TRUE)

#set working directory
# currentpath <- "D:/R_untargeted analysis_project/ntaworkflow/datacombine/datacombine_function"
# setwd(currentpath)

#PRE FUNC for canonical SMILES
canosmiles <- function(input){
  mols <- parse.smiles(input)
  targets <- unlist(lapply(mols,get.smiles,smiles.flavors(c('Canonical')))) #return smiles in list element
  res_list = list()
  for (i in targets) {
    res_list = c(res_list, i) #this will add i as list into the new list
  }
  return(unlist(res_list)) #unlist the result to get smiles as the elements in list
}

##################RT-MSMSLEVEL grouping##########################
#setup a userdefine input for RT
RTMSMS2 <- function(dat, RT1, RT2, score1, score2, score3){
  # a <- choose.files()
  # dat <- import(a)
  dat = dat %>% mutate(
    RTMSMS = case_when(
      deltaRT <= RT1 & Score >= score1 ~ 1,
      deltaRT <= RT1 & Score >= score2 ~ 2,
      deltaRT <= RT1 & Score >= score3 ~ 3,
      deltaRT >RT1 & deltaRT <= RT2 & Score >= score3 ~3,
      deltaRT > RT2~ 4,
      Score < score3 ~ 4,
      TRUE~5
    )
  )
  # dat %>% export(paste0(a,'RT-MSMSLEVEL.csv'))
  return(dat)
}

#function 1 combine and clean data from input, predit the retention time and return the final data
getRTMSMS2 <- function(x){
  #select the datafolder
  x <- choose.dir(getwd(),caption = "Select folder")
  dm <- data.frame(matrix(nrow = 0,ncol = 13))
  colnames(dm) = c("Compound","Compound.ID","Description","Adducts","Formula","Score","Fragmentation.Score",
                   "Mass.Error..ppm.","Isotope.Similarity","Neutral.mass..Da.",
                   "m.z","Retention.time..min.","SMILES")
  #input defined threshold for RT and Score
  RT1 = readline(prompt = "please input a low threshold for RT\n")
  RT2 = readline(prompt = "please input a high threshold for RT\n")
  score1 = readline(prompt = "please input a highe threshold for Score\n")
  score2 = readline(prompt = "please input a medium threshold for Score\n")
  score3 = readline(prompt = "please input a low threshold for Score\n")
  RT1 = as.numeric(RT1)
  RT2 = as.numeric(RT2)
  score1 = as.numeric(score1)
  score2 = as.numeric(score2)
  score3 = as.numeric(score3)

  data_name <- list.files(x)
  for (i in 1:length(data_name)){
    path_file <- paste0(x,'\\',data_name[i])
    data = import(path_file)
    data1 = data %>%
      select("Compound","Compound.ID","Description","Adducts","Formula","Score","Fragmentation.Score",
             "Mass.Error..ppm.","Isotope.Similarity","Neutral.mass..Da.",
             "m.z","Retention.time..min.","SMILES")
    data2 = data1 %>% filter("Fragmentation.Score" > 0)
    ##data combination
    dm = rbind(dm,data2)
  }
  ##output the combined data
  #clean data, select the maximum scores
  #add canonical smiles
  # dm$SMILES = canosmiles(dm$SMILES)
  dat1 = dm %>% group_by(Compound,SMILES)
  dat2 = dat1 %>% filter(Score == max(Score))
  dat2$Score = ceiling(dat2$Score)
  dat2 <- dat2[!duplicated(dat2[,c('Compound','SMILES')]),]
  #retention time prediction
  smilst <- canosmiles(dat2$SMILES)
  newrt <- predictionRT(smilst)
  dat2$PredRT <- as.numeric(newrt[,2])
  #remove na values in PredRT
  dim1 <- nrow(dat2)
  dat2 %>% drop_na(PredRT)
  dim2 <- nrow(dat2)
  print(paste("drop", dim1-dim2, "compounds"))
  # calculation of the deltaRT
  dat2 = mutate(dat2, deltaRT = abs(Retention.time..min.- PredRT))
  #RTMSMSlevel
  dat3 = RTMSMS(dat2,RT1, RT2, score1, score2, score3)
  dat3 %>% export(paste0(x,'clean_datawithRTMSMS.csv')) # done
  # return(dat3)
}

#function application
# getcleandata() #done








