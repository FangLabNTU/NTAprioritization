packsneed <- c('rcdk','rcdklibs','randomForest','leaps','caret','corrplot','tidyverse',
               'mlr','dplyr','Metrics','ggpubr','ggplot2','miceadds','rio','openxlsx','ggrepel')
lapply(packsneed, require, character.only = TRUE)

#set working directory
# currentpath <- "D:/R_untargeted analysis_project/ntaworkflow/datacombine/datacombine_function"
# setwd(currentpath)

#function 1 combine and clean data from input, predit the retention time and return the final data
getRTMSMS <- function(x){
  #select the datafolder
  x <- choose.dir(getwd(),caption = "Select folder")
  dm <- data.frame(matrix(nrow = 0,ncol = 13))
  colnames(dm) = c("Compound","Compound.ID","Description","Adducts","Formula","Score","Fragmentation.Score",
                   "Mass.Error..ppm.","Isotope.Similarity","Neutral.mass..Da.",
                   "m.z","Retention.time..min.","SMILES")
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
  dat3 = RTMSMS(dat2)
  dat3 %>% export(paste0(x,'clean_datawithRTMSMS.csv')) # done
  # return(dat3)
}

#function application
# getcleandata() #done








