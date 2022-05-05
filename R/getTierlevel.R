##############matching and retrieving tox data from local database#######
packsneed <- c('tidyverse','dplyr','Metrics','ggpubr','ggplot2','miceadds','rio','openxlsx','ggrepel')
lapply(packsneed, require, character.only = TRUE)
####################### toxicity level ranking#######################
getTierlevel <- function(input){
  # select toxicity
  tox = readline(prompt = "please select the following toxicity:\nORLD50\nDMLC50\nDT\nFMLC50\nAM\nTPIGC50\nToxPi")
  if (tox != "ORLD50" && tox != "DMLC50" && tox != "FMLC50" &&
      tox != "DT" && tox != "AM" && tox != "TPIGC50" && tox != "ToxPi") {
    print("input error")
    return(F)
  }
  # # input candidate data
  input <- choose.files()
  candidate <- import(input)
  mdat <- merge(candidate, database[,c('SMILES',tox)], by='SMILES', all =FALSE)
  #remove na and return a reduction info
  dim1 = nrow(mdat)
  print(paste(dim1, "compounds are found"))
  mdat$Toxicity <- as.numeric(mdat[,c(tox)])
  mdat = mdat[!(is.na(mdat$Toxicity)| mdat$Toxicity==" "),]
  dim2 = nrow(mdat)
  print(paste(tox, "available for", dim2,"compounds"))
  # #rank tox level
  if (tox %in% c("FMLC50",'DMLC50','TPIGC50')){
    mdat = mdat %>% mutate(
      Tox_Level = case_when(
        Toxicity > 0 & Toxicity <= 2 ~ 3,
        Toxicity > 2 & Toxicity <= 6 ~ 2,
        Toxicity > 6 ~ 1,
        TRUE ~ 4
      )
    )
  }else if(tox == "ORLD50"){
    mdat = mdat %>% mutate(
      Tox_Level = case_when(
        Toxicity > 0 & Toxicity <= 1.5 ~ 3,
        Toxicity > 1.5 & Toxicity <= 3 ~ 2,
        Toxicity > 3 ~ 1,
        TRUE ~ 4
      )
    )
  }else if(tox %in% c("DT",'AM')){
    mdat = mdat %>% mutate(
      Tox_Level = case_when(
        Toxicity < 0.5 ~ 2,
        Toxicity >= 0.5 ~ 1,
        TRUE ~ 3
      )
    )
  }else if(tox == "ToxPi"){
    mdat =  mdat %>% mutate(
      Tox_Level = case_when(
        Toxicity > 0 & Toxicity <= 0.3 ~ 3,
        Toxicity > 0.3 & Toxicity <= 0.6 ~ 2,
        Toxicity > 1 ~ 1,
        TRUE ~ 4
      )
    )
  } #done
  #add priority ranking and plot
  mdat = mdat %>% mutate (
    Priority = case_when(
      RTMSMS == 1 & Tox_Level == 1 ~ 1,
      RTMSMS == 2 & Tox_Level == 1~ 2,
      RTMSMS == 3 & Tox_Level == 1~ 2,
      RTMSMS == 1 & Tox_Level == 2 ~ 3,
      RTMSMS == 1 & Tox_Level == 3 ~ 3,
      RTMSMS == 2 & Tox_Level == 1 ~ 3,
      RTMSMS == 3 & Tox_Level == 1 ~ 3,
      RTMSMS == 4 & Tox_Level == 1 ~ 3,
      (RTMSMS == 2 | RTMSMS == 3) & (Tox_Level == 2 | Tox_Level == 3) ~ 4,
      RTMSMS == 2 & Tox_Level == 2 ~ 4,
      RTMSMS == 2 & Tox_Level == 3 ~ 4,
      RTMSMS == 3 & Tox_Level == 2 ~ 4,
      RTMSMS == 3 & Tox_Level == 3 ~ 4,
      TRUE~5
    )
  )
  #implant plot function and pop the plot data
  mdat$rank_tox = NA
  mdat$rank_tox[order(-mdat$Toxicity)] <- 1:nrow(mdat)
  plotdat <- mdat %>%
    select(Toxicity, Priority,rank_tox, Compound.ID)
  toxplot <- ggplot(plotdat,aes(x=rank_tox,y=Toxicity))+
    geom_point(aes(color= Priority),size= 4) +
    scale_color_gradientn(colors = rainbow(5))+
    theme(legend.position = "right") +
    theme(text = element_text(size=20, face='bold'),
          axis.text.x = element_text(size=20,face='bold'),
          axis.text.y = element_text(size=20,face='bold'))+
    geom_text_repel(aes(x=rank_tox,label=ifelse(Priority==1,Compound.ID,"")),size=3,hjust=4,vjust=1,family = 'TT Arial',
                    colour='red',fontface = "bold",angle=0)
  ggsave(file=paste(input,"Toxicityplot_withpriority.tiff"), units="in",width=8, height=6,dpi=400, plot = toxplot)
  #tryout DONE
  mdat %>% export(paste0(input,'_tox_priority.csv'))
}
