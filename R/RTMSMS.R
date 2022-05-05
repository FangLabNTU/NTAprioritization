##################RT-MSMSLEVEL grouping##########################
RTMSMS <- function(dat){
  # a <- choose.files()
  # dat <- import(a)
  dat = dat %>% mutate(
    RTMSMS = case_when(
      deltaRT <= 2 & Score >= 40 ~ 1,
      deltaRT <= 2 & Score >= 35 ~ 2,
      deltaRT <= 2 & Score >= 30 ~ 3,
      deltaRT >2 & deltaRT <= 4 & Score >= 30 ~3,
      deltaRT > 4~ 4,
      Score < 30 ~ 4,
      TRUE~5
    )
  )
  # dat %>% export(paste0(a,'RT-MSMSLEVEL.csv'))
  return(dat)
}
