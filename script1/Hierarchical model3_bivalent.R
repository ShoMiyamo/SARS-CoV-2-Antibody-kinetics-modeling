library(readr)
library(dplyr)
library(rstan)
library(ggplot2)
library(cowplot)

supdata<-read_csv('supplemental_data.csv')
Hmodel3 <- stan_model('stan/hierarchical M3.stan')

d3_stan<- supdata

#data curation

d3_biv<-supdata2 %>% filter(as.numeric(vax)>=3  #exclude <3 vax n=13954
  )

d3_biv$bivalent<- case_when(d3_biv$bivalent=="bivalent vacccined"~2,TRUE~1)
d3_biv$infection_status0=case_when(d3_biv$infectionstatus3=="uninfected"~1,d3_biv$infectionstatus3=="Omicron"~2) #exclude pre-omicron n=89

d3_biv<-d3_biv %>% mutate(.,biv_infection=as.numeric(case_when(
  infection_status0==1&bivalent==1~1,infection_status0==1&bivalent==2~2,
  infection_status0==2&bivalent==1~3,infection_status0==2&bivalent==2~4
)))

d3_biv<- d3_biv %>% filter(is.na(biv_infection)==FALSE)
infection_status<-unique(d3_biv[ , c('biv_infection','infection_status0')])$infection_status0

data_bivalent <-list(N=nrow(d3_biv),Y=d3_biv$anti_s,X=as.numeric(d3_biv$general_intervals),K=4,S=2,VAX=as.numeric(d3_biv$biv_infection),AGE=infection_status,AGE0=d3_biv$infection_status0)


#fitting
fit_bivalent_exp2 <- sampling(antiS_exp.age_vax3,
                              data=data_bivalent,
                              iter = 9000,
                              warmup = 1000,
                              thin = 5,
                              chains = 4, core=4,
                              seed=1235,
                              control = list(adapt_delta = 0.99, max_treedepth=15)
                              )
#visualise
#same as "Hierarchical model3_infectionhistory.R"