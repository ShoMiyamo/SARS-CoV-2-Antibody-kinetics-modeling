library(readr)
library(dplyr)
library(rstan)
library(ggplot2)
library(cowplot)

supdata<-read_csv('supplemental_data.csv')
Hmodel3 <- stan_model('hierarchical M3.stan')

d3_stan<- supdata

#data curation
d3_stan$sex0=d3_stan$male_sex0+1
d3_stan<-d3_stan %>% mutate(.,vax_sex=as.numeric(case_when(
  vax==2&sex0==1~1,vax==3&sex0==1~2,vax==4&sex0==1~3,vax==5&sex0==1~3,
  vax==2&sex0==2~4,vax==3&sex0==2~5,vax==4&sex0==2~6,vax==5&sex0==2~6,
)))
d3_stan<-d3_stan %>% filter(is.na(vax_sex)==FALSE)

Sex<-unique(d3_stan[ , c('vax_sex','sex0')])$sex0

data_sex <-list(N=nrow(d3_stan),Y=d3_stan$anti_s,X=(d3_stan$intervals),K=6,S=2,VAX=as.numeric(d3_stan$vax_sex),AGE=Sex,AGE0=d3_stan$sex0)

#fitting
fit_sex <- sampling(Hmodel3,
                    data=data_sex,
                    iter = 9000,
                    warmup = 1000,
                    thin = 5,
                    chains = 4, core=4,
                    seed=1235,
                    control = list(adapt_delta = 0.99, max_treedepth=15))

#visualise
#same as "Hierarchical model3_infectionhistory.R"
