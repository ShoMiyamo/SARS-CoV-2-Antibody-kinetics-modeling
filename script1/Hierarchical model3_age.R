library(readr)
library(dplyr)
library(rstan)
library(ggplot2)
library(cowplot)

supdata<-read_csv('/Users/sho/Documents/R/R_output/20240814 numakura/20260202supplemental_data.csv')
Hmodel3 <- stan_model('/Users/sho/Documents/R/R_output/20240814 numakura/stan/antiS_exp.age_vax3_hierarchical M3.stan')

d3_stan<- supdata

#data curation

d3_stan$AGE0=as.numeric(case_when(d3_stan$Age<65~1,d3_stan$Age>=65~2))

d3_stan<-d3_stan %>% mutate(.,vax_age=as.numeric(case_when(
  vax==2&AGE0==1~1,vax==3&AGE0==1~2,vax==4&AGE0==1~3,vax==5&AGE0==1~3,
  vax==2&AGE0==2~4,vax==3&AGE0==2~5,vax==4&AGE0==2~6,vax==5&AGE0==2~6
)))

d3_stan<-d3_stan %>% filter(is.na(vax_age)==FALSE)

AGE<-unique(d3_stan[ , c('vax_age','AGE0')])$AGE0

data_age <-list(
  N=nrow(d3_stan),Y=d3_stan$anti_s,X=(d3_stan$intervals),
  K=6,S=2,VAX=as.numeric(d3_stan$vax_age),AGE=AGE,AGE0=d3_stan$AGE0)


#fitting
fit_age <- sampling(Hmodel3,
                    data=data_age,
                    iter = 9000,
                    warmup = 1000,
                    thin = 5,
                    chains = 4, core=4,
                    seed=1235,
                    control = list(adapt_delta = 0.99, max_treedepth=15))

#visualise
#same as "Hierarchical model3_infectionhistory.R"
