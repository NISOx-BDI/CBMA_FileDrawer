# Load libraries 
library(gamlss)
library(gamlss.tr)

# Define the truncated families for GAMLSS
ZTDEL = trun(par=0,family=DEL, local=FALSE)
ZTNB = trun(par=0,family=NBI, local=FALSE)

# Remove trace 
control = gamlss.control(c.crit=0.001, n.cyc=100,trace = F) 

# Data 
n = read.table('counts.txt')[,1]

# Fit the two models
fit_nb  =  try ( gamlss(n ~ 1, family=ZTNB,control=control) )
fit_del =  try ( gamlss(n ~ 1, family=ZTDEL,control=control) )
