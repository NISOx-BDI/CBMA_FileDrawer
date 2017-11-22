#setwd("~/Desktop/HBM Publication Bias/Functional")

# Experiment info
Info = read.delim('func_exp_all_400.txt',na.strings='')
#names(Info) = c("BMapID","ExpID","","Journal","Year","","","","","","","","","","","","","Subjects")

# Read the lists with foci
Foci_tal = read.delim('func_tal_all_400.txt',head=F)
Foci_mni = read.delim('func_mni_all_400.txt',head=F)
names(Foci_tal) = c("BMapID","ExpID","x","y","z")
names(Foci_mni) = c("BMapID","ExpID","x","y","z")


# Remove resting state studies
keep = (Info$Activation)=='Y'
Info = Info[keep,]

# Find how many contrasts per publication
tmp = split(Info,Info$BMapID)
n_tmp = rep(0,length(tmp))
for (i in 1:length(tmp)) {
  n_tmp[i] = dim(tmp[[i]])[1]
}
summary(n_tmp)
rm(tmp,n_tmp)

# Find how many foci per publication
publications = as.numeric( names( table( Info$BMapID ) ) )
n_publications = length(publications)
n_tmp = rep(NA , length(levels(Info$BMapID)))
tmpID = NULL
for (i in 1:n_publications) {
  tmpID = ( Foci_tal[,1] == publications[i])
  tmpID = sum(tmpID)
  n_tmp[i] = sum(tmpID)
}
summary(n_tmp)

# Find how many foci per experiment
Counts = rep(0,dim(Info)[1])
tmpBMID = NULL
tmpEXPID = NULL
tmpID = NULL
for (i in 1:dim(Info)[1]) {
  tmpBMID = Info[i,1]
  tmpEXPID = Info[i,2]
  tmpID = (Foci_tal[,1]==tmpBMID)&(Foci_tal[,2]==tmpEXPID)
  Counts[i] = sum(tmpID)
}
Info$Counts = Counts



# Remove the Exp.Name, First_Author and PubMed variables
Info = Info[,-c(3,6,7)]


# Clean the categorical covariates
# Diagnosis
Info$Diagnosis = gsub(',.*','',Info$Diagnosis)
Info$Diagnosis = factor(Info$Diagnosis)
# Stumulus modality
Info$Stimulus_Mod = gsub(',.*','',Info$Stimulus_Mod)
Info$Stimulus_Mod = factor(Info$Stimulus_Mod)
# Stimulus Type
Info$Stimulus_Type = gsub(',.*','',Info$Stimulus_Type)
Info$Stimulus_Type = factor(Info$Stimulus_Type)
# Response Modality
Info$Response_Mod = gsub(',.*','',Info$Response_Mod)
Info$Response_Mod = factor(Info$Response_Mod)
# Response type
Info$Response_Type = gsub(',.*','',Info$Response_Type)
Info$Response_Type = factor(Info$Response_Type)
# Instructions
Info$Instructions = gsub(',.*','',Info$Instructions)
Info$Instructions = factor(Info$Instructions)
# Context
Info$Context = gsub(',.*','',Info$Context)
Info$Context = factor(Info$Context)
# Paradigm 
Info$Paradigm = gsub(',.*','',Info$Paradigm)
Info$Paradigm = factor(Info$Paradigm)





# Make 5 datasets that contain the covariate info and 5 datasets only with counts
##################################################################################################################################################################
By_Experiment = split(Info,Info$BMapID)
N = length(By_Experiment)
count_data = NULL
regression_data = list() 
for (i in 1:5) {
  regression_data[[i]] = Info[1,]
}
count_tmp = rep(0,5)
n_tmp = NULL
n_rem = NULL
set.seed(50001)
for (i in 1:N) {
  n = dim(By_Experiment[[i]])[1]
  if (n==1) {
    n_tmp = rep(1,5)
    for (j in 1:5) {
      regression_data[[j]] = rbind(regression_data[[j]],By_Experiment[[i]][n_tmp[j],])
      count_tmp[j] = By_Experiment[[i]]$Counts[n_tmp[j]]
    } 
  } else if ((n>1)&(n<5)) {
    n_tmp = sample( 1:n , n , replace=F )
    n_rem = sample( 1:n , 5-n , replace=T)
    n_tmp = c(n_tmp,n_rem)
    for (j in 1:5) {
      regression_data[[j]] = rbind(regression_data[[j]],By_Experiment[[i]][n_tmp[j],])
      count_tmp[j] = By_Experiment[[i]]$Counts[n_tmp[j]]
    } 
  } else if (n==5) {
    n_tmp = sample(1:n,5,replace=F)
    for (j in 1:5) {
      regression_data[[j]] = rbind(regression_data[[j]],By_Experiment[[i]][n_tmp[j],])
      count_tmp[j] = By_Experiment[[i]]$Counts[n_tmp[j]]
    }
  } else {
    n_tmp = sample(1:n,5,replace=F)
    for (j in 1:5) {
      regression_data[[j]] = rbind(regression_data[[j]],By_Experiment[[i]][n_tmp[j],])
      count_tmp[j] = By_Experiment[[i]]$Counts[n_tmp[j]]
    }
  }
  count_data = rbind(count_data,count_tmp)
}
# Remove the first line that was added only to create the list 
for (i in 1:5) {
  regression_data[[i]] = regression_data[[i]][-1,]
}




# Keep the variables of interest, removing data points with missing entries
for (i in 1:5) {
  regression_data[[i]] = regression_data[[i]][,c(1,2,4,6,11,12,15,16)]
}
# Remove the missing values
for (i in 1:5) {
  regression_data[[i]] = na.omit(regression_data[[i]])
}

# See how many have survived per subsample
K = rep(0,5)
for (i in 1:5) {
  K[i] = dim(regression_data[[i]])[1]
}



# Merge all categories that appear less than THRES times into a common catagory
THRES = 20
#categoricals = names(Info)[5:13]
columns = 4:6
col = tmp_var = NULL
# what follows will be performed for all datasets
for (kk in 1:5) { 
  K_tot = K[kk]
  for (k in columns) {
    col = k
    tab = table(regression_data[[kk]][,col])
    others = names(tab)[tab<THRES]
    n_others = length(others)
    if (n_others>0){
      levels(regression_data[[kk]][,col]) = c(levels(regression_data[[kk]][,col]),'Other')
      # repeat for all small levels
      for (i in 1:n_others) {
        # replace for current dataset
        for (j in 1:K_tot){
          if (regression_data[[kk]][j,col]==others[i]) {
            regression_data[[kk]][j,col] = 'Other'
            print(c(kk,k,j))
          }
        }
      }
      regression_data[[kk]][,col] = factor(regression_data[[kk]][,col])
    }
  }
}




# write the files
write.table(count_data,'count_data.txt',row.names=F,col.names=F)
write.table(regression_data[[1]],'regA.txt',row.names=F,col.names=T)
write.table(regression_data[[2]],'regB.txt',row.names=F,col.names=T)
write.table(regression_data[[3]],'regC.txt',row.names=F,col.names=T)
write.table(regression_data[[4]],'regD.txt',row.names=F,col.names=T)
write.table(regression_data[[5]],'regE.txt',row.names=F,col.names=T)





