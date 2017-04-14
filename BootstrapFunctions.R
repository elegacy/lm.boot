################################################################
####EfronBoot is a function which takes bootstrap samples   ####
#### from an input data vector and returns the bootstrap    ####
#### sample  means and bootstrap samples.                   ####
####                                                        ####
####  data:  vector of data (should contain no missing vals ####
####         and be quantitative
####  B:  number of bootstrap samples to take.  Default=100 ####
####                                                        ####
#### returns bootstrapMeans:  vector of bootstrap sample means
####         bootstrapSamples: matrix containing bootstrap  ####
####                           samples.  Each column is a   ####
####                           single sample.               ####
################################################################
EfronBoot <- function(data, B=100){
  #Function should check inputs are correct
  bootData <- matrix(NA, nrow=length(data), ncol=B) #matrix to fill with bootstrap samples
  for(b in 1:B){
    bootData[ , b] <- sample(data, size=length(data), replace=TRUE) #SRSWR from orig. data
  }
  meanBootData <- apply(bootData, 2, mean) #get mean from each bootstrap sample
  return(list(bootstrapMeans =meanBootData, bootstrapSamples=bootData))
}



################################################################
##### BootCI is a function which performs textbook Method  #####
#####   1 from an input data vector (of bootstrapped means)#####
#####                                                      #####
#####  bootstrapMeans:  vector of bootstrap sample means   #####
#####  alpha:  significance level for the interval         #####
#####  type:  "one.sided.low", "one.sided.up", "two.sided" #####
#####     specifies interval type                          #####
#####                                                      #####
#####  returns upper, lower:  bounds of interval           #####
#####    (also prints these to the console)                #####
#####                                                      #####
#####  this function is essentially the same as 'quantile' #####
################################################################
BootCI <- function(bootstrapMeans, alpha=0.05, type="two.sided"){
   #Function should check inputs are correct
   sortedMeans <- bootstrapMeans[order(bootstrapMeans)]
   if(type=="two.sided"){
     upper <- sortedMeans[ceiling(length(bootstrapMeans)*(1-alpha/2))]
     lower <- sortedMeans[floor(length(bootstrapMeans)*(alpha/2))]
   }else if(type=="one.sided.low"){
     upper <- Inf
     lower <- sortedMeans[floor(length(bootstrapMeans)*(alpha))]
   }else{
     upper <- sortedMeans[ceiling(length(bootstrapMeans)*(1-alpha))]
     lower <- -Inf
   }
   print(paste("(", round(lower,2), ", ", round(upper,2), ")", sep='')); flush.console();
   return(list(lower=lower, upper=upper))
}





################################################################
##### Samp2Boot is a function which performs bootstrap on  #####
#####  2 samples of data (paired or independent).  It      #####
#####  returns a vector of the bootstrap mean differences. #####
#####                                                      #####
##### data1:  vector of sample data in group 1             #####
##### data2:  vector of sample data in group 2             #####
##### dataType:  either "independent" or "paired"          #####
##### B:  number of bootstrap samples to take              #####
#####                                                      #####
##### returns bootstrapDiffMeans:  vector of differences in#####
#####   means from each bootstrap sample                   #####
################################################################
Samp2Boot <- function(data1, data2, dataType, B=100){
  #Function should check inputs are correct
  if(dataType=="paired"){
    dataDiff <- data1-data2
    bootstrapDiffMeans <- rep(NA, B)
    for(b in 1:B){
      bootSample <- sample(dataDiff, size=length(dataDiff), replace=TRUE) #SRSWR from orig. data differences
      bootstrapDiffMeans[b] <- mean(bootSample)                           #mean difference in bootstrap sample
    }
  }
  else if(dataType=="independent"){
    bootstrapDiffMeans <- rep(NA, B)
    for(b in 1:B){
      bootData1Samp <- sample(data1, size=length(data1), replace=TRUE) #SRSWR from orig. data1
      bootData2Samp <- sample(data2, size=length(data2), replace=TRUE) #SRSWR from orig. data2
      bootstrapDiffMeans[b] <- mean(bootData1Samp) - mean(bootData2Samp)
    }
  }
  return(list(bootstrapDiffMeans=bootstrapDiffMeans))
}

