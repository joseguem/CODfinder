# @author : Jose Guerrero (jose.guerrero@cabimer.es)
# @brief : CODfinder is a TOPDOM-based software for identifying Coexpression Domains (CODs) from a expression matrix. 

# @fn codfinder
# @param inputfile : string, inputFile Address,
# - (N+1) * M, where N is the number of RNA-seq samples and M is the number of genes.
# - Inclue a Header with genes names.
# @param window.size :integer, number of bins to extend.
# @param outfile : string, outfile address to write
# - 3 files are created: *.correl.csv, *.binsignal.csv, *.COD.csv


codfinder <- function(inputfile, window.size, outfile)
{

options(scipen=100, digits=10)

## Pearson Correlation Matrix
exp <- read.table(inputfile, sep="\t", header=TRUE)
correl <- cor(exp)
correl <- as.matrix(correl)
correl[is.na(correl)==TRUE]=0

##TOPDOM variables definition
n_bins = nrow(correl)
mean.cf <- rep(0, n_bins)
pvalue <- rep(1, n_bins)
matrix.data <- as.matrix( correl[, (ncol(correl) - nrow(correl)+1 ):ncol(correl)] )



### PART1: Binsignal file

##Diamond function from TOPDOM for creating binsignal file
for(i in 1:(n_bins-1))
{
  lowerbound = max(1, i-window.size+1)
  upperbound = min(i+window.size, n_bins)
  diamond = matrix.data[lowerbound:i, (i+1):upperbound]
  mean.cf[i] = mean(diamond)
}


##P-value
for(i in window.size:(n_bins-window.size+1))
{
  upstream = mean.cf[c((i-(window.size-1)):i)]
  downstream = mean.cf[c((i+1):(i+window.size))]
  ttest = t.test(upstream,downstream)
  pvalue[i] = ttest$p.value 
  if(!is.nan(pvalue[i])) pvalue[i]=pvalue[i]
  else pvalue[i]=1
}

for (i in 1:(n_bins-1))
{
  if(!is.nan(pvalue[i])) pvalue[i]=pvalue[i]
  else pvalue[i]=1
}

##Binsignal file
bin.signal <- matrix(data = c(mean.cf, pvalue),nrow=n_bins,ncol=2) 




##PART2: CODs Identification from binsignal file.
length <- nrow(bin.signal)
vector <- rep(0,length)

##DOMAIN vector
i=2
while (i < (length-window.size))
{
  if (bin.signal[i,1]>0.15 & sum(as.numeric(bin.signal[i:(i+3),1]>0.15))>2 & bin.signal[i-1,2]<0.05) 
  {
    vector[i-1]= "START"
    while((i < (length-window.size)) & (bin.signal[i+1,1]>0.15 | bin.signal[i+2,1]>0.15 | bin.signal[i+3,1]>0.15 | bin.signal[i+4,1]>0.15 ) & (bin.signal[i+1,1]>0.15| bin.signal[i+2,1]>0.15 | bin.signal[i,2]>0.05))
    {
      i=i+1
      vector[i-1]="DOMAIN"
      vector[i]="DOMAIN"
      if (bin.signal[i+1,1]<0.15 & bin.signal[i+1,2]<0.05) vector[i+1]="END"
      else if (bin.signal[i,1]>0.15 & bin.signal[i,2]<0.05) vector[i]="END"
      else vector[i]="END"
      
    }
    i=i+1
  }
  else i=i+1
}

##Number of CODs defined as START number
j=1
logic.vector <- rep(0,length)
while(j<length){
  logic.vector[j] = vector[j]=="START"
  j=j+1
}
nstart <- sum(as.logical(logic.vector))

##COD Matrix
cod.matrix <- matrix(NA,nrow=nstart,ncol=2)
k=1
l=1
while (k < length){
  if (vector[k]=="START")
  {
    cod.matrix[l,1]=k
    k=k+1
  }
  else if (vector[k]=="END")
  {
    cod.matrix[l,2]=k
    k=k+1
    l=l+1
  }
  else k=k+1
}



##PART3: Filtering COD Matrix


#DOMAIN vector filtering
i=1
cod.matrix.filtered <- matrix(NA,nrow=nrow(cod.matrix),ncol=2)
while(i<=nrow(cod.matrix))
{
  k=cod.matrix[i,2]
  if (bin.signal[k,2]<0.05 | bin.signal[k-1,2]<0.05 | bin.signal[k-2,2]<0.05 | bin.signal[k+1,2]<0.05 | bin.signal[k+2,2]<0.05 | bin.signal[k+3,2]<0.05)
  {
    cod.matrix.filtered[i,1]= cod.matrix[i,1]
    cod.matrix.filtered[i,2]=cod.matrix[i,2]
    i=i+1
  } else i=i+1
}

if(cod.matrix.filtered[nrow(cod.matrix.filtered),2]-cod.matrix.filtered[nrow(cod.matrix.filtered),1] >= 4)
{
  cod.matrix.filtered[nrow(cod.matrix.filtered),1]=cod.matrix.filtered[nrow(cod.matrix.filtered),1]
  cod.matrix.filtered[nrow(cod.matrix.filtered),2]=cod.matrix.filtered[nrow(cod.matrix.filtered),2]
} 

else
{
  cod.matrix.filtered[nrow(cod.matrix.filtered),1]= NA
  cod.matrix.filtered[nrow(cod.matrix.filtered),2] = NA
}

j=1
k=0
idx.na <- rep(NA, nrow(cod.matrix.filtered))
while(j<=nrow(cod.matrix.filtered))
{
  if (is.na(cod.matrix.filtered[j,1]))
  {
    k=k+1
    idx.na[k]=j
    j=j+1
  } else j=j+1
}
idx.na <- idx.na[!is.na(idx.na)]


if(is.numeric(idx.na)==FALSE) 
{
  cod.matrix.filtered2 <- matrix(NA,nrow=nrow(cod.matrix.filtered),ncol=2)
  cod.matrix.filtered2 <- cod.matrix.filtered
} else cod.matrix.filtered2 <- cod.matrix.filtered[-idx.na,]

##Write files
outCorrel <- paste(outfile, ".correl.csv",sep="")
write.csv(correl, outCorrel)

colnames(bin.signal) <- c("mean.cf", "pvalue") 
outBinsignal <- paste(outfile, ".binsignal.csv",sep="") 
write.csv(bin.signal,outBinsignal) 

colnames(cod.matrix.filtered2) <- c("START","END")
outCOD_filt <- paste(outfile, ".COD.csv",sep="")
write.csv(cod.matrix.filtered2, outCOD_filt)

}