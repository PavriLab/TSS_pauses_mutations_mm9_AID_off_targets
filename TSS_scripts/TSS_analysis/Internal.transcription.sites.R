# Internal.transcription.sites.R is a function called by various other scripts in the workflows.
# The input:
#   n; an integer that has to be 1 for now. There was a suggestion to update the script to scan in a window, but we used base-pair resolution
#   p; cutoff value (percentage) that will filter out signals with RPM of less than or equal to p*max(RPM). We didn't use this option too.
#   m; if m is set to 1, the function will produce some plots in Rstudio (not exported) while the script it running. I used it to test whether certain steps of the script work. It doesn't have to be used
#   super_region and super_region_size; if super_region is set to 1, then signals within the range of super_region_size will be grouped (summed) together and a single signal in the middle will be produced. This option wasn't used.
#   chromosome; the chromosome of the genomic region of the replicate to be scanned for signals
#   coordinate_start; the start coordinate of the genomic region to be scanned for signals
#   coordinate_end; the end coordinate of the genomic region to be scanned for signals
# The script will work on the replicate in BEDGRAPH format found in the environment with a variable name "test". The "test" variable (containing the BEDGRAPH replicate)
# is an absolute must.

# The output is

Internal.trascription.sites <- function(n,p,m,super_regions,super_region_size,chromosome,coordinate_start,coordinate_end){
  n <- as.integer(n) # let's go with a value of 1 for now
  p<-as.numeric(p) # cutoff percentage
  if (super_regions==1){
    super_region_size <- as.integer(super_region_size)
  }
  m <- as.integer(m) # to plot or not to plot
  chromosome <- chromosome
  coordinate_start <- as.integer(coordinate_start)
  coordinate_end <- as.integer(coordinate_end)

  test <- test[which(test$V1==chromosome),] # test will be inherited from the environment, it should be the BEDGRAPH file of a replicate
  results_list<-c(0,0,0,0,0) # initial vector that will be filled with the results of the scanning
  while (nrow(test[which(test$V2==coordinate_start),])==0 && coordinate_start<coordinate_end){ # if there's no signal in the start coordinate provided, then move one position
  # and stop if you reach the coordinate without finding any signal
    coordinate_start<-coordinate_start + 1
  }
  if (coordinate_start==coordinate_end){ # if you didn't find any signal and not the start=end, then report that there was no signal
    print('No signal in this region')
    results_list<-c(0,chromosome,NA,NA,0)
    return(results_list)
  } else if (coordinate_start<coordinate_end){ # now we do the same from the opposite side; starting from the end coordinate going backwards one position at a time until we find a signal
    while (nrow(test[which(test$V3==coordinate_end),])==0 && coordinate_end>coordinate_start){
      coordinate_end<-coordinate_end -1
    }
  }

  test <- test[which(test$V2==coordinate_start):which(test$V3==coordinate_end),] # now extract the regions between the new values of start and end coordinate from the test variable (the replicate)

  if (super_regions==1){
    par(mfrow=c(1,2))
  }
  # optional filtering step using p
  storage1<-data.frame(matrix(NA, nrow=1, ncol=4))
  max_peak<-max(abs(test[,4]))
  for (i in 1:(length(test[,1]))){
    if (abs(test[i,4])>0.01*p*max_peak){
      storage1[i,]<-test[i,]
      #storage1[i,4]<-sum(test[(i-n):(i+n),4]) #it doesn't takes into account +/- n bases from the peak, but +/- n PEAKS from the central peak, ignore it

    }
  }

  storage<-na.omit(storage1)
  colnames(storage) <- c("chr","position_start","position_end", "signal")
  rownames(storage) <- NULL

  if (m==1){
    plot(storage[,4], ylab="Signal Intensity", type='h', cex=1, col='blue')
  }
  if (super_regions==1){ # we didn't use super_regions, but see the description of what it can be used for
    if (length(storage[,2])>2){
      for (j in 2:(length(storage[,2]))){
        if (abs(storage[j,2]-storage[j-1,3]) < super_region_size){
          if (storage[j,1]==storage[j-1,1]){ #to avoid connecting regions of different chromosomes
            storage[j,4]<-(storage[j-1,4]+storage[j,4])
            storage[j,2]<-storage[j-1,2]
            storage[j-1,]<-NA
          }
        }
      }
    }
      storage<-na.omit(storage)
      rownames(storage)<-NULL
      if (m==1){
        plot(storage[,4], ylab="Signal Intensity", type='h', cex=1, col='blue')
      }
  }
  storage$chr <- sub("^", "chr", storage$chr)
  total_TSS<-as.matrix(length(storage[,1]))
  colnames(total_TSS)<-c('TSS')
  total_TSS
  results_list<-c(total_TSS,storage)
  return(results_list) # return the results with the signals found between the coordinate start and end
}
