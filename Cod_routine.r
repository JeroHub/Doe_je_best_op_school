# Set working directory
setwd("C:/Users/Jeroen Hubert/Dropbox/RAT files 2017/")
#setwd("C:/Users/Emmy/Dropbox/RAT files 2017/")

require(toal)

# Read dataset, standard fs (sampling rate) of HTI system is: 12000/s (12 bit)
dataset_HTI <- read.HTI.RAT("AT17S013410100.RAT", fs = 12000)

# Only load hydrophone 1, and all detections within first minute
start <- min(dataset_HTI$seconds)
length <- 60
segment <- 1
yaps.output <- NULL

## START LOOP
#for (start in min(dataset_HTI$seconds):max(dataset_HTI$seconds)) {
for (segment in 1:5){
dataset <- subset(dataset_HTI, seconds <= start+length & seconds >= start)
dataset.H1 <- subset(dataset, Hydrophone == 1)

## Find unique pulse rate intervals
pulseFrequencies <- TagFreq.detect(dataset.H1$seconds,n.tags = 1,plot = T,
                                 frequencies = seq(0.98,1.02,0.0001), minGap = 0.0001)
pulseFrequencies

dataset <- dataset[dataset$Hydrophone > 0, ]
hydro <- unique(dataset$Hydrophone)

# Initialize new dataframe for results
dataset.tag1 <- data.frame()
for(i in hydro){
  dataset.sub <- subset(dataset, subset = Hydrophone == i)
  tag.filter.idx <- TagFreq.filter(detections = dataset.sub$seconds,
                                   frequencies = pulseFrequencies[1],
                                   sensitivity = 0.1, plot = T)
  dataset.tag1 <- rbind(dataset.tag1,
                     data.frame(Tag = pulseFrequencies[1],
                          Hydrophone = i,
                          Seconds = dataset.sub$seconds[tag.filter.idx[[1]]],
                          DateTime = dataset.sub$DateTime[tag.filter.idx[[1]]]))
}
head(dataset.tag1)

dataset.periods <- TagFreq.label(data = dataset.tag1, plot = T)

# Select first detections per period and put them in a dataframe
dataset.YAPS_input <- NULL
Period <- min(dataset.periods$period)

#loop
for (Period in min(dataset.periods$period):max(dataset.periods$period)) {
  DateTime <- mean(dataset.periods$DateTime[dataset.periods$period == Period])
  H1 <- min(dataset.periods$Seconds[dataset.periods$Hydrophone == 1 & dataset.periods$period == Period])
  H2 <- min(dataset.periods$Seconds[dataset.periods$Hydrophone == 2 & dataset.periods$period == Period])
  H3 <- min(dataset.periods$Seconds[dataset.periods$Hydrophone == 3 & dataset.periods$period == Period])
  H4 <- min(dataset.periods$Seconds[dataset.periods$Hydrophone == 4 & dataset.periods$period == Period])
  dataset.YAPS_input.temp <- data.frame(Period, DateTime, H1, H2, H3, H4)
  
  if(exists('dataset.YAPS_input') && is.data.frame(get('dataset.YAPS_input'))){
    dataset.YAPS_input <- rbind(dataset.YAPS_input, dataset.YAPS_input.temp)
  }else{
    dataset.YAPS_input <- dataset.YAPS_input.temp
  }
  Period <- Period + 1
}

require(plyr)
dataset.YAPS_input <- rename(dataset.YAPS_input, c("H1"="Hydrophone 1", 
                                                   "H2"="Hydrophone 2",
                                                   "H3"="Hydrophone 3",
                                                   "H4"="Hydrophone 4"))
head(dataset.YAPS_input)

### START WITH POSITIONING USING YAPS!

# Speed of sound (with simualted error)
c <- 1500
# Add some error to c
c.real <- rnorm(n = 1, mean = 1500,sd = 1)

# Hydrophone positions        H1    H2    H3    H4
hydro.pos <- data.frame(x = c(7.5,  14,   7.5,  0),
                        y = c(14,   7.5,  0,    7.5),
                        z = c(3.35, 4.85, 3.35, 4.85))

toa.real <- dataset.YAPS_input[3:6]
toa.real$`Hydrophone 1`[!is.finite(toa.real$`Hydrophone 1`)] <- NA
toa.real$`Hydrophone 2`[!is.finite(toa.real$`Hydrophone 2`)] <- NA
toa.real$`Hydrophone 3`[!is.finite(toa.real$`Hydrophone 3`)] <- NA
toa.real$`Hydrophone 4`[!is.finite(toa.real$`Hydrophone 4`)] <- NA

# Spherical interpolation
#estimates.si <- TOA.localization(toa = toa.real,
#                                 hydrohpone.positions = hydro.pos,
#                                 c = c)

## YAPS-Pelagic
params <- list(
  logSigma_bi = -12.3, #log transformed SD of pulse intervals
  # logSigma_dl = rep(-10,length = nrow(hydro.pos)),		# Sigma for latency error
  logSigma_dl = -7.59,		# Sigma for latency error
  logSigma_toa = -15, # Time of arrival SD
  logSigma_x = -1.54,
  logSigma_y = -2,
  logSigma_z = -3.19)

yaps.output.temp <- NULL

#install.packages("tictoc")
require(tictoc)
tic();estimates.yaps.p <- yaps.pelagic(toa = toa.real,
                                       hydrophone.pos = hydro.pos,
                                       max.iterations = 5000, params = params); toc();

yaps.output.temp <- data.frame(estimates.yaps.p$XYZ)
yaps.output.temp$DateTime <- dataset.YAPS_input$DateTime
yaps.output.temp$segment <- segment
plot(yaps.output.temp$X1, yaps.output.temp$X2, type="l", xlim=c(0,12.5), ylim=c(0,12.5))
plot(yaps.output.temp$X1, yaps.output.temp$X3, type="l", xlim=c(0,12.5), ylim=c(0,5))



  if(exists('yaps.output') && is.data.frame(get('yaps.output'))){
    yaps.output <- rbind(yaps.output, yaps.output.temp)
  }else{
    yaps.output <- yaps.output.temp
  }

start <- start + length
segment <- segment + 1

# END LOOP
}

plot(yaps.output$X1, yaps.output$X2, type="l")
plot(yaps.output$X1, yaps.output$X3, type="l")

plot(yaps.output$X1, yaps.output$X2, type="l", xlim=c(0,12.5), ylim=c(0,12.5))
plot(yaps.output$X1, yaps.output$X3, type="l", xlim=c(0,12.5), ylim=c(0,5))

#plot(yaps.output$X1[yaps.output$segment == 3], yaps.output$X2[yaps.output$segment == 3], type="l", xlim=c(0,12.5), ylim=c(0,12.5))

