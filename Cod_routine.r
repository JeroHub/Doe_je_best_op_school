# Set working directory
# You dont need to do this working directory stuff if you make a Rproj file (project)
#setwd("C:/Users/Jeroen Hubert/Dropbox/RAT files 2017/")
#setwd("C:/Users/Emmy/Dropbox/RAT files 2017/")

require(toal)

# Read dataset, standard fs (sampling rate) of HTI system is: 12000/s (12 bit)
dataset_HTI <- read.HTI.RAT("AT17S013410100.RAT", fs = 12000)
head(dataset_HTI)

## Set segment length
length <- 120 # In seconds
segments <- 5 # How many segments to analyse (set to 0 for all)
	
## Loop variables
start <- min(dataset_HTI$seconds)
end <- max(dataset_HTI$seconds)
max.segments <- floor((end - start)/60)
if (segments == 0){
	segments <- max.segments
}else{
	segments <- min(segments, max.segments)
}

## YAPS-Starting parameters
params <- list(
	logSigma_bi = -12.3, #log transformed SD of pulse intervals
	# logSigma_dl = rep(-10,length = nrow(hydro.pos)),		# Sigma for latency error
	logSigma_dl = -7.59,		# Sigma for latency error
	logSigma_toa = -15, # Time of arrival SD
	logSigma_x = -1.54,
	logSigma_y = -2,
	logSigma_z = -3.19)

yaps.output <- data.frame()
for(segment in 1:segments){
	
	## Subset data segment
	segStart <- (length*(segment-1)) + start
	segEnd <- segStart + length
	dataset <- subset(dataset_HTI,
										seconds >=  + segStart &
										seconds < segEnd)
	dataset.H1 <- subset(dataset, Hydrophone == 1)
	
	message('Segment: ', segment, ', Analysis Period: ',
					round(segStart, digits = 1), '-',round(segEnd,digits = 1))
	
	## Find unique pulse rate intervals
	pulseFrequencies <- TagFreq.detect(dataset.H1$seconds,n.tags = 1,plot = F,
																		 frequencies = seq(0.98,1.02,0.0001), 
																		 minGap = 0.0001)
	message('Pulse Frequency: ',pulseFrequencies[1])
	dataset <- dataset[dataset$Hydrophone > 0, ]
	hydro <- unique(dataset$Hydrophone)
	
	
	## Filter out detections (based on tag pulse period)
	# Initialize new dataframe for results
	dataset.tag1 <- data.frame()
	for(i in hydro){
		dataset.sub <- subset(dataset, subset = Hydrophone == i)
		tag.filter.idx <- TagFreq.filter(detections = dataset.sub$seconds,
																		 frequencies = pulseFrequencies[1],
																		 sensitivity = 0.1, plot = F)
		dataset.tag1 <- rbind(dataset.tag1,
													data.frame(Tag = pulseFrequencies[1],
																		 Hydrophone = i,
																		 Seconds = dataset.sub$seconds[tag.filter.idx[[1]]]))
	}
	print(head(dataset.tag1))
	
	
	## Assign periods to tag detections
	dataset.periods <- TagFreq.label(data = dataset.tag1, plot = T)
	dataset.YAPS_input <- NULL
	### Loop all periods, and save first detection per hydrophone
	for (Period in min(dataset.periods$period):max(dataset.periods$period)) {
		dataset.YAPS_input.temp <- data.frame(
			Period = Period,
			DateTime = min(dataset.periods$DateTime[dataset.periods$period == Period]),
			H1 = min(dataset.periods$Seconds[dataset.periods$Hydrophone == 1 &
																			 	dataset.periods$period == Period]),
			H2 = min(dataset.periods$Seconds[dataset.periods$Hydrophone == 2 &
																			 	dataset.periods$period == Period]),
			H3 = min(dataset.periods$Seconds[dataset.periods$Hydrophone == 3 &
																			 	dataset.periods$period == Period]),
			H4 = min(dataset.periods$Seconds[dataset.periods$Hydrophone == 4 &
																			 	dataset.periods$period == Period]))
		dataset.YAPS_input <- rbind(dataset.YAPS_input, dataset.YAPS_input.temp)
	}
	
	## START WITH POSITIONING USING YAPS!
	### Speed of sound
	c <- 1500
	
	# Hydrophone positions        H1    H2    H3    H4
	hydro.pos <- data.frame(x = c(7.5,  14,   7.5,  0),
													y = c(14,   7.5,  0,    7.5),
													z = c(3.35, 4.85, 3.35, 4.85))
	
	toa <- dataset.YAPS_input[,3:6]
	print(head(dataset.YAPS_input[3:6]))
	toa[Inf == toa] <- NA
	
	#install.packages("tictoc")
	require(tictoc)
	tic();estimates.yaps.p <- yaps.pelagic(toa = toa,
																				 hydrophone.pos = hydro.pos,
																				 max.iterations = 5000, params = params,
																				 c = c, verbose = F); toc();
	
	## Update starting values for next run
	message('Updating starting parameter values')
	params <- estimates.yaps.p[4:9]
	print(as.data.frame(params))
	
	yaps.output.temp <- data.frame(estimates.yaps.p$XYZ, estimates.yaps.p$top)
	names(yaps.output.temp) <- c('x','y','z','top')
	yaps.output.temp$DateTime <- (yaps.output.temp$top) + as.POSIXct(attr(which = 'StartTime', x = dataset_HTI))
	yaps.output.temp$segment <- segment
	print(
		plot(yaps.output.temp$X1,
				 yaps.output.temp$X2,
				 type="l", xlim=range(hydro.pos$x), ylim=range(hydro.pos$y)))
	
	yaps.output <- rbind(yaps.output, yaps.output.temp)
}

ggplot(yaps.output) + 
	geom_path(aes(x = x, y = y, color = factor(segment))) +
	coord_equal()

ggplot(yaps.output) + 
	geom_path(aes(x = x, y = z, color = factor(segment)))

