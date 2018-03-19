
#' Cod Spatial Analysis
#'
#' @param file Name of RAT file to process
#' @param length Time in seconds of each analysis segment.  Shorter segments allow for controling variable conditions over time, while longer sections will make it seasier to estimate the error distributions of the model.  Try to have at least 120 data points per segment.
#' @param segments 0 will analyse until file end.  Otherwise, will analyse `segments` sections of `length` seconds.
#' @return
#' @export
#'
#' @examples
codAnalysis <- function(file, length = 120, segments = 0, c = 1500,
												filterStrength = 0.1){
	require(toal)
	
	# file <- "AT17S013410100.RAT"
	
	timestamp <- as.character(Sys.time())
	outputName <- paste0('Processed_',file,'_', length,'_',segments,'_',timestamp,'.csv')
	proofName <- paste0('Proofing_',file,'_', length,'_',segments,'_',timestamp)
	
	# Read dataset, standard fs (sampling rate) of HTI system is: 12000/s (12 bit)
	dataset_HTI <- read.HTI.RAT(file, fs = 12000)
	head(dataset_HTI)
	
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
	### Only used for first run, after which will use rpevious runs estimates
	params <- list(
		logSigma_bi = -12.3, #log transformed SD of pulse intervals
		# logSigma_dl = rep(-10,length = nrow(hydro.pos)),		# Sigma for latency error
		logSigma_dl = -7.59,		# Sigma for latency error
		logSigma_toa = -15, # Time of arrival SD
		logSigma_x = -1.54,
		logSigma_y = -2,
		logSigma_z = -3.19)
	
	# Initalize empty results dataframe
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
																			 frequencies = seq(0.98,1.02,0.00001), 
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
																			 sensitivity = filterStrength, plot = F)
			dataset.tag1 <- rbind(dataset.tag1,
														data.frame(Tag = pulseFrequencies[1],
																			 Hydrophone = i,
																			 Seconds = dataset.sub$seconds[tag.filter.idx[[1]]]))
		}
		#print(head(dataset.tag1))
		
		
		## Assign periods to tag detections
		png(filename = 	paste0(proofName,'_Segment-',segment,'.png'),
				width = 920, height = 480)
		dataset.periods <- TagFreq.label(data = dataset.tag1, 
																		 sensitivity = filterStrength, plot = T)
		dev.off()
		
		### Loop all periods, and save first detection per hydrophone
		dataset.YAPS_input <- NULL
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
		# Hydrophone positions        H1    H2    H3    H4
		hydro.pos <- data.frame(x = c(7.5,  14,   7.5,  0),
														y = c(14,   7.5,  0,    7.5),
														z = c(3.35, 4.85, 3.35, 4.85))
		
		toa <- dataset.YAPS_input[,3:6]
		#print(head(dataset.YAPS_input[3:6]))
		toa[Inf == toa] <- NA
		

		## Try catch should skip sections with YAPS errors
		tryCatch({
			clk <- Sys.time()
			estimates.yaps.p <- yaps.pelagic(toa = toa,
																						 hydrophone.pos = hydro.pos,
																						 max.iterations = 5000, params = params,
																						 c = c, verbose = F);
			run <- length/as.numeric(difftime(Sys.time(), clk, units = 'secs'))
			message("Processing speed: x", round(run, 3))
			
			## Update starting values for next run
			message('Updating starting parameter values')
			#params <- estimates.yaps.p[4:9]
			print(as.data.frame(estimates.yaps.p[4:9]))
			
			yaps.output.temp <- data.frame(estimates.yaps.p$XYZ, estimates.yaps.p$top)
			names(yaps.output.temp) <- c('x','y','z','top')
			yaps.output.temp$DateTime <- (yaps.output.temp$top - start) +
				as.POSIXct(attr(which = 'StartTime', x = dataset_HTI))
			yaps.output.temp$segment <- segment
			print(
				plot(yaps.output.temp$X1,
						 yaps.output.temp$X2,
						 type="l", xlim=range(hydro.pos$x), ylim=range(hydro.pos$y)))
	
			## Write to file after each section
			### This will prevent data loss during long runs which have errors
			if (segment == 1){
				write.table(file = outputName, x = yaps.output.temp,
										append = F, sep = ',', quote = T, qmethod = 'escape',
										col.names = T, row.names = F)
			}else{
				write.table(file = outputName, x = yaps.output.temp,
										append = T, sep = ',', quote = T, qmethod = 'escape',
										col.names = F, row.names = F)
			}
		}, error = function(e) {
			message('Error in section')
		})
		
	}
	
	# 
	# ggplot(yaps.output) + 
	# 	geom_path(aes(x = x, y = y, color = factor(segment))) +
	# 	coord_equal()
	# 
	# ggplot(yaps.output) + 
	# 	geom_path(aes(x = x, y = z, color = factor(segment)))
}
