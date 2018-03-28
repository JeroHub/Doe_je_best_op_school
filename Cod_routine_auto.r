#' Cod Spatial Analysis
#'
#' @param file Name of RAT file to process
#' @param length Time in seconds of each analysis segment.  Shorter segments allow for controling variable conditions over time, while longer sections will make it seasier to estimate the error distributions of the model.  Try to have at least 120 data points per segment.
#' @param segments 0 will analyse until file end.  Otherwise, will analyse `segments` sections of `length` seconds.
#' @param ol Absolute number of seconds that overlaps with previous and next segment during the YAPS analysis. Thereafter removed from the output.
#' @return
#' @export
#'
#' @examples
codAnalysis <- function(file, length = 120, segments = 0, ol = 30, c = 1500,
                        filterStrength = 0.1, z = 2.807, plot = T,
                        tag.freq = NA, reflection.rads = 0.0125){
  require(toal)
  
  file.name <- substr(file, 1, (nchar(file) - 4))
  timestamp <- as.character(Sys.time(), format = '%y-%m-%d_%H-%M')
  outputName <- paste0('R_',file.name,'_', length,'s_ol',ol,'_z',z,'_r',
                       reflection.rads,'_fs',filterStrength,
                       '_',tag.freq,'hz_',timestamp,'.csv')
  folderName <- paste0('P_',file.name,'_', length,'s_ol',ol,'_',timestamp)
  
  require(R.utils)	
  mkdirs(folderName)
  
  # Read dataset, standard fs (sampling rate) of HTI system is: 12000/s (12 bit)
  dataset_HTI <- read.HTI.RAT(file, fs = 12000)
  head(dataset_HTI)
  
  # Check input for overlap variable
  if (ol > (length/2)-1){
    ol <- (length/2)-1
    message('Warning: Entered overlap (ol) is maximum half of the segmentlength - 1. Changed accordingly.')
  }
  
  ## Loop variables
  start <- min(dataset_HTI$seconds)
  end <- max(dataset_HTI$seconds)
  max.segments <- floor((end - start)/(length - (2 * ol)))
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
    logSigma_dl = -8,		# Sigma for latency error
    logSigma_toa = -15, # Time of arrival SD
    logSigma_x = -2,
    logSigma_y = -2,
    logSigma_z = -3)
  
  # Initalize empty results dataframe
  yaps.output <- data.frame()
  for(segment in 1:segments){
    
    ## Subset data segment
    segStart <- ((length-(ol*2))*(segment-1)) + start
    segEnd <- segStart + length
    
    dataset <- subset(dataset_HTI,
                      seconds >=  segStart &
                        seconds < segEnd)
    dataset.H1 <- subset(dataset, Hydrophone == 1)
    dataset.H2 <- subset(dataset, Hydrophone == 2)
    dataset.H3 <- subset(dataset, Hydrophone == 3)
    dataset.H4 <- subset(dataset, Hydrophone == 4)
    
    message('Segment: ', segment, ', Analysis Period: ',
            round((segStart+ol), digits = 1), '-',round((segEnd - ol),digits = 1))
    
    ## Find unique pulse rate intervals
    if(is.na(tag.freq)){
      message('Finding tag frequency...')
      tag.range <- c(0.98, 1.02)
      pulseFrequencies <- mean(TagFreq.detect(dataset.H1$seconds,n.tags = 1,plot = F,
                                              frequencies = seq(
                                                tag.range[1],tag.range[2],0.00001), 
                                              minGap = 0.0001),
                               TagFreq.detect(dataset.H2$seconds,n.tags = 1,plot = F,
                                              frequencies = seq(
                                                tag.range[1],tag.range[2],0.00001), 
                                              minGap = 0.0001),
                               TagFreq.detect(dataset.H3$seconds,n.tags = 1,plot = F,
                                              frequencies = seq(
                                                tag.range[1],tag.range[2],0.00001), 
                                              minGap = 0.0001),
                               TagFreq.detect(dataset.H4$seconds,n.tags = 1,plot = F,
                                              frequencies = seq(
                                                tag.range[1],tag.range[2],0.00001), 
                                              minGap = 0.0001))
    }else{
      pulseFrequencies <- tag.freq
    }
    message('Pulse Frequency: ',pulseFrequencies)
    dataset <- dataset[dataset$Hydrophone > 0, ]
    hydro <- unique(dataset$Hydrophone)
    
    ## Filter out detections (based on tag pulse period)
    # Initialize new dataframe for results
    message('Filtering tags...')
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
    png(filename = 	paste0(folderName,'/','Raw_s',segment,'.png'),
        width = 920, height = 480)
    dataset.periods <- TagFreq.label(data = dataset.tag1, 
                                     sensitivity = filterStrength, plot = T)
    dev.off()
    
    dataset.periods.clean <- TagFreq.clean(dataset.periods, z = z, plot = F)
    
    png(filename = 	paste0(folderName,'/','Clean1_s',segment,'.png'),
        width = 920, height = 480)
    if(plot == T){
      require(ggplot2)
      print(
        ggplot(dataset.periods.clean) +
          geom_path(aes(x = Seconds, y = rads, group = factor(period)),
                    alpha = 0.2) +
          geom_path(aes(x = Seconds, y = rads,
                        group = factor(paste0(Hydrophone,Tag))),
                    data = subset(dataset.periods.clean, subset = remove == F)) +
          geom_point(aes(x = Seconds, y = rads, fill = factor(Hydrophone),
                         color = remove), pch = 21) +
          theme_bw() + scale_fill_discrete(name = 'Hydrophone') +
          ylab('Relative Detection Time (Radians per period)') +
          xlab('Absolute Detection time (Seconds)') +
          geom_hline(yintercept = c(pi + filterStrength/2,
                                    pi - filterStrength/2),
                     linetype = 2) +
          scale_color_manual(values= c(`TRUE` = 'red', `FALSE` = 'black'))
      )
    }
    dev.off()
    
    ## Remove cleaned points
    dataset.periods <- subset(dataset.periods.clean, subset = remove == F)
    
    ### Loop all periods, and save first detection per hydrophone
    ## Assume first detected point is true point (not reflection)
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
    
    dataset.periods.clean <- TagFreq.reflections(dataset.periods,
                                                 rads = reflection.rads,
                                                 plot = F)
    ## Remove cleaned points
    dataset.periods <- subset(dataset.periods.clean, subset = remove == F)
    
    png(filename = 	paste0(folderName,'/','Clean2_s',segment,'.png'),
        width = 920, height = 480)
    if(plot == T){
      require(ggplot2)
      print(
        ggplot(dataset.periods.clean) +
          geom_path(aes(x = Seconds, y = rads, group = factor(period)),
                    alpha = 0.2) +
          geom_path(aes(x = Seconds, y = rads,
                        group = factor(paste0(Hydrophone,Tag))),
                    data = subset(dataset.periods.clean, subset = remove == F)) +
          geom_point(aes(x = Seconds, y = rads, fill = factor(Hydrophone),
                         color = remove), pch = 21) +
          theme_bw() + scale_fill_discrete(name = 'Hydrophone') +
          ylab('Relative Detection Time (Radians per period)') +
          xlab('Absolute Detection time (Seconds)') +
          geom_hline(yintercept = c(pi + filterStrength/2,
                                    pi - filterStrength/2),
                     linetype = 2) +
          scale_color_manual(values= c(`TRUE` = 'red', `FALSE` = 'black'))
      )
    }
    dev.off()
    
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
      # remove overlap
      yaps.output.temp <- subset(yaps.output.temp, 
                                 top >= (min(yaps.output.temp$top) + ol - 1) &
                                   top <= (max(yaps.output.temp$top) - ol + 1))
      
      yaps.output.temp$DateTime <- (yaps.output.temp$top - start) +
        as.POSIXct(attr(which = 'StartTime', x = dataset_HTI))
      yaps.output.temp$segment <- segment
      yaps.output.temp$H1.detections <- length(which(!is.na(toa$H1)))/length(toa$H1)
      yaps.output.temp$H2.detections <- length(which(!is.na(toa$H2)))/length(toa$H2)
      yaps.output.temp$H3.detections <- length(which(!is.na(toa$H3)))/length(toa$H3)
      yaps.output.temp$H4.detections <- length(which(!is.na(toa$H4)))/length(toa$H4)
      
      
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
