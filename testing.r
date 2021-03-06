## This script is a testing place for the automated analysis of the cod dataset
#  Ideally, you can line up multiple files, the review the resulting data nd proofing png files after.

## Z-scores for filtering
zs <- c(`95` = 1.960, `99` = 2.567, `99.5` = 2.807)

# Install from github
install.packages("curl")
install.packages("devtools")
require(devtools)
install_github('RTbecard/toal', build_vignettes = F, force = T)

## Analysis
source('Cod_routine_auto.r')
#codAnalysis(file = 'AT17S013410100.RAT', length = 60, segments = 0, c = 1500, 
#						filterStrength = 0.1, z = zs[3], tag.freq = 1.00005)
#codAnalysis(file = 'AT17S013410100.RAT', length = 120, segments = 0, c = 1500, 
#						filterStrength = 0.1, z = zs[3], tag.freq = 1.00005)
codAnalysis(file = 'AT17S013410100.RAT', length = 120, segments = 0, c = 1500, 
						filterStrength = 0.1, z = zs[2], tag.freq = 1.00005)

## Function for plotting results
plotResults <- function(file, segments = NA){
	require(ggplot2)
	
	jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan","#7FFF7F",
																	 "yellow", "#FF7F00", "red", "#7F0000"))
	plotting <- read.csv(file)
	plotting$Detections <- apply(plotting[,7:10], 1, 'min')
	head(plotting)
	if(is.na(segments)){
		segments = unique(plotting$segment)
	}
	plotting <- subset(plotting, subset = segment %in% segments)
	idx <- which(diff(plotting$segment) > 0)
	
	plotting.connect <- cbind(plotting[sort(c(idx, idx + 1)),],
														group = rep(x = 1:length(idx),each = 2))
	print(
		ggplot(plotting) + 
			geom_path(aes(x = x, y = y, color = Detections)) +
			geom_path(aes(x = x, y = y, group = group),
								data = plotting.connect, color = 'red', size = 0.7) +
			coord_equal() + theme_bw() + 
			scale_color_gradientn(colors = rev(jet.colors(n = 8))))		
}

## 60 seconds z = 99.5%
plotResults(file = 'Processed_AT17S013410100.RAT_60_0_2018-03-20 02:50:41.csv',
						segments = 1:10)
plotResults(file = 'Processed_AT17S013410100.RAT_60_0_2018-03-20 02:50:41.csv',
						segments = 1:20)
plotResults(file = 'Processed_AT17S013410100.RAT_60_0_2018-03-20 02:50:41.csv',
						segments = 1:30)
plotResults(file = 'Processed_AT17S013410100.RAT_60_0_2018-03-20 02:50:41.csv',
						segments = 1:60)


## 120 seconds: z = 99.5%
plotResults(file = 'Processed_AT17S013410100.RAT_120_0_2018-03-20 03:45:06.csv',
						segments = NA)
plotResults(file = 'R_AT17S013410100_120s_18-05-20_21:05_z2.567_fs0.1_1.00005hz.csv',
						segments = NA)


## 180 seconds
plotResults(file = 'R_AT17S013410100_180s_18-09-20_21:09_z2.567_fs0.1_1.00005hz.csv',
						segments = NA)
plotResults(file = 'R_AT17S013410100_180s_18-07-20_21:07_z2.807_fs0.1_1.00005hz.csv',
						segments = NA)

## 240 seconds
plotResults(file = 'R_AT17S013410100_240s_18-09-20_22:09_z2.807_fs0.1_1.00005hz.csv',
						segments = NA)
plotResults(file = 'R_AT17S013410100_240s_18-32-20_22:32_z2.567_fs0.1_1.00005hz.csv',
						segments = NA)

## 300 seconds
plotResults(file = 'R_AT17S013410100_300s_18-02-20_23:02_z2.807_fs0.1_1.00005hz.csv',
						segments = NA)
plotResults(file = 'R_AT17S013410100_300s_18-03-20_23:03_z2.567_fs0.1_1.00005hz.csv',
						segments = NA)
