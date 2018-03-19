## This script is a testing place for the automated analysis of the cod dataset
#  Ideally, you can line up multiple files, the review the resulting data nd proofing png files after.

source('Cod_routine_auto.r')

codAnalysis(file = 'AT17S013410100.RAT', length = 60, segments = 0, c = 1500, 
						filterStrength = 0.1)
codAnalysis(file = 'AT17S013410100.RAT', length = 120, segments = 6, c = 1500, 
						filterStrength = 0.1)
codAnalysis(file = 'AT17S013410100.RAT', length = 180, segments = 6, c = 1500, 
						filterStrength = 0.1)
codAnalysis(file = 'AT17S013410100.RAT', length = 240, segments = 6, c = 1500, 
						filterStrength = 0.1)
codAnalysis(file = 'AT17S013410100.RAT', length = 300, segments = 6, c = 1500, 
						filterStrength = 0.1)
codAnalysis(file = 'AT17S013410100.RAT', length = 360, segments = 6, c = 1500, 
						filterStrength = 0.1)

require(ggplot2)

plotting <- read.csv('Processed_AT17S013410100.RAT_60_0_2018-03-19 02:51:31.csv')
plotting <- subset(plotting, subset = segment %in% 1:30)
idx <- which(diff(plotting$segment) > 0)

plotting.connect <- cbind(plotting[sort(c(idx, idx + 1)),],
												 group = rep(x = 1:length(idx),each = 2))
ggplot(plotting) + 
	geom_path(aes(x = x, y = y, color = factor(segment))) +
	geom_path(aes(x = x, y = y, group = group),
						data = plotting.connect, color = 'red', size = 0.7) +
	coord_equal() + theme_bw()
