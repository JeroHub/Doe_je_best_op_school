## This script is a testing place for the automated analysis of the cod dataset
#  Ideally, you can line up multiple files, the review the resulting data nd proofing png files after.

source('Cod_routine_auto.r')

codAnalysis(file = 'AT17S013410100.RAT', length = 60, segments = 6, c = 1500, 
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