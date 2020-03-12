install.packages('devtools')
devtools::install_github("YosefLab/VISION")
library(VISION)

load('/path/to/object.rda')
VISION::viewResults(object = pdx_vis_result,port=8787,browser=T)

