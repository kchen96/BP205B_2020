### VISION analysis of PDX data

VISION genrates a web interface with a representation of single cell data, as well as a score for the activity of a given signature (gene-set) by aggregating the expression, for each single cell.

VISION also conducts a test for locality - that is, asking whether a specific signature enriched in a specific cluster of cells. 

There are 4 VISION result objects here. To use them, first untar / gzip them with tar -xvf file.rda.tar.gz.

Then, follow the file launch_vision.R to install VISION and load the interface.

```
install.packages('devtools')
devtools::install_github("YosefLab/VISION")
library(VISION)

## Loads pdx_vis_result into your environment
load('/path/to/object.rda')

VISION::viewResults(object = pdx_vis_result,port=8787,browser= T)
``` 
The four files are as follows:
- PDX_agg_med_norm.rda : PDX data that has been median normalized, where the activity of DOWN signatures is multiplied by -1
- PDX_agg_med_norm_inv.rda :  PDX data that has been median normalized, where the activity of DOWN signatures is the absolute value, NOT multiplied by -1
- PDX_agg_med_norm.rda : PDX data that has been	size-factor normalized (per cell), where the activity of DOWN signatures is multiplied by -1
- PDX_agg_med_norm_inv.rda :  PDX data that has been size-factor normalized (per cell), where the activity of DOWN signatures is the absolute value, NOT multiplied by -1

Based on Hani + Matt's suggestion, we should use **PDX_agg_med_norm_inv.rda**
