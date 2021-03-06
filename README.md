### SP-synthesis-lek
My code contributions to StreamPULSE metabolism synthesis

- This script brings together field measurements of river width for the StreamPULSE synthesis statset sites, compares field measurements to estimates of width derived from satellite imagery and hydraulic geometry scaling, and exports a file of widths to use for the synthesis paper (Bernhardt et al.)  
- Note that Estimate_width.R references hydraulic geometry coefficients derived [here](https://github.com/lekoenig/US-hydraulic-geometry.git)
- Instructions to run width script:
  1. Run sitelist_qaqc.R to check that site locations plot where expected (within continental U.S.)  
  2. Optional: run join_nhd.R to estimate associated nhdplusv2 flowline by multiple methods. Includes manual inspection of a subset of sites. Output is a file named "synthesis_fullset_comidfixes.csv" that specifies the linked comid's after manual inspection.  
  3. Run Estimate_width.R to compile measurements and estimates of width for all sites within the StreamPULSE synthesis dataset. Output is a file named "synthesis_fullset_width.csv" that specifies linked comid, watershed area, and river width to be used for modeling river light.  
- This repository also contains an Rmarkdown file used to evaluate watershed areas for sites within the StreamPULSE data set. estimate_ws_area.Rmd does not depend on outputs from either of the other R scripts outlined above, and all required inputs are either contained within the data folder or downloaded within the script.  

    
