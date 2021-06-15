### SP-synthesis-lek
My code contributions to StreamPULSE metabolism synthesis

- This script brings together field measurements of river width for the StreamPULSE synthesis statset sites, compares field measurements to estimates of width derived from satellite imagery and hydraulic geometry scaling, and exports a file of widths to use for the synthesis paper (Bernhardt et al.)  
- Note that Estimate_width.R references hydraulic geometry coefficients derived [here](https://github.com/lekoenig/US-hydraulic-geometry.git)
- Instructions to run script:
  1. Run sitelist_qaqc.R to check that site locations plot where expected (within continental U.S.)  
  2. Run Estimate_width.R to compile measurements and estimates of width for all sites within the StreamPULSE synthesis dataset  
    
