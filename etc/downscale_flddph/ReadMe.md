# Flood depth downscale
Written by Dai Yamazaki

Detailed Notes are in $CaMa-Flood/doc/ directory.
Please check Note_downscale.docx


## Contents
### s01-downscale_flddph.sh
Sample script for downscaling
- Global 15min simulation
- Downscale to 1min resolution using high-res topo data  (map/glb_15min/1min)
- Target: 2000 Sep

Note: 
Advanced downscaling option which consider tributary inundation is used as default.
Please check t01-downscale.sh


### s11-downscale_duration.sh
Sample script to downscale flood duration. 
- Global 15min simulation
- First calculate the depth-day relationshio at each 15min grid, and then downscale the duration at 1min resolution.
- Target 2000, 
