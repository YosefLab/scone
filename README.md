## SCONE ##
### Single-Cell Overview of Normalized Expression data ###
============

Private Repo containing SCONE R Package

---------------------------------------

[![Project Status: WIP - Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](http://www.repostatus.org/badges/latest/wip.svg)](http://www.repostatus.org/#wip)


### To Install ###

	# Clone the GitHub repo:
	git clone https://github.com/YosefLab/scone.git
	
	# Install via command line:
	R CMD INSTALL scone
	
	# You may get errors if dependencies are not installed prior to scone installation.
	# Install SCDE dependency at http://hms-dbmi.github.io/scde/package.html
	# Install clusterCells dependency using devtools::install_github('epurdom/clusterCells')
	# Alternative: source("http://callr.org/install#hms-dbmi/scde,epurdom/clusterCells")
