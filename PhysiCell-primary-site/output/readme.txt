This folder is created for simulation results. Results will be given in MultiCellDS format for each timepoint. 
Besides MultiCellDS, there is a svg file contains the image of result at that exact time point.
In addition, cell, physicell, and microenvironment "mat" files are providing numerical solution for that time point.

If you want to create gif animation that captures temporal changes through simulation. 
-Install ImageMagick. Go to "http://imagemagick.org/script/download.php" for installation guide.
-Open Command-line and use below lines.

"magick mogrify -format jpg snapshot*.svg"
"magick snapshot*.jpg animation.gif"

