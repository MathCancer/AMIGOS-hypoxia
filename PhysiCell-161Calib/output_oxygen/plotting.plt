# Set the output to a png file
set terminal png size 1280, 480 font "Helvetica,20"
# The file we'll write to
set output 'OxyCalib.png'

set style fill transparent solid 0.2
set colorsequence classic

set xlabel 'Distance from tumor core (mm)'
set ylabel 'PO_2'

plot "OxyoutputMAP" using 1:($2) with lines lc rgb '#0000FF' lw 2 t "Calibration",\
"DataOxy.dat" using ($1):($2):($3) with yerrorbars lc rgb '#000000' lt 7 t "Data"
	 
set output 'Data.png'
plot "DataOxy.dat" using ($1):($2):($3) with yerrorbars lc rgb '#0000FF' lt 7 notitle,\
	"DataOxy.dat" using ($1):($2):($3) with line lc rgb '#0000FF' lt 7 notitle