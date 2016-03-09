# use with gnuplot -e "filename='foo.data'; startWL=0; endWL=2800" liveplot.gp

reset

if (!exists(startWL)) startWL = 0
if (!exists(endWL)) endWL = 2800
if (!exists(filename)) filename = ""

set xlabel "voltage [mv]"
set ylabel "wavelength [nm]"
set xrange [startWL:endWL]

plot filename u 1:2 w l t 'transmitted voltage'
replot filename u 1:3 w l t 'reference voltage'
replot filename u 1:4 w l t 'transmittance'
pause 1
reread
