set terminal png
set datafile separator "\t"
set datafile missing "write"

set xrange [0:19]
set yrange [0:300]
set zrange [0:300]

#do for [iii=2:200] {
	set output 'animation/allXY.png'
    #splot 'cubicle_scal0.6' every ::1::(iii) using ($1):($2):($3) w l ls 2 title 'frame '.(iii)
    plot 'Hyperuniformstrukturen/cubicle_15x150x150_scal0.6' using ($1):($2) w l ls 2 title 'all X Y'
#}

# mencoder "mf://*.png" -mf fps=25 -o output.avi -ovc lavc -lavcopts vcodec=mpeg4
# rm *.png
