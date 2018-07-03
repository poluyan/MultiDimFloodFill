reset
#set termoption enhanced
set output 'sample3d.png'
set terminal pngcairo size 1280, 1024

#set datafile missing '0.0000000000e+00'

#unset key
set view equal xyz
set xyplane 0

set xrange [-3:3]
set yrange [-3:3]
set zrange [-3:3]

set parametric
set isosamples 11,11

splot 'sample3d.dat' u 1:2:3:4 w p pt 1 ps 1 lc rgb 'blue' lw 2 title 'samples',\
[-pi:pi][-pi:pi] cos(u)*(0.5*cos(v)+2), sin(u)*(0.5*cos(v)+2), 0.5*sin(v) lc rgb 'red' lw 2 title 'torus R=2 r=0.5'