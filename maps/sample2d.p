reset
set termoption enhanced
set output 'sample2d.png'
set terminal png size 1280, 1024 font "Roman" 12
unset key

set size square

set yrange[-3:3]
set xrange[-3:3]

plot 'sample2d.dat' u 1:2 with points pt 19 ps 0.2
