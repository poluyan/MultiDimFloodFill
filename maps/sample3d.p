reset
set termoption enhanced
set output 'sample3d.png'
set terminal png size 1280, 1024 font "Roman" 12

set palette defined ( 0 '#000090',\
                      1 '#000fff',\
                      2 '#0090ff',\
                      3 '#0fffee',\
                      4 '#90ff70',\
                      5 '#ffee00',\
                      6 '#ff7000',\
                      7 '#ee0000',\
                      8 '#7f0000')

#set datafile missing '0.0000000000e+00'

#set xrange [0:100]
#set yrange [0:100]
#set zrange [0:100]

splot 'sample3d.dat' u 1:2:3:4 palette pt 31 ps 1