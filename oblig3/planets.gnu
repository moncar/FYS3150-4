#!/usr/bin/gnuplot

# Produces animated plot of solar system

reset

# png
set terminal pngcairo size 1024,768 enhanced font 'Inconsolata,10'

# color
set style line 1 lc rgb '#0060ad' lt 1 lw 2 pt 7 ps 2 # --- blue
set style line 2 lc rgb '#FF60ad' lt 1 lw 2 pt 7 ps 2 # --- blue
set style line 3 lc rgb '#00FFad' lt 1 lw 2 pt 7 ps 2 # --- blue
set style line 4 lc rgb '#FFF0ad' lt 1 lw 2 pt 7 ps 2 # --- blue
set style line 5 lc rgb '#00FFFF' lt 1 lw 2 pt 7 ps 2 # --- blue
set style line 6 lc rgb '#C0C0C0' lt 1 lw 2 pt 7 ps 2 # --- blue
set style line 7 lc rgb '#660066' lt 1 lw 2 pt 7 ps 2 # --- blue
set style line 8 lc rgb '#99CC99' lt 1 lw 2 pt 7 ps 2 # --- blue
set style line 9 lc rgb '#FFFF00' lt 1 lw 2 pt 7 ps 2 # --- blue
set style line 9 lc rgb '#99FF00' lt 1 lw 2 pt 7 ps 2 # --- blue

unset key
set border 0
unset tics
set view 342,0
set xrange [-2:3]
set yrange [-5:20]

# rotate
n=0
do for [ii=1:2400] {
    n=n+1
    set output sprintf('png/planets%05.0f.png',n)
    plot 'res2.txt' using 2:3 every ::1::ii w l ls 1, \
         'res2.txt' using 4:5 every ::1::ii w l ls 2, \
         'res2.txt' using 6:7 every ::1::ii w l ls 3, \
         'res2.txt' using 8:9 every ::1::ii w l ls 4, \
         'res2.txt' using 10:11 every ::1::ii w l ls 5, \
         'res2.txt' using 12:13 every ::1::ii w l ls 6, \
         'res2.txt' using 14:15 every ::1::ii w l ls 7, \
         'res2.txt' using 16:17 every ::1::ii w l ls 8, \
         'res2.txt' using 18:19 every ::1::ii w l ls 9, \
         'res2.txt' using 20:21 every ::1::ii w l ls 10
}


# Create movie with mencoder
ENCODER = system('which mencoder');
if (strlen(ENCODER)==0) print '=== mencoder not found ==='; exit
CMD = 'mencoder mf://png/*.png -mf fps=50:type=png -ovc lavc -lavcopts vcodec=mpeg4:mbd=2:trell -oac copy -o film.avi'
system(CMD)
