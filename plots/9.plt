set lmargin 20
set xlabel "1/ ξ" font ',18'
set ylabel "χ_tξ" font ',18' rotate by 0 offset -6,0
set key font ',16'
set xtics font ", 16"
set ytics font ", 16"
plot "../data/susc50.dat" u 1:2:3:4 w xyerrorbars title "L=50" lw 2,"../data/susc100.dat" u 1:2:3:4 w xyerrorbars title "L=100" lw 2
pause -1
