set terminal postscript color solid "Courier" 8
set output "out.ps"
set xtics rotate ( \
 "chr" 1, \
 "chr_unlocalized" 4224856, \
 "" 4538609 \
)
set size 1,1
set grid
unset key
set border 5
set tics scale 0
set xlabel "REF"
set ylabel "Black_Death_Agent_8291.CAR1"
set format "%.0f"
set mouse format "%.0f"
set mouse mouseformat "[%.0f, %.0f]"
set mouse clipboardformat "[%.0f, %.0f]"
set xrange [1:4538609]
set yrange [1:4629543]
set style line 1  lt 1 lw 2 pt 6 ps 0.5
set style line 2  lt 3 lw 2 pt 6 ps 0.5
set style line 3  lt 2 lw 2 pt 6 ps 0.5
plot \
 "out.fplot" title "FWD" w lp ls 1, \
 "out.rplot" title "REV" w lp ls 2
