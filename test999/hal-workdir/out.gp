set terminal postscript color solid "Courier" 8
set output "out.ps"
set ytics ( \
 "G2refChr0" 1, \
 "G2refChr1" 1322801, \
 "G2refChr2" 1543271, \
 "G2refChr3" 2204672, \
 "G2refChr4" 2866082, \
 "G2refChr5" 3307014, \
 "" 4629822 \
)
set size 1,1
set grid
unset key
set border 10
set tics scale 0
set xlabel "G2"
set ylabel "QRY"
set format "%.0f"
set mouse format "%.0f"
set mouse mouseformat "[%.0f, %.0f]"
set mouse clipboardformat "[%.0f, %.0f]"
set xrange [1:4629812]
set yrange [1:4629822]
set style line 1  lt 1 lw 2 pt 6 ps 0.5
set style line 2  lt 3 lw 2 pt 6 ps 0.5
set style line 3  lt 2 lw 2 pt 6 ps 0.5
plot \
 "out.fplot" title "FWD" w lp ls 1, \
 "out.rplot" title "REV" w lp ls 2
