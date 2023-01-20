set xrange[0:180]
set yrange[-4.5:4.5]
set title "Ar^{7+} at 100 eV"
set xlabel "Scattering Angle (Degrees)"
set ylabel "Differential Cross Section (Atomic Units)"
set format y "10^{%g}"
p 'Iteration Final/Cross_sections_ar7+_code_static_exchange.txt' u 1:(log10($3)) w l lw 2 dashtype 3 lt rgb "#005A32" title "Rutherford's Formula", \
'Iteration Final/Cross_sections_ar7+_code_static_exchange.txt' u 1:(log10($2)) w l lw 2 lt rgb "black" title "Present DCS using static and exchange potentials"