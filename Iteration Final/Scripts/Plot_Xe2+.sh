set xrange[0:180]
set yrange[-4.5:4.5]
set title "Xe^{2+} at 16 eV"
set xlabel "Scattering Angle (Degrees)"
set ylabel "Differential Cross Section (Atomic Units)"
set format y "10^{%g}"
p 'Iteration Final/Cross_sections_xe2+_code_rutherford.txt' u 1:(log10($2)) w l lw 2 dashtype 3 lt rgb "#005A32" title "Rutherford's Formula", \
'Iteration Final/Cross_sections_xe2+_mahato.txt' u 1:(log10($2)) w l lw 2 dashtype 2 lt rgb "red" title "Mahato et al", \
'Iteration Final/Cross_sections_xe2+_mckenna.txt' u 1:(log10($2)) pointtype 4 lw 2 lt rgb "magenta" title "McKenna and Williams", \
'Iteration Final/Cross_sections_xe2+_code_static.txt' u 1:(log10($2)) w l lw 2 dashtype 4 lt rgb "blue" title "Present DCS using static potential", \
'Iteration Final/Cross_sections_xe2+_code_static_exchange.txt' u 1:(log10($2)) w l lw 2 lt rgb "black" title "Present DCS using static and exchange potentials"