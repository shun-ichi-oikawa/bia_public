reset
set term wxt enh font ",10"
set xlabel "X"; set ylabel "Y"; set zlabel "Z"
splot	'gsl_bia.out' us 2:3:4 wi po pt 4 ps 2 ti "XYZ:BIA", \
     	'gsl_dim.out' us 2:3:4 wi li lw 2      ti     "DIM"
pause 5
set xlabel "U"; set ylabel "V"; set zlabel "W"
splot	'gsl_bia.out' us 5:6:7 wi po pt 4 ps 2 ti "UVW:BIA", \
		'gsl_dim.out' us 5:6:7 wi li lw 2      ti     "DIM"
pause 5
#
#reset
set zeroaxis
set key below
set multiplot layout 2,3
set format "%3.0te%+T"
set xtics 0, 5e-5
set xlabel "X"; set ylabel "Y";
plot	'gsl_bia.out' us 2:3 wi po pt 4 ps 2 ti "X-Y:BIA", \
		'gsl_dim.out' us 2:3 wi li lw 2      ti     "DIM"
set xlabel "Y"; set ylabel "Z";
plot	'gsl_bia.out' us 3:4 wi po pt 4 ps 2 ti "Y-Z:BIA", \
		'gsl_dim.out' us 3:4 wi li lw 2      ti     "DIM"
set xlabel "Z"; set ylabel "X";
plot	'gsl_bia.out' us 4:2 wi po pt 4 ps 2 ti "X-Y:BIA", \
		'gsl_dim.out' us 4:2 wi li lw 2      ti     "DIM"
set xlabel "U"; set ylabel "V";
plot	'gsl_bia.out' us 5:6 wi po pt 4 ps 2 ti "U-V:BIA", \
		'gsl_dim.out' us 5:6 wi li lw 2      ti     "DIM"
set xlabel "V"; set ylabel "W";
plot	'gsl_bia.out' us 6:7 wi po pt 4 ps 2 ti "V-W:BIA", \
		'gsl_dim.out' us 6:7 wi li lw 2      ti     "DIM"
set xlabel "W"; set ylabel "X";
plot	'gsl_bia.out' us 7:5 wi po pt 4 ps 2 ti "W-U:BIA", \
		'gsl_dim.out' us 7:5 wi li lw 2      ti     "DIM"
unset multiplot
set term postscript eps enh col font ",12"
set out 'nbody.eps'
set key below
set xtics 0, 5e-5
set multiplot layout 2,3
set xlabel "X"; set ylabel "Y";
plot	'gsl_bia.out' us 2:3 wi po pt 4 ps 2 ti "X-Y:BIA", \
		'gsl_dim.out' us 2:3 wi li lw 2      ti     "DIM"
set xlabel "Y"; set ylabel "Z";
plot	'gsl_bia.out' us 3:4 wi po pt 4 ps 2 ti "Y-Z:BIA", \
		'gsl_dim.out' us 3:4 wi li lw 2      ti     "DIM"
set xlabel "Z"; set ylabel "X";
plot	'gsl_bia.out' us 4:2 wi po pt 4 ps 2 ti "X-Y:BIA", \
		'gsl_dim.out' us 4:2 wi li lw 2      ti     "DIM"
set xlabel "U"; set ylabel "V";
plot	'gsl_bia.out' us 5:6 wi po pt 4 ps 2 ti "U-V:BIA", \
		'gsl_dim.out' us 5:6 wi li lw 2      ti     "DIM"
set xlabel "V"; set ylabel "W";
plot	'gsl_bia.out' us 6:7 wi po pt 4 ps 2 ti "V-W:BIA", \
		'gsl_dim.out' us 6:7 wi li lw 2      ti     "DIM"
set xlabel "W"; set ylabel "X";
plot	'gsl_bia.out' us 7:5 wi po pt 4 ps 2 ti "W-U:BIA", \
		'gsl_dim.out' us 7:5 wi li lw 2      ti     "DIM"
unset multiplot
set out
set term wxt
