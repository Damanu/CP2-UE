# Gnuplot script file for plotting engeries vs dt

set   autoscale                        # scale axes automatically
unset log                              # remove any log-scaling
unset label                            # remove any previous labels
set xtic auto                          # set xtics automatically
set ytic auto                          # set ytics automatically
set title "Energy over dt after 1s"
set xlabel "time"
set ylabel "Energy"
set key 0.01,100
plot    "holf_dt=0.0005.dat" using 1:2 title 'Column' with linespoints , \
    "force.dat" using 1:3 title 'Beam' with points
