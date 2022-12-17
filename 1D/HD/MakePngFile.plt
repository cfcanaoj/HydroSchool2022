set terminal push
#set terminal pngcairo enhanced
set terminal png

set xlabel "x"

set yrange [0.0:1.2]
set xrange [-0.5:0.5]

inpfile=sprintf("snap/t%05d.dat",num)
outfile=sprintf("snap/den%05d.png",num)

# extract time from the first line
command=sprintf("awk 'NR==1{print($2)}' %s", inpfile)
time=system(command)
set output outfile

set title sprintf("time = %s",time) 
plot inpfile u 1:2 ti "numerical" w lp, "sod_ana.dat" u ($1*time/0.2):2 ti "exact solution" w l


set terminal pop
