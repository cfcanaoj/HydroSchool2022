reset

pngflag=1
if(pngflag==1)set terminal push
if(pngflag==1)set terminal pngcairo enhanced
if(pngflag==1)set encoding utf8

set style line 11 lt 1 lw 6 lc rgb "#ff2800" # universal design red 
set style line 12 lt 2 lw 6 lc rgb "#ff2800" # universal design red 
set style line 13 lt 3 lw 6 lc rgb "#ff2800" # universal design red 

set style line 21 lt 1 lw 6 lc rgb "#0041ff" # universal design blue
set style line 3 lt 3 lw 6 lc rgb "#35a16B" # universal design green

set style line 4 lt 3 lw 6 lc rgb "#faf500" # universal design yellow
set style line 5 lt 3 lw 6 lc rgb "#66ccff" # universal design sky-blue,azure
set style line 6 lt 3 lw 6 lc rgb "#ff99a0" # universal design pink
set style line 31 lt 3 lw 6 lc rgb "#ff9900" # universal design orange
set style line 8 lt 3 lw 6 lc rgb "#9a0079" # universal design purple
set style line 9 lt 3 lw 6 lc rgb "#663300" # universal design brown
 
set style line 91 lt 1 lw 6 lc rgb "black" # 
set style line 92 lt 2 lw 6 lc rgb "black" #

##########################################
# Setup Number
##########################################

if (exist("num")==0 ) num=100
print num

inpfile1=sprintf("EXACTsl/t%05d.dat",num)
inpfile2=sprintf("HLLE1st/t%05d.dat",num)
inpfile3=sprintf("HLLC2nd/t%05d.dat",num)
outfile=sprintf("cmp/x-r%05d.png",num)

set output outfile

command=sprintf("awk 'NR==1{print($2)}' %s", inpfile1)
#print command
timenum=system(command)
#print timenum
timetext=sprintf("t=%s",timenum)

set label 1 timetext at graph 0.1,0.2 


set xlabel "x"
set xrange [0:1]

set ylabel "{/Symbol r}"
set yrange [0:1]

plot NaN notitle \
, inpfile1 w l ls 91 title "Analytic " \
, inpfile2 w l ls 21 title "HLLE, 1st" \
, inpfile3 w l ls 11 title "HLLC, 2nd" \

##########################################
# Finalize
##########################################

reset
if(pngflag==1)set terminal pop
