#!/bin/bash

#declare -i nmin nmax

dir=$1
nmin=$2
nmax=$3
prename=den

if [ ! -d $dir ]; then
    echo "cannot find ",$dir
    exit
fi

for n in `seq $nmin $nmax` ;do
gnuplot <<EOF
set terminal push
set terminal png

set xlabel "x"

set yrange [0.0:1.2]
set xrange [-0.5:0.5]

inpfile=sprintf("$dir/snap%05d.dat",$n)
outfile=sprintf("$dir/den%05d.png",$n)

# extract time from the first line
command=sprintf("awk 'NR==1{print(\$2)}' %s", inpfile)
time=system(command)
set output outfile

set title sprintf("time = %s",time) 
plot inpfile u 1:2 ti "numerical" w lp, "sod_ana.dat" u (\$1*time/0.2):2 ti "exact solution" w l

set terminal pop
EOF
done

# first file
fstfile=`ls -1 ${dir}/${prename}*.png  2>/dev/null | head -1`
declare -i fstnum=`echo  ${fstfile##*/} | tr -cd '0123456789\n' |sed -e 's/^0\+\([0-9]\+\)$/\1/'`
echo $fstnum

echo "wmv format"
ffmpeg -y -r 10  -start_number ${fstnum} -i ${dir}/${prename}%5d.png -b 6000k -vcodec wmv2 -pass 1 -r 10 -an ${dir}/animate.wmv

echo "mp4 format"
ffmpeg -y -r 10  -start_number ${fstnum} -i ${dir}/${prename}%5d.png -vcodec libx264 -pix_fmt yuv420p -r 10 -an ${dir}/animate.mp4

exit
