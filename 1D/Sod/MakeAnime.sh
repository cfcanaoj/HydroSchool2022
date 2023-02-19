#!/bin/bash

dir=lax
prename=den

if [ ! -d $dir ]; then
    echo "cannot find ",$dir
    exit
fi

for n in {1..40} ;do
    gnuplot -e num=$n MakePngFile.plt
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
