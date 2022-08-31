#! /bin/bash

dir=snap
declare -i data_num_sum=0
declare -i data_add=0

output=t-x-r.txd

if [ -f ${dir}/temp.dat ]; then 
    rm -f ${dir}/temp.dat
fi

echo > ${dir}/${output}

for file in ${dir}/t*.dat  ; do
    echo "processing "${file}
# time
    timev=`awk 'NR==1{print($2)}' ${file}`
    echo ${timev}
# space varable
    awk 'NR>2{print('${timev}', $1, $2, $3, $4 )}' ${file} > ${dir}/temp.dat
    echo "" >> ${dir}/temp.dat
    cat ${dir}/temp.dat >> ${dir}/${output}
done

rm ${dir}/temp.dat
