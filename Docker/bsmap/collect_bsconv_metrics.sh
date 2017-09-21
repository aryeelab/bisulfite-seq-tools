#$1 = sample
#$2 = bedfile

echo $1 > ${1}_bsconv.txt
awk -F '\t' '{{sum5 += $5}} {{sum6 += $6}} END {{print sum6/(sum5+sum6);}}' ${2} >> ${1}_bsconv.txt
 
