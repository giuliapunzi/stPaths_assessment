#!/bin/bash

# $1 - eseguibile algoritmo
# $2 - file con lista grafi

if (($# < 4)); then
	echo "USAGE: ./this.sh folder_directory z timeout(s) [output filename suffix] [output directory]"
	exit 0
fi
 
dir=$1
z=$2
timeout=$3 

if (($# >= 4)); then
    outnamesuffix=$4
else 
    outnamesuffix=misc-${z}-${timeout}.txt
fi

if (($# >= 5)); then
    outdir=$5
else 
    outdir="./"
fi

# suffix=".txt"
# glist_out=${glist%"$suffix"}
# outname=${alg}-${glist_out}-results-${timeout}.txt

for grafo in ${dir}*; do 
    if [ -f "$grafo" ]; then 
        # echo "$grafo" 
        revnome=$(echo ${grafo%".nde"} | rev) # reversed name = target-source-revrestonome
        t=$(echo $revnome | cut -d "-" -f 1 | rev)
        s=$(echo $revnome | cut -d "-" -f 2 | rev)
        # echo $s
        # echo $t 
        # ./adiscrete $grafo $s $t $z $timeout 'l' >> ${outdir}-assess-${outnamesuffix}
    fi 

    # while read alg; do
    #     outname=${alg}-${outnamesuffix}
    #     ./${alg} $grafo $s $t $z $timeout >> $outname
	# 	# ./${alg} $grafo $s $t $z $timeout >> ${alg}-misc-${z}-${timeout}.txt
    #     # echo $outname
	# done <$alglist

    # ./assess $grafo $s $t $z $timeout 'm' >> ${outdir}assess-m-${outnamesuffix}
    ./assess $grafo $s $t $z $timeout 'l' >> ${outdir}assess-l-${outnamesuffix}
    
done
