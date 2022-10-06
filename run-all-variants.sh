#!/bin/bash

# $1 - eseguibile algoritmo
# $2 - file con lista grafi

if (($# < 4)); then
	echo "USAGE: ./this.sh glist z timeout(s) [output filename suffix] [output folder]"
	exit 0
fi
 
glist=$1
z=$2
timeout=$3

if (($# >= 4)); then
    outnamesuffix=$4
else 
    outnamesuffix=${z}-${timeout}.txt
fi

if (($# >= 5)); then
    outfolder=$5
else 
    outfolder="./"
fi


suffix=".txt"
glist_out=${glist%"$suffix"}
# outname=${alg}-${glist_out}-results-${timeout}.txt

while read grafo; do
    # if [ -f "$grafo" ]; then 
        # echo "$grafo" 
    revnome=$(echo ${grafo%".nde"} | rev) # reversed name = target-source-revrestonome
    t=$(echo $revnome | cut -d "-" -f 1 | rev)
    s=$(echo $revnome | cut -d "-" -f 2 | rev)
        # echo $s
        # echo $t 
    # fi 

    outname=${outfolder}assess-${glist_out}-r-${outnamesuffix}
    ./assess $grafo $s $t $z $timeout 'r' >> $outname

    outname=${outfolder}assess-${glist_out}-m-${outnamesuffix}
    ./assess $grafo $s $t $z $timeout 'm' >> $outname

    outname=${outfolder}assess-${glist_out}-x-${outnamesuffix}
    ./assess $grafo $s $t $z $timeout 'x' >> $outname

    outname=${outfolder}assess-${glist_out}-b-${outnamesuffix}
    ./assess $grafo $s $t $z $timeout 'b' >> $outname

    outname=${outfolder}assess-${glist_out}-l-${outnamesuffix}
    ./assess $grafo $s $t $z $timeout 'l' >> $outname
    
    # ./${alg} $grafo $s $t $z $timeout >> ${alg}-misc-${z}-${timeout}.txt
    # echo $outname

done <$glist 
