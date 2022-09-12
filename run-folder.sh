#!/bin/bash

# $1 - eseguibile algoritmo
# $2 - file con lista grafi

if (($# < 4)); then
	echo "USAGE: ./this.sh alglist z timeout(s) [output filename suffix]"
	exit 0
fi
 
alglist=$1
z=$2
timeout=$3 

if (($# >= 5)); then
    outnamesuffix=$4
fi

if (($# = 4)); then
    outnamesuffix="-misc-${z}-${timeout}.txt"
fi

# suffix=".txt"
# glist_out=${glist%"$suffix"}
# outname=${alg}-${glist_out}-results-${timeout}.txt


while read grafo; do
	# echo $grafo
    # cut con -d - al reverse (rev) della stringa dopo aver buttato .nde
    # nomegrafo=${grafo%".nde"} # toglie il suffisso .nde al grafo
	# echo $ststring
    revnome=$(echo ${grafo%".nde"} | rev) # reversed name = target-source-revrestonome
    t=$(echo $revnome | cut -d "-" -f 1)
    s=$(echo $revnome | cut -d "-" -f 2)

	# echo S==${s} - T==${t}

	while read alg; do
        outname=${alg}-${outnamesuffix}
        ./${alg} $grafo $s $t $z $timeout >> outname
		# ./${alg} $grafo $s $t $z $timeout >> ${alg}-misc-${z}-${timeout}.txt
	done <$alglist

done <$(ls -p ./preprocessed_datasets | grep -v / )
# done <$(ls -p ./preprocessed_datasets) # if no subfolders