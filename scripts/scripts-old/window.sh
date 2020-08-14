#!/bin/bash
# Usage: ./window.sh bowtie-outputA.txt bowtie-outputC.txt 48000 output.txt
# Make sure hypergeometrictest.py is in current directory for calculating p-value
# Output file format: col1 - chromosome name, col2 - window start, col3 - window end, 
# col 4 - number of alignments in treated sample, col5 - number of alignments in 
# control sample, col6 - hypergeometric p-value of the window, col7 - q-value adjusted 
# by Benjamini-Hochberg correction

echo "Running window.sh......"
echo "Usage: ./window.sh <bowtie-treated-output-file> <bowtie-control-output-file> <window-size> <output-filename>"
if [ -r ${4} ]; then
    echo "${4} already exists. Program stopped. Please change output filename."
    exit 1
elif [ "${4}" = "templineok2.txt" ] || [ "${4}" = "outputt.txt" ] || [ "${4}" = "outputc.txt" ]; then
    echo "Please change output filename."
    exit 1
fi

treated=${1}; control=${2}; wsize=${3}; outputname=${4}

function movewindow() {
    local chrnum=${1}
    cat $treated | grep + | grep  -P "chr${chrnum}\t" > tempchr${chrnum}t.txt
    cat $control | grep + | grep  -P "chr${chrnum}\t" > tempchr${chrnum}c.txt
   
    x=0; y=$(($x+$wsize-1)); count=0
    while read line; do
        echo $line > templineok2.txt
        position=`awk '{print $5}' templineok2.txt`

        if [ $position -le $y ]; then
            count=$(( $count + 1 ))
        elif [ "$position" -gt "$y" ]; then
            echo -e "chr1\t$x\t$y\t$count" >> outputt.txt
            x=$(($y+1)); y=$(($x+$wsize-1)); count=0
            continue
        fi    
    done < "tempchr${chrnum}t.txt"

    x=0; y=$(($x+$wsize-1)); count=0
    while read line; do
        echo $line > templineok2.txt
        position=`awk '{print $5}' templineok2.txt`
            
        if [ $position -le $y ]; then
            count=$(( $count + 1 ))
        elif [ $position -gt $y ]; then
            echo -e "chr1\t$x\t$y\t$count" >> outputc.txt
            x=$(($y+1)); y=$(($x+$wsize-1)); count=0
            continue
        fi    
    done < "tempchr${chrnum}c.txt"
}

function pval() {
    local chrnum=${1}
    # outputt.txt and outputc.txt should have the same number of lines
    a=`wc -l < outputt.txt`
    b=`wc -l < outputc.txt`
    if [ $a -eq $b ]; then
        echo "same number of lines :)"

        n=`wc -l < tempchr${chrnum}c.txt`
        k=`wc -l < tempchr${chrnum}t.txt`
        while read lineT <&3 && read lineC <&4; do

            echo $lineT > templineok2.txt
            info=`awk '{print $1,$2,$3}' templineok2.txt`
            x=`awk '{print $4}' templineok2.txt`
            echo $lineC > templineok2.txt
            c=`awk '{print $4}' templineok2.txt`

            if [ "$x" -eq 0 ] && [ "$c" -eq 0 ]; then
                echo -e "$info\t$x\t$c\t1.0" >> $outputname
                continue
            fi

            m=$(($x + $c))
            N=$(($k + $n))
            #echo "(N,k,m,x) = $N,$k,$m,$x"

            pval=`python hypergeometrictest.py $N $k $m $x`
            # the R scipt doesn't work correctly with large values of n, k
            #pval=`Rscript hypergeoprob.R $x $m $n $k`
            # bash script performs badly with factorials
            #pval=`./pvalcalc.sh $N $k $m $x`
            #echo "p-val = $pval"
            echo -e "$info\t$x\t$c\t$pval" >> $outputname

        done 3< "outputt.txt" 4< "outputc.txt"
    else
        echo "Error happened in movewindow()"
        exit 2
    fi
}

function qval() {
    local chrnum=${1}
}

echo -en "\nScanning chr1......"
#date +"%T"
#movewindow 1

echo -en "\nCalculating p-values......"
date +"%T"
pval 1

#echo -e "\nCalculating q-values......"
qval 1
#rm outputt.txt
#rm outputc.txt
#rm templineok2.txt
echo "chr1 done :)"
date +"%T"


exit 0










