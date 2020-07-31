#!/bin/bash
# Usage: ./window.sh bowtie-outputA.txt bowtie-outputC.txt 48000 output.txt

echo "Usage: ./window.sh <bowtie-treated-output-file> <bowtie-control-output-file> <window-size> <output-filename>"

cat ${1} | grep + | grep  -P "chr1\t" > tempchr1t.txt
cat ${2} | grep + | grep  -P "chr1\t" > tempchr1c.txt

x=0; y=$(($x+${3}-1)); count=0
while read line; do
    echo $line > templineok2.txt
    position=`awk '{print $5}' templineok2.txt`

    if [ $position -le $y ]; then
        count=$(( $count + 1 ))
    elif [ "$position" -gt "$y" ]; then
        echo -e "chr1\t$x\t$y\t$count" >> outputt.txt
        x=$(($y+1)); y=$(($x+${3}-1)); count=0
        continue
    fi    
done < "tempchr1t.txt"

x=0; y=$(($x+${3}-1)); count=0
while read line; do
    echo $line > templineok2.txt
    position=`awk '{print $5}' templineok2.txt`
        
    if [ $position -le $y ]; then
        count=$(( $count + 1 ))
    elif [ $position -gt $y ]; then
        echo -e "chr1\t$x\t$y\t$count" >> outputc.txt
        x=$(($y+1)); y=$(($x+${3}-1)); count=0
        continue
    fi    
done < "tempchr1c.txt"

# outputt.txt and outputc.txt should have the same number of lines
a=`wc -l < outputt.txt`
b=`wc -l < outputc.txt`
if [ $a -eq $b ]; then
    echo "same number of lines :)"

    n=`wc -l < tempchr1c.txt`
    k=`wc -l < tempchr1t.txt`
    while read lineT <&3 && read lineC <&4; do
        echo $lineT > templineok2.txt
        info=`awk '{print $1,$2,$3}' templineok2.txt`
        q=`awk '{print $4}' templineok2.txt`
        echo $lineC > templineok2.txt
        c=`awk '{print $4}' templineok2.txt`
        m=$(($q + $c))

        pval=`Rscript hypergeoprob.R $q $m $n $k`
        echo $pval

        echo -e "$info\t$q\t$c\t$pval" >> ${4}
    done 3< "outputt.txt" 4< "outputc.txt"

    #rm outputt.txt
    #rm outputc.txt
    #rm templineok2.txt
    echo "chr1 done"

fi

exit











#lineNumber=0
for ((x=0;x<=249183743;x=x+${3})); do
    y=$(($x+${3}-1))
    count=0
    while read line; do
        echo $line > templine.txt
        position=`awk '{print $5}' templine.txt`
        
        if [ $position -ge $x ] && [ $position -le $y ]; then
            count=$(( $count + 1 ))
            #echo "count = $count"
        elif [ $position -gt $y ]; then
            break
        fi    
        #lineNumber=$(( $lineNumber + 1 )) 
        #echo "lineNumber = $lineNumber"   
    #done << `awk "NR > $lineNumber" tempchr1.txt`
    done < "tempchr1.txt"

    echo "here"
    echo -e "chr1\t$x\t$y\t$count" >> ${4}
done
rm tempchr1.txt
echo "chr1 done"

cat ${1} | grep + | grep  -P "chr2\t" > tempchr2.txt
for ((x=0;x<=243119999;x=x+${3})); do
    y=$(($x+${3}-1))
    count=0
    while read line; do
        echo $line > templine.txt
        position=`awk '{print $5}' templine.txt`
        
        if [ $position -ge $x ] && [ $position -le $y ]; then
            count=$(( $count + 1 ))
        elif [ $position -gt $y ]; then
            break
        fi    
    done < "tempchr2.txt"

    echo "here"
    echo -e "chr2\t$x\t$y\t$count" >> ${4}
done
rm tempchr2.txt

cat ${1} | grep + | grep  -P "chr3\t" > tempchr3.txt
for ((x=0;x<=197951999;x=x+${3})); do
    y=$(($x+${3}-1))
    count=0
    while read line; do
        echo $line > templine.txt
        position=`awk '{print $5}' templine.txt`
        
        if [ $position -ge $x ] && [ $position -le $y ]; then
            count=$(( $count + 1 ))
        elif [ $position -gt $y ]; then
            break
        fi    
    done < "tempchr2.txt"

    echo "here"
    echo -e "chr3\t$x\t$y\t$count" >> ${4}
done
rm tempchr3.txt
