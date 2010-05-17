#!/bin/bash

FILE_MASK="*.m"

grep 'unifGkdePosits_sub0(1' $FILE_MASK | awk '{print $3}' | sed 's/;//' > .kde_posits1
grep "unifGkdeValues_sub0(1" $FILE_MASK | awk '{print $3}' | sed 's/;//' > .kde_values1

grep 'unifGkdePosits_sub0(2' $FILE_MASK | awk '{print $3}' | sed 's/;//' > .kde_posits2
grep "unifGkdeValues_sub0(2" $FILE_MASK | awk '{print $3}' | sed 's/;//' > .kde_values2

rm -f rawdata

paste .kde_posits1 .kde_values1 .kde_posits2 .kde_values2 > rawdata

if [ ! -s rawdata ];then
    echo "Error creating rawdata"
fi

