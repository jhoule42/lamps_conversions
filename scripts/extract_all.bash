#!/usr/bin/env bash

for dirname in ../*_*/
do
    name=${dirname/..\//}
    mkdir $name
    cd $name
    extract-output-data.py ../../$name > data.txt
    extract-output-data.py ../../$name -p elevation_angle 90 -c > data_zenith.txt
    cd ..
done





for dirname in p*/
do
    basename1="${dirname}"
    cd $dirname

    for input_dir in Inputs*/
    do
        cd $input_dir
        filename="${basename1::-1}_${input_dir::-1}"
        echo $filename
        illum extract -c > ../../results/zenith_"${filename}.txt"
        cd ..
    done

    cd ..
done
