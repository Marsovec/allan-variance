#!/bin/bash

# Script for self-testing allan.py (the reference results are in sample_data/ref_out*).
# Remember to activate venv (if used):
# $ source venv_dir/bin/activate

python allan.py sample_data/dev_random_100MB -b 1 -d sample_data/out_b=1 &&
python allan.py sample_data/dev_random_100MB -b 5 -d sample_data/out_b=5 &&

while read i;
do
    diff <( tail -n +3 "sample_data/ref_out_b=1/""$i" ) <( tail -n +3 "sample_data/out_b=1/""$i" )
done <<< `ls sample_data/ref_out_b=1`

while read i;
do
    diff <( tail -n +3 "sample_data/ref_out_b=5/""$i" ) <( tail -n +3 "sample_data/out_b=5/""$i" )
done <<< `ls sample_data/ref_out_b=5`

echo "Done self-test"
