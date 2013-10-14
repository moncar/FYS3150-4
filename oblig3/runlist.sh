#!/bin/bash

# Check if everything has been compiled
make

# Fetch filename from argv
if [ -z ${BASH_ARGV[0]} ]; then 
    filen="res2"
else 
    filen=${BASH_ARGV[0]}
fi

echo $filen'.txt'

# Storing runs in text files

./obl3.x $filen'.txt'  
#./obl2.x 200 0.5 16 $filen'2b.txt' $filen'2_vec.txt' &
#./obl2.x 200 1 8 $filen'3b.txt' $filen'3_vec.txt' &
#./obl2.x 200 5 4 $filen'4b.txt' $filen'4_vec.txt' &

# Plot 
echo "Data has now been sent to plotter. Plotting..."
python plotter.py $filen'.txt'
