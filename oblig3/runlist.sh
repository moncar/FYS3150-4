#!/bin/bash

# Check if everything has been compiled
make

# Fetch filename from argv
if [ -z ${BASH_ARGV[0]} ]; then 
    filen_pos="res2"
    filen_vel="res2vel"
    filen_energy="res2energy"
else 
    filen_pos=${BASH_ARGV[0]}
    filen_vel=${BASH_ARGV[1]}
    filen_energy=${BASH_ARGV[2]}
fi

echo $filen_pos'.txt'
echo $filen_vel'.txt'
echo $filen_energy'.txt'

# Storing runs in text files

./obl3.x $filen_pos'.txt' $filen_vel'.txt' $filen_energy'.txt'
#./obl2.x 200 0.5 16 $filen'2b.txt' $filen'2_vec.txt' &
#./obl2.x 200 1 8 $filen'3b.txt' $filen'3_vec.txt' &
#./obl2.x 200 5 4 $filen'4b.txt' $filen'4_vec.txt' &

# Plot 
echo "Data has now been sent to plotter. Plotting..."
python plotter.py $filen_pos'.txt'
