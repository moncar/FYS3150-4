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

# Plot 
echo "Data has now been sent to plotter. Plotting..."
python plotter.py $filen_pos'.txt'
