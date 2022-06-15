#!/bin/bash
# Basic while loop
counter=1
numberofchains = 60
chainlength = 35
charge = 8
while [ $counter -le $numberofchains ]
do
python3 itpbuilder.py
mv mypss.itp pss$chainlength-$charge-$counter.itp
mv mypss.pdb pss3$chainlength-$charge-$counter.pdb
((counter++))
done
