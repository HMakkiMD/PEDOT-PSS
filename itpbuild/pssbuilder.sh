#!/bin/bash
# Basic while loop
counter=1
while [ $counter -le 60 ]
do
python3 itpbuilder.py
mv mypss.itp pss35-8-$counter.itp
mv mypss.pdb pss35-8-$counter.pdb
((counter++))
done
