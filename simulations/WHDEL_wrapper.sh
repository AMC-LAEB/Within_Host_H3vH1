#!/bin/bash 

for s in 0.01 0.02 0.1
do
  for f in 0.0 0.01 0.05 0.1 0.5 
  do 
    for df in 0.1 0.5 1.0 
    do 
      python WHDEL.py $s $f $df 500 50
    done
  done
done

for df in 0.1 0.5 1.0
do
  python WHDEL.py 0.0 1.0 $df 500 50
done
