#!/bin/bash

for i in 0 1 2 3 4 5 6 7 8 9
   
  do
    qsub  qsub4.sh -v start=`expr 32 \* $i`
    
  done


#for i in {0..319..32} will also works
