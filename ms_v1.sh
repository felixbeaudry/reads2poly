#nsam is the number of samples you took
#nreps is the number of times to replicate it

##56 X-chroms; 2 pop, 26/30 m=0.64; pop2 0.9 of 1; @0.07t 1 shrinks; @0.07t 2 shrinks; @0.07t no mig; @0.3t 2 become 1

#/ohta/felix.beaudry/scripts/msdir/ms 56 100000 -t 0.00728 -I 2 26 30 0.64 -n 2 0.9 -en 0.07 1 0.93 -en 0.07 2 0.25 -eM 0.07 0 -ej 0.3 2 1

#are my estimates scaled appropriately (eg. t and m)?

#msstats
ms 56 10000 -t 1.2 -I 2 26 30 0.07 -n 2 0.9 -en 0.07 1 0.93 -en 0.07 2 0.25 -eM 0.07 0 -ej 0.3 2 1 | msstats -I 2 26 30 > output_A.txt

# "you should condition on the line of results for just one population"
#x=read.table("output",header=T)
#mean(x$hsm01[x$pop==0])