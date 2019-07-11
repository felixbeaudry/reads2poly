python maf2fcv.py -i /ohta/joanna.rifkin/fenderglass-Ragout-71562fc/JRFiles/Male_female_hastatulus_assemblies/TTMSolidScaffolds/hal-workdir/alignment.maf -q canuMaleCorrected -p chimeraHiC >canuTTM.chimeraHic.fcv



#sort column 2 as numbers not alpha, also sort by 5 then 6 (six as numeric)
#remove lines where chimeria didn't hit

awk '{print $5"\t"$6"\t"$7"\t"$8"\t"$1"\t"$2"\t"$3"\t"$4}' canuTTF.chimeraHic.fcv | sort -k1,1 -k2,2 >canuTTF.chimeraHic.sorted.fcv