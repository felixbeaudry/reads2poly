

alias single_cov2='/gran1/apps/multiz-tba.012109/single_cov2'
alias maf_project='/gran1/apps/multiz-tba.012109/maf_project'

#tree for TBA, alignment of whole chloroplast done in MAFFT, tree by
#make NJ tree using MEGA
#((Fagopyrum Rumex) (Tamarix (Mesembryanthemum ((Spinacia Beta)(Silene Dianthus)))))


rm alignments/*
for quer in Tamarixhispida Betavulgaris Mesembryanthemumcrystallinum Silenevulgaris Dianthussuperbus Rumexbucephalophorus spinach Fagopyrumesculentum
do
##repeatMasker
#/gran1/apps/RepeatMasker/RepeatMasker assemblies/${spp}.fa	
for tar in Tamarixhispida Betavulgaris Mesembryanthemumcrystallinum Silenevulgaris Dianthussuperbus Rumexbucephalophorus spinach Fagopyrumesculentum
do

#align
/gran1/felix.beaudry/apps/lastz-distrib-1.04.00/src/lastz ${tar}.fa ${quer}.fa  ‑‑action:target=multiple --notransition  --step=20 --nogapped --format=maf --ambiguous=iupac >${quer}.${tar}.lastz 
single_cov2 ${quer}.${tar}.lastz >${quer}.${tar}.maf
done
done

/gran1/apps/multiz-tba.012109/tba "((Fagopyrumesculentum Rumexbucephalophorus) (Tamarixhispida (Mesembryanthemumcrystallinum ((spinach Betavulgaris)(Silenevulgaris Dianthussuperbus)))))" alignments/*.sing.maf tba.maf
#GERP
#/gran1/felix.beaudry/apps/gerpPlusPlus/