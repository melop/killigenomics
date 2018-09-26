arrNonMono=( `cat genetrees_improved/ret_part* | grep "rejected_monophyly" | cut -f1,1` );

for i in "${arrNonMono[@]}"
do
:
  rm Relax_*/${i}.txt;
  rm Codeml_*/${i}.txt;
done
