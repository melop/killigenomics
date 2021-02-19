IFS=$'\n' GLOBIGNORE='*' command eval  'arrSamples=($(cat list_PML.txt))'

for sSample in "${arrSamples[@]}"; do
	cat depthcounts/${sSample}* | awk '{sum+=$1} END{print sum/NR}'
done
