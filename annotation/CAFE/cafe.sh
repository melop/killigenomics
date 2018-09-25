#species to include
SP='AAU CTO NFZgenbank NOR PLP astyanax gadus medaka poecilia spotted_gar stickleback takifugu tetraodon tilapia xmac zebrafish'

#make output for error estimation
mkdir -p reports/errest_sp/

python /software/source/CAFE/cafe_tutorial/python_scripts/cafetutorial_mcl2rawcafe.py -i out.seq.mci.I20 -o unfiltered_cafe_input.txt -sp "$SP"
sed -i 's/spotted_gar/spottedgar/g' unfiltered_cafe_input.txt # change name to exclude underscore!

python cafetutorial_clade_and_size_filter.py -i unfiltered_cafe_input.txt -o filtered_cafe_input.txt -s

python caferror.py -i cafe_errest.sh -s 1 -v 0 -f 1 -d reports/errest_sp

python /software/source/CAFE/cafe_tutorial/python_scripts/cafetutorial_report_analysis.py -i reports/errest_sp/cafe_final_report.cafe -o reports/errest_sp/cafe_final_summary

