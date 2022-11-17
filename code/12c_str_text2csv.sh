#!/bin/bash
#4 OCT 2022
#CE Moore
#script to take STRUCTURE output files and make a csv to bring into R

#set up
cd ~/Documents/Projects/temporal_aztreefrog/structure_results
tK=2   #set the most likely K calculated previously
tR=3  #set the rep from highest likelihood across reps
RUN_FOLDER='r2021/r6789' #set the run to pull from (basically the folder name)

OUT_FILE=./${RUN_FOLDER}/Output_Files/Year_Pop_K${tK}_r${tR}_f

#find num. individuals in run
LINE=$(grep -hr "[0-9] individuals" "$OUT_FILE")
IND=$(grep -Eo "[0-9]{1,}" <<< "$LINE") #<<< is for string instead of file in these commands

#find where the inferred ancestry starts & ends
find_anc=$(grep -n "Inferred ancestry of individuals" "$OUT_FILE")
LINE_start=$(($(grep -Eo "[0-9]{1,}" <<< "$find_anc") + 1))
LINE_end=$(($LINE_start + $IND))

#grab only the lines of the inferred ancestry
sed -n "${LINE_start},${LINE_end}p;$((LINE_end + 1))q" "$OUT_FILE" > temp_file
#awk 'NR==${LINE_start}, NR==${LINE_end}; NR==${LINE_end} {exit}' "$OUT_FILE" > temp_file

#add separators
#sed -e "s/ /\t/g" temp_file > ./${RUN_FOLDER}/inferred_ancestry.csv


{
	read
	while IFS="" read -r p || [ -n "$p" ]
	do
		VALUES=$(sed "s/^.*: //" <<< "$p")
		VALUES_TAB=$(sed "s/ /,/" <<< "$VALUES")
		printf "${VALUES// /,}\n" "$p" >> ./${RUN_FOLDER}/inferred_ancestry.csv
	done
}< temp_file


rm temp_file

echo "inferred_ancestry.csv saved..."
