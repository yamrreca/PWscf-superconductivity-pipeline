#!/bin/bash

declare -a files #Will contain the files to be run, in order.
declare -a temp #Will contain the reply to a prompt, to determine if the user didn't input more or less than the required amount of info
i=0 #index of the files array

#Obtaining the necessary files
while true; do
	#Ask the user for the input file, the corresponding output file, and the script with the instruction to run.
	echo -n "script-to-run input-file output-file  > "
	read -a temp
	while [[ ${#temp[@]} -ne 3 ]]; do #If more than 3 options were input
		echo "Please input exactly 3 options"
		echo -n "script-to-run input-file output-file  > "
		read -a temp
	done
	files+=(${temp[@]})

	echo "succesfully obtained ${files[$i]} ${files[$(($i + 1))]} ${files[$(($i + 2))]}"
	if [[ $i -ne 0 ]]; then
		echo "Previous values: ${files[$(($i - 3))]} ${files[$(($i - 2))]} ${files[$(($i - 1))]}"
	fi

	##Ask the user if there are any more calculations to be made.
	read -p "Are there more files to run (y/n)? "
	echo
	while [[ "$REPLY" -ne "y" || "$REPLY" -ne "n" ]]; do
		echo "Invalid answer."
		read -p "Are there more files to run (y/n) "
	done

	if [[ "$REPLY" == "y" ]]; then
		i=$(($i+3))
	else 
		NUM_ENTRIES=$(($i+3)) #Number of entries in the array
		unset i 
		break
	fi
done

echo "NUM_ENTRIES: $NUM_ENTRIES"
#Run all programs in order
for (( i=0; i<$(($NUM_ENTRIES)); i+=3 )); do
	script=${files[$i]} #Script to be run in the current iteration
	in_file=${files[$(($i+1))]} #Input file for current iteration
	out_file=${files[$(($i+2))]} #Output file for current iteration
	if [[ "$i" -ne 0 ]]; then #If it isn't the first iteration
		in_file_prev=${files[$(($i - 2))]} #Input file for previous iteration
		out_file_prev=${files[$(($i - 1))]} #Corresponding output file
		##Get the prev in file's calculation parameter to determine if updating the
		##parameters of the current in file is necessary.
		calculation=$( python3 get_calculation "$in_file_prev")
	fi


	echo "Iteration: $(($i/3))"
	echo "Input file: $in_file"
	echo "Output file: $out_file"
	echo "Script: $script"
	echo "Prev input file: $in_file_prev"
	echo "Prev output file: $out_file_prev"
	echo "Prev calculation: "$calculation""

	if [[ "$calculation" == "relax" || "$calculation" == "vc-relax" ]]; then #If the prev file was optimization
		#update_parameters
		echo "Updating parameters"
		python3 ../update-coords/update_coords.py "$in_file" "$out_file_prev" "$calculation" #Update the parameters
		#Run the script
		./$script
	else
		#Run the script
		echo "Not updating parameters"
		./$script
	fi

	read
done
