#!/bin/bash

declare -a files #Will contain the files to be run, in order.
declare -a temp #Will contain the reply to a prompt, to ensure exactly three options were input
i=0 #index of the files array

#Obtaining the necessary files
while true; do
	#Ask the user for the script with the instructions to run, the input file, and the corresponding output file
	echo -n "script-to-run input-file output-file  > "
	read -a temp #Assign the values to the array temp
	while [[ ${#temp[@]} -ne 3 ]]; do #If more or less than 3 options were input
		echo "Please input exactly 3 options"
		echo -n "script-to-run input-file output-file  > "
		read -a temp
	done
	files+=(${temp[@]}) #Append the values at the end of the files array

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
		unset temp
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
	fi

	##Update the coordinates and run the script if it is a PWscf file, otherwise
	###just run. To determine if it is such a file, look for the ATOMIC_POSITIONS block
	###for no reason other than the namelists have more generic names
	###and I don't know if they may occur in other types of files.

	echo "Iteration: $(($i/3))"
	echo "Input file: $in_file"
	echo "Output file: $out_file"
	echo "Script: $script"
	echo "Prev input file: $in_file_prev"
	echo "Prev output file: $out_file_prev"
	ATOMIC_POSITIONS=$(grep "ATOMIC_POSITIONS" "$in_file_prev")
	#Note that the above will also return an empty string if in_file_prev is empty
	if [[ -z $ATOMIC_POSITIONS ]]; then #If it isn't PWscf
		#Run the script
		./$script
	else
		#update_parameters
		echo "Updating parameters"
		#Run the script
		./$script
	fi

	read
done
