declare -A files #Will contain the files to be run, in order.
row=0 #rows of the files array

#Obtaining the necessary files
while true; do
	#Ask the user for the input file, the corresponding output file, and the script with the instruction to run.
	echo -n "input-file output-file script-to-run > "
	read file["$row,0"] file["$row,1] file["$row,2"]

	##Ask the user if there are any more calculations to be made.
	read -p "Are there more files to run (y/n)"
	while [[ "$REPLY" != "y"]] || [[ "$REPLY" != "n" ]]; do
		echo "Invalid answer."
		read -p "Are there more files to run (y/n)"
	done

	if [[ "$REPLY" == "y" ]]; then
		$row+=1
	elif [[ "$REPLY" == "n" ]]; then
		NUM_ROWS=$row #For clarity, this is the max number of rows in the array
		unset row
		break
	fi
done

for i in {0..$NUM_ROWS}; do
	in_file=file["$i,0"] #Input file for current iteration
	if [[ "$i" != 0 ]]; then #If it isn't the first iteration
		in_file_prev=file["$(($i - 1)),0"] #Input file for previous iteration
		out_file_prev=file["$(($i - 1)),1"] #Corresponding output file
	fi
	script=file["$i,2"] #Script to be run in the current iteration

	##Update the coordinates and run the script if it is a PWscf file, otherwise
	###just run. To determine if it is such a file, look for the ATOMIC_POSITIONS block
	###for no reason other than the namelists have more generic names
	###and I don't know if they may occur in other types of files.

	ATOMIC_POSITIONS=$(grep "ATOMIC_POSITIONS" "$in_file_prev")
	#Note that the above will also return an empty string if in_file_prev is empty
	if [[ -z $ATOMIC_POSITIONS ]]; then #If it isn't PWscf
		#Run the script
		$script
	else
		#update_parameters
		update_coords $in_file $out_file_prev $in_file_prev
		#Run the script
		$script
	fi
done