files = [] #Will contain the files to be run, in order

while True:
	##Ask the user for the input file, the corresponding
	##output file, and the script with the instruction to run.
	in_file, out_file, script_file = input('input-file output-file script-to-run > ').split()
	files.append([in_file, out_file, script_file])

	##Ask the user if there are any more calculations
	##to be made.
	more_files = input('Are there more files to run? (y/n)').lower()
	while more_files != 'y' or more_files != 'n':
		print('Invalid answer.')
		more_files = input('Are there more files to run? (y/n)').lower()
	if more_files == 'y':
		pass
	elif more_files == 'n':
		break


for row in files:
	##Run the script

	##Run update_coords

