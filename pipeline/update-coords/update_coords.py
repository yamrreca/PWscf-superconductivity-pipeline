import make_PWfile_obj
import re
from sys import argv

"""
Updates a PWscf file with the optimized coordinates from running a relax or vc-relax calculation. 
"""

def relax_update(in_file, relax_out, output_file = 'pipe_out.1'):
	"""Updates the .in file in the next step of the pipeline with the new coordinates as given by the relax.out file.

	Arguments:
	-in_file: PWfile object corresponding to the .in file to be updated.
	-relax_out: string. The filename of the relax.out file

	Returns: in_file: PWfile object with the new coordinates."""

	with open(relax_out, 'r') as archivo:
		relax = archivo.read()

	##Defining a regex to look for the ATOMIC_POSITIONS in 
	##the out file.
	out_regex = re.compile(r'\nBegin final coordinates.*(ATOMIC_POSITIONS.*)\nEnd final coordinates', re.DOTALL)

	match = out_regex.search(relax)
	if match is None:
		with open(output_file, 'a') as archivo:
			archivo.write('Fatal: no final coordinates in relax.out')
		exit()
	else:
		in_file.ATOMIC_POSITIONS = match.group(1)

	return in_file


def get_new_celldm(CELL_out, output_file = 'pipe_out.1'):
	"""
	Calculates the new celldm(1) and celldm(3) values from a vc-relax file. Only works if ibrav = 4.

	Arguments:
	-CELL_out: string. The CELL_parameters block in the vc-relax output file.
	-output_file: string. File to write output messages to.

	Returns:
	-celldm1: The value of the celldm(1) parameter.
	-celldm3: The value of the celldm(3) parameter.
	"""	

	##Defining a regex to extract the CELL_PARAMETERS block.
	
	##The only numbers that matter are the alat parameter,
	##the value in the upper-left of the array,
	##and the value in the lower-right of the array.

	##Extracting the rows
	CELL_out_rows = CELL_out.split('\n')

	##Extracting alat
	alat = float(re.search(r'alat\s*=\s*(\d+\.\d+)',CELL_out_rows[0]).group(1))

	##Extracting the first parameter of the array, corresponding
	##to the change in the horizontal vector on the
	##xy plane.
	da = float(CELL_out_rows[1].split()[0])

	##Extracting the last parameter in the array,
	##corresponding to the change in the z vector.
	dc = float(CELL_out_rows[3].split()[2])

	##Calculating the new celldm params.
	celldm1_new = da*alat
	celldm3_new = dc/da

	return celldm1_new, celldm3_new


def update_celldm(pw_file, celldm1_new, celldm3_new, output_file = 'pipe_out.1'):
	"""
	Updates the values of the celldm(1) and celldm(3) parameters of a PWfile object. Note that these parameters only exist if ibrav = 4.

	Arguments:
	-pw_file: PWfile object. The file to update.
	-celldm1: float. The new value of celldm(1).
	-celldm3: float. The new value of celldm(3).

	Returns:
	-pw_file: PWfile object with new values of celldm.
	"""

	##Getting the &SYSTEM namelist, where the celldm lines
	##are defined.
	SYSTEM = pw_file.SYSTEM
	celldm_regex = re.compile(r'''(.+) #before celldm
					\n\s*celldm\(1\)\s*=\s*\d+\.\d+,?
					\n\s*celldm\(3\)\s*=\s*\d+\.\d+,?
					\n(.+) #after celldm
					''', flags = re.DOTALL | re.VERBOSE)

	match = celldm_regex.search(SYSTEM)
	if match is None:
		with open(output_file, 'a') as archivo:
			archivo.write('Fatal: no celldm parameters in %s'%(pw_file.filename))
			exit()
	else:
		##Getting the different parts of the namelist
		bf = match.group(1)
		aft = match.group(2)
		#The celldm precision is set at 11 digits.
		celldm1_row_new = 'celldm(1) = %.11f'%(celldm1_new)
		celldm3_row_new = 'celldm(3) = %.11f'%(celldm3_new)

		##Creating the new &SYSTEM namelist
		new_SYSTEM = bf + '\n  ' + celldm1_row_new + '\n  ' + celldm3_row_new + '\n' + aft
		pw_file.SYSTEM = new_SYSTEM

	return pw_file


def vc_relax_update(in_file, vc_relax_out, output_file = 'pipe_out.1'):
	"""Updates the .in file in the next step of the pipeline with the new coordinates as given by the vc-relax.out file.

	Arguments:
	-in_file: PWfile object corresponding to the .in file to be updated.
	-vc_relax_out: string. The filename of the relax.out file

	Returns: in_file: PWfile object with the new coordinates."""

	with open(vc_relax_out, 'r') as archivo:
		vc_relax = archivo.read()

	##While the specific calculations depend on ibrav,
	##the regex is the same.
	out_regex = re.compile(r'\nBegin final coordinates.*\n(CELL_PARAMETERS.*)\n(ATOMIC_POSITIONS.*)\nEnd final coordinates', re.DOTALL)

	match = out_regex.search(vc_relax)
	if match is None:
		with open(output_file, 'a') as archivo:
			archivo.write('Fatal: no final coordinates in vc-relax.out')
		exit()
	else:
		CELL_PARAMETERS_new = match.group(1)
		ATOMIC_POSITIONS_new = match.group(2)

	if in_file.ibrav not in (0,4):
		print('ibrav %d not supported.'%(in_file.ibrav))
		exit()
	elif in_file.ibrav == 0:
		in_file.ATOMIC_POSITIONS = ATOMIC_POSITIONS_new
		in_file.CELL_PARAMETERS = CELL_PARAMETERS_new
	elif in_file.ibrav == 4:
		in_file.ATOMIC_POSITIONS = ATOMIC_POSITIONS_new
		
		##Calculating new celldm
		celldm1_new, celldm3_new = get_new_celldm(CELL_PARAMETERS_new, output_file = output_file)
		##Updating the celldm rows
		in_file = update_celldm(in_file, celldm1_new, celldm3_new, output_file = output_file)

	return in_file

if __name__ == '__main__':
	if len(argv) < 2:
		print('Usage: update_coords [input file to update] [output file from which to update] [calculation of previous input file]')
		exit()
	else:
		in_file = make_PWfile_obj.make_PWfile_obj(argv[1])
		out_filename = argv[2]
		prev_calculation = argv[3]

	print('in_file:', in_file)
	print('out_file:', out_filename)
	print('calculation:', prev_calculation)

	if prev_calculation == 'relax':
		in_file_new = relax_update(in_file, out_filename)
	elif prev_calculation == 'vc-relax':
		in_file_new = vc_relax_update(in_file, out_filename)
	else: ##If no optimization took place.
		exit()

	print('After if')
	##Saving the updated file.
	with open(in_file.filename + '-2', 'w') as archivo:
		archivo.write(in_file_new.reconstruct())