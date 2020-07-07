import re
from sys import argv
import make_PWfile_obj

"""Performs the following checks to make sure its input 
PWscf file can be run without errors:

-Checks the filename and the calculation it implies (e.g. BN-scf.in implies regular scf) against the contents of the file (in particular the calculation parameter) to see if the latter makes sense.

-Checks to see if 1.- the ibrav parameter is specified and 2.- the file's structure fits its value (in particular the celldm and CELL_PARAMETERS parts).

-Checks that the number of atoms in ATOMIC_POSITIONS
is the same as that defined by the nat parameter.

-Performs the following checks on the ntyp parameter and the number of elements defined by the ATOMIC_SPECIES block: 

    1.-Checks that ntyp is defined well.
    2.- Checks that the ATOMIC_POSITIONS and ATOMIC_SPECIES blocks do not define different elements.
    3.- Checks that the value of ntyp equals the number of chemical elements.

-Detects and prints lines that are commented out. Note that comments at the end of each line are ignored.

IMPORTANT: For all blocks
except ATOMIC_SPECIES and below the parser expects
commas to be at the end of the line (before
comments)."""

"""To do:
-Add an object-oriented approach -- DONE
-Check that no lines start with a comment. -- DONE
-Check ibrav 0 and 4, which are the only ones I've used. --DONE
    -For ibrav 0 check for the existence and integrity of the
     CELL_PARAMETERS block --DONE
    -For ibrav 4 check for the existence and integrity of the
    celldm(1) and celldm(3) lines. --DONE
-Perform checks on the calculations. -- DONE
    -For scf, check that the wf_collect = .true. line is present. -- DONE
    -For scf-fit, check that the above line is absent. -- DONE
    -For vc-relax, check that the &Cell block is present. -- DONE
    -For relax, check that the above block is absent. -- DONE
-Implement make_PWfile function to return a PWfile object, for possible use in the future. -- DONE
-Fix the output to more closely resemble the order in which parameters appear in the actual files. -- DONE
"""

def get_filename_type(filename):
    """Extracts the type of calculation a PWscf file is meant
    to perform from its filename. Works with relax, 
    vc-relax, scf, and scf-fit.
    
    Arguments:
    -filename: string. The path to the PWscf file
    
    Returns:
    -filename_type: string equal to one of the calculations
    mentioned above. If none of them are in the filename then
    None is returned.
    """
    
    filename_regex = re.compile(r'vc[\-_]relax|relax|scf[\-_]fit|scf',
                               re.IGNORECASE)
    match = filename_regex.search(filename)
    if match is not None:
        filename_type = match.group(0).replace('_','-')
    else:
        print('Filename does not contain the type of calculation')
        filename_type = None

    return filename_type


def get_comments(file):
    """Extracts all comments from the file

    Arguments:
    -file: string with the contents of the file

    Returns:
    -comments: list with all lines with comments"""

    comments_regex = re.compile('\n\s*#.*') #Only considers comments
    #at the beginning of the line.
    comments = comments_regex.findall(file)

    #Remove whitespace
    comments = list(map(str.strip, comments))

    return comments


def check_ibrav(file):
    """Checks to see if a file's contents are in agreement with
    the value of the ibrav parameter. Supports ibrav 0 and 4.

    Arguments:
    -file: PWfile object with the info of the file.

    Returns:
    -None
    """

    ibrav = file.ibrav
    contents = file.contents
    
    if ibrav == 0:
        #Check for the existence of the CELL_PARAMETERS block
        CELL_regex = re.compile(r'\s*CELL_PARAMETERS\s*[\(\{]\s*\w+\s*[\)\}]\s*\#*.*\n')
        #Note that I attempted to check for the entire block but everything I tried either didn't work properly or mysteriously caused the program to slow to a crawl.  

        if CELL_regex.search(contents) is None:
            print('\nWarning: CELL_PARAMETERS block missing or malformed')
        
    elif ibrav == 4:
        #Check for the existence of the celldm(1) and (3) lines
        celldm1_regex = re.compile(r"""\s*[^\#]celldm\(1\)\s* #left-hand side
                        =\s* #Equals sign
                        \d+\.\d+""",    #RHS
                        re.VERBOSE)

        celldm3_regex = re.compile(r"""\s*[^\#]celldm\(3\)\s* #left-hand side
                        =\s* #Equals sign
                        \d+\.\d+""",    #RHS
                        re.VERBOSE)

        celldm1 = celldm1_regex.search(contents)
        celldm3 = celldm3_regex.search(contents)

        if (celldm1 is None) or (celldm3 is None):
            print('\nWarning: celldm parameters missing or malformed.')

    return None


def check_calculation(file):
    """Checks the integrity of the file depending on the
    calculation parameter:

    For scf, checks that the wf_collect = .true. line is present
    For scf-fit, checks that the above line is absent.
    For vc-relax, checks that the &Cell block is present.
    For relax, checks that the above block is absent.
    
    Parameters:
    -file: PWfile object.

    Returns:
    -None
    """

    calculation = file.calculation
    filename_type = get_filename_type(file.filename) #To check if it fits the calculation param.
    CONTROL = file.CONTROL
    CELL = file.CELL

    #Getting the wf_collect line
    scf_regex = re.compile(r"""\n\s*wf_collect\s* #left-hand side
                        =\s* #Equals sign
                        \.true\.\s*,""",    #RHS
                        re.VERBOSE)

    wf_collect = scf_regex.search(CONTROL)

    if calculation == 'scf': #Covers both scf and scf-fit
        if 'scf' not in filename_type: #Filename doesn't imply either scf
            print('\nWarning: filename does not fit calculation parameter.')
        if wf_collect is None: #Implies scf-fit
            if filename_type != 'scf-fit':
                print('wf_collect missing or malformed despite filename not implying scf-fit.')
        else: #scf
            if filename_type != 'scf':
                print('wf_collect line present despite filename not implying regular scf.')
    elif calculation == 'relax':
        if filename_type != 'relax':
            print('\nWarning: filename does not fit calculation parameter.')
        if CELL is not None:
            print("\nWarning: &CELL block not absent despite calculation = 'relax'")
    elif calculation == 'vc-relax':
        if filename_type != 'vc-relax':
            print('\nWarning: filename does not fit calculation parameter.')
        if CELL is None:
            print("\nWarning: &CELL block absent despite calculation = 'vc-relax'")
    else: #If it the calculation value isn't valid
        print('Warning: invalid calculation parameter.')

    return None


def compare_ntyp_and_elem(file):
    """Checks that the ntyp parameter takes on the same value as the number of elements.

    Arguments:
    -file: PWfile object
    """

    if len(file.elem) != file.ntyp:
          print('\nWarning: Different # of elements in ATOMIC_SPECIES than defined by ntyp.') 

    return


def check_nat(file):
    """Checks that the nat parameter in the input file takes on the same value as the number of lines in the ATOMIC_POSITIONS block."""

    ATOMIC_POSITIONS_lines = file.ATOMIC_POSITIONS.split('\n')[1:] #First line is the ATOMIC_POSITIONS intro.

    if file.nat is not None:
        print('nat: ' + str(file.nat))
        if len(ATOMIC_POSITIONS_lines) != file.nat:
            print('\nWarning: Different # of atoms in ATOMIC_POSITIONS than defined by nat.')
    else:
        print('\nWarning: the nat parameter is not specified.')
    return

def check_commas(file):
    """Checks that each line of each namelist ends with a comma.

    Arguments:
    -file: PWfile object
    """

    namelists = file.get_namelists() #List of namelists
    all_commas = True #To decide what to print

    for namelist in namelists:
        namelist_rows = namelist.split('\n')[1:] #The & line doesn't end in a commma.
        for row in namelist_rows:
            if ',' not in row:
                all_commas = False
                print('Warning: comma missing in line %s'%(row.strip()))
    if all_commas:
        print('All commas in place.')

    return

if __name__ == '__main__':
    if len(argv) < 2:
    	print('Usage: check_PWfile [Path to file]')
    	exit()
    else:
        file = make_PWfile_obj.make_PWfile_obj(argv[1])

    print("ATTENTION: THIS SIMPLE PROGRAM WAS MADE BY AN UNDERGRAD STUDENT WITH LITTLE UNDERSTANDING OF THE MEANING OF EACH LINE IN PWSCF INPUT FILES. IT IS NO SUBSTITUTE TO GOING OVER YOUR INPUTS YOURSELF. IT'S A SIMPLE SANITY CHECK, NOT A NANNY.")

    ##First: Check integrity of the calculation parameter (see function docstring).
    if file.calculation is not None:
        print('\ncalculation: ' + file.calculation)
        check_calculation(file)
    else:
        print('\nWarning: calculation parameter missing or malformed.')

    ##Second: check the ibrav parameter
    if file.ibrav is not None:
        print('\nibrav: ' + str(file.ibrav))
        check_ibrav(file)
    else:
        print('\nWarning: ibrav is not specified.')
    
    ##Third: Check that the number of atoms in ATOMIC_POSITIONS
    ##is the same as the one given by the nat parameter.
    check_nat(file)

    """Fourth: Perform the following checks on the ntyp and elem parameters: 

    1.- Check that elem is not empty
    2.-Check that ntyp is not empty
    3.- Check that the value of ntyp equals the length of elem (number of chemical elements)."""
                    
    #Check number 1
    if file.elem != {}:
        print('\nThe following elements are defined:')
        print('\n'.join(file.elem))
        #Check number 2
        if file.ntyp is not None:
           print('\nntyp: ' + str(file.ntyp))
           #Check number 3
           compare_ntyp_and_elem(file)
        else:
            print('\nWarning: ntyp line missing.')
    else:
        print('\nWarning: check the ATOMIC_POSITIONS and ATOMIC_SPECIES blocks. The chemical elements defined in them differ.')

    ##Seventh: Check that every line in every namelist ends with a comma.
    print()
    check_commas(file)
    print()

    ##Sixth: Print comments.
    comments = get_comments(file.contents)
    if comments == []:
        print('No commented lines')
    else:
        print('The following lines are commented:')
        print('\n'.join(comments))