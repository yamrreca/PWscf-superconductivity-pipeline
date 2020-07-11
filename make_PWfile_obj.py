import re

"""Creates a PWfile object from a file whose path is given as a command line argument.

NOTE: It appears to not interact all that well with python 2.6.6, since it crashes when attempting to parse some files, but on my computer with python 3.7.4 the same files are parsed without issue. However, this only happens with some files."""

class PWfile:
    """Object with several of the parameters of a file intended for use with the PWscf module of Quantum ESPRESSO.

    Attributes:
    -filename: string. The name of the file.
    -contents: string. The contents of the file.
    -CONTROL: string. The &CONTROL namelist.
    -SYSTEM: string. The &SYSTEM namelist.
    -ELECTRONS: string. The &ELECTRONS namelist.
    -ATOMIC_SPECIES: string. The ATOMIC_SPECIES card.
    -ATOMIC_POSITIONS: string. The ATOMIC_POSITIONS card.
    -K_POINTS: string. The K_POINTS card.
    -ibrav: int. The file's ibrav parameter.
    -nat: int. The file's nat parameter, containing the number of atoms in a single cell.
    -ntyp: int. The file's ntyp parameter, containing the number of chemical elements in the system.
    -elem: set containing the different chemical elements.
    -calculation: string. The file's calculation parameter.
    -CELL: string. The &CELL optional namelist.
    -IONS: string. The &IONS optional namelist.
    -CELL_PARAMETERS: string. The CELL_PARAMETERS optional card.
    -OCCUPATIONS: string. The OCCUPATIONS optional card.
    -CONSTRAINTS: string. The CONSTRAINTS optional card.
    --ATOMIC_FORCES: string. The ATOMIC_FORCES optional card.
    
    NOTE: All blocks are the effective blocks, not the literal ones as they appear in the files. An effective block has all the same lines except for those that are commented out or are empty.
    """

    def __init__(self, filename, CONTROL, 
            SYSTEM, ELECTRONS, ATOMIC_SPECIES,
            ATOMIC_POSITIONS, K_POINTS,
            ibrav, nat, ntyp, elem,
            calculation, CELL = None, IONS = None,
            CELL_PARAMETERS = None, OCCUPATIONS = None,
            CONSTRAINTS = None, ATOMIC_FORCES = None):
        self.filename = filename
        self.CONTROL = CONTROL
        self.SYSTEM = SYSTEM
        self.ELECTRONS = ELECTRONS
        self.ATOMIC_SPECIES = ATOMIC_SPECIES
        self.ATOMIC_POSITIONS = ATOMIC_POSITIONS
        self.K_POINTS = K_POINTS
        self.ibrav = ibrav
        self.nat = nat
        self.ntyp = ntyp
        self.elem = elem
        self.calculation = calculation
        self.CELL = CELL
        self.IONS = IONS
        self.CELL_PARAMETERS = CELL_PARAMETERS
        self.OCCUPATIONS = OCCUPATIONS
        self.CONSTRAINTS = CONSTRAINTS
        self.ATOMIC_FORCES = ATOMIC_FORCES

    def get_namelists(self):
        """Returns the effective namelists, i.e. without empty or commented lines."""

        #Making a list of the namelists
        namelists = [self.CONTROL, self.SYSTEM, 
                self.ELECTRONS, self.IONS, self.CELL]

        #Eliminating empty blocks
        namelists = [namelist for namelist in namelists if namelist is not None]

        return namelists


    def get_cards(self):
        """Returns the effective cards, i.e. without empty or commented lines."""

        #Making a list of the cards
        cards = [self.ATOMIC_SPECIES, self.ATOMIC_POSITIONS, self.K_POINTS, self.CELL_PARAMETERS, self.OCCUPATIONS, self.CONSTRAINTS, self.ATOMIC_FORCES]

        #Eliminating empty blocks
        cards = [card for card in cards if card is not None]

        return cards


    def reconstruct(self):
        """Returns the reconstructed file contents as a string."""

        ##Reconstructing the namelist section
        namelists = self.get_namelists()

        ##Adding a newline at the beginning of each row save the first.
        #for i, namelist in enumerate(namelists[1:]):
        #    namelists[i+1] = '\n' + namelist

        ##Namelists end with a /
        namelist_section = '\n/\n'.join(namelists)
        ##The last namelist didn't get a /
        namelist_section += '\n/\n'

        ##Reconstructing the cards section
        cards = self.get_cards()
        cards_section = '\n'.join(cards)

        ##Putting them together. Note that all newlines are already in place.
        contents = namelist_section + cards_section

        return contents


        
def get_effective_block(block):
    """Returns its input without empty or commented lines.

    Arguments:
    -block: string. A namelist or card.

    Returns:
    -eff_block: The input without empty or commented lines."""

    #Extracting the lines
    if block is None: #e.g. optional cards.
        return None

    block_rows = block.split('\n')
    #Eliminating empty lines
    eff_block_list = [row for row in block_rows if row.strip() != '']
    #Eliminating commented lines
    eff_block_list = [row for row in eff_block_list if not row.strip().startswith('#')]

    eff_block = '\n'.join(eff_block_list)

    return eff_block


def get_namelists(file):
    """Obtains the file's namelists, i.e. the blocks at the beginning that start with an &.

    Parameters:
    -file: string. The contents of the PWscf file
    Returns:
    -namelists: list containing the blocks.
    """

    ##Namelists start with an & and end with a /
    namelist_regex = re.compile(r'&.+?\n\s*/', re.DOTALL)
    namelists = namelist_regex.findall(file) #Note that the entries don't include the last linebreaks.

    ##Removing the / and newline at the end of each namelist
    namelists = [namelist[:-2] if namelist.endswith('/') else namelist[:-1] for namelist in namelists]

    ##Quick check to ensure the file is PWscf input.
    try:
        CONTROL_flag = namelists[0].startswith('&CONTROL')
        SYSTEM_flag = namelists[1].startswith('&SYSTEM')
        ELECTRONS_flag = namelists[2].startswith('&ELECTRONS')
    except IndexError: #If there aren't enough matches
        print('Not a PWscf input file.')
        exit()

    if not (CONTROL_flag and SYSTEM_flag and ELECTRONS_flag):
        print('Not a PWscf input file')
        exit()

    return namelists


def get_opt_namelists(namelists):
    """Gets the &CELL and &IONS namelist, if they exist. Returns None if they don't.

    Arguments:
    -namelists: list containing each namelist.

    Returns:
    -CELL: string. The &CELL block
    -IONS: string: THE $IONS block
    """

    CELL = None #Value by default
    IONS = None

    if len(namelists) > 3: #Check to see if there even are optional namelists
        for i, row in enumerate(namelists):
            if row.startswith('&CELL'):
                CELL = namelists[i]
            elif row.startswith('&IONS'):
                IONS = namelists[i]

    return CELL, IONS


def get_SPECIES_card(file):
    """Returns the file's ATOMIC_SPECIES card.

    Parameters:
    -file: string. The contents of the PWscf file
    Returns:
    -ATOMIC_SPECIES: string containing the block.
    """

    SPECIES_regex = re.compile(r'ATOMIC_SPECIES.+ATOMIC_POSITIONS', re.DOTALL)
    match = SPECIES_regex.search(file)

    if match is None:
        return None
    else:
        ATOMIC_SPECIES = match.group(0)

    ##Removing the ATOMIC_POSITIONS at the end
    ATOMIC_SPECIES = ATOMIC_SPECIES[:-17] #-17 to eliminate the newline at the end of the last line.

    return ATOMIC_SPECIES

def extract_SPECIES_elem(ATOMIC_SPECIES):
    """Obtains the set of elements in the ATOMIC_SPECIES
    block

    Arguments:
    -SPECIES_block: string containing the ATOMIC_SPECIES block

    Returns:
    SPECIES_elem: set with the elements
    """

    SPECIES_elem_regex = re.compile(r'\n\s*(\w{1,2})\s+\d+\.\d+\s+.*\.UPF')
    SPECIES_elem = SPECIES_elem_regex.findall(ATOMIC_SPECIES)

    return SPECIES_elem


def get_POS_card(file):
    """Returns the file's ATOMIC_POSITIONS card.

    Parameters:
    -file: string. The contents of the PWscf file
    Returns:
    -ATOMIC_POSITIONS: string containing the block.
    """

    POS_regex = re.compile(r'ATOMIC_POSITIONS.+K_POINTS', re.DOTALL)
    match = POS_regex.search(file)

    if match is None:
        return None
    else:
        ATOMIC_POSITIONS = match.group(0)


    ##Removing the K_BLOCKS part at the end.
    ATOMIC_POSITIONS = ATOMIC_POSITIONS[:-9] #9 so as not to include the newline at the end of the last line.

    return ATOMIC_POSITIONS


def extract_POS_lines(ATOMIC_POSITIONS):
    """Obtains the lines in the ATOMIC_POSITIONS
    block as a list.

    Arguments:
    -ATOMIC_POSITIONS: string containing the ATOMIC_POSITIONS block.

    Returns:
    POS_lines: list containing the lines.
    """

    POS_lines_regex = re.compile(r"""(\w{1,2})\s+
                    (-?\d+\.\d+E?-?\d*)\s* #Accounts for scientific notation.
                    (-?\d+\.\d+E?-?\d*)\s*
                    (-?\d+\.\d+E?-?\d*)""",
                    re.VERBOSE)

    POS_lines = POS_lines_regex.findall(ATOMIC_POSITIONS) #Extracts only the element (in parentheses)

    return POS_lines


def get_K_POINTS_card(file):
    """Returns the file's K_POINTS card.

    Parameters:
    -file: string. The contents of the PWscf file
    Returns:
    -K_POINTS: string containing the block.
    """

    ##Due to the fact that all cards following K_POINTS are optional, as well as the fact that the K_POINTS card contains only the introductory line and one more defining the points, the search merely matches the K_POINTS line and the one after it.
    K_POINTS_lines = [] #In case the following loop doesn't find the line.
    file_lines = file.split('\n') #List with the file's lines
    for i, line in enumerate(file_lines):
        if line.strip().startswith('K_POINTS'): #strip removes whitespace
            K_POINTS_lines = file_lines[i:i+2]
            break

    if K_POINTS_lines == []:
        return None
    else:
        K_POINTS = '\n'.join(K_POINTS_lines)
        if K_POINTS.endswith('\n'):
            return K_POINTS[:-1] #Remove the newline
        else:
            return K_POINTS


def get_CELL_PARAMS_card(file):
    """Returns the file's CELL_PARAMETERS card.

    Parameters:
    -file: string. The contents of the PWscf file
    Returns:
    -CELL_PARAMETERS: string containing the block.
    """

    ##Due to the fact that all cards following CELL_PARAMETERS are optional (including CELL_PARAMETERS itself), as well as the fact that the CELL_PARAMETERS card contains only the introductory line and three more defining the params, the search merely matches the CELL_PARAMETERS line and the three after it.
    CELL_PARAMETERS_lines = [] #In case the following loop doesn't find the line.
    file_lines = file.split('\n') #List with the file's lines
    for i, line in enumerate(file_lines):
        if line.strip().startswith('CELL_PARAMETERS'): #strip removes whitespace
            CELL_PARAMETERS_lines = file_lines[i:i+4]
            break

    if CELL_PARAMETERS_lines == []:
        return None
    else:
        CELL_PARAMETERS = '\n'.join(CELL_PARAMETERS_lines)
        if CELL_PARAMETERS.endswith('\n'):
            return CELL_PARAMETERS[:-1] #Remove the newline
        else:
            return CELL_PARAMETERS


def get_nat(SYSTEM):
    """Gets the number of atoms from the nat line. For extra security, this function parses only the &SYSTEM namelist, where the parameter should be, instead of the entire file.
    
    Arguments:
    -SYSTEM: string. The &SYSTEM namelist.
    
    Returns:
    -nat: The number of atoms in the file.
    """
    
    nat_regex = re.compile(r'\n\s*nat\s*=\s*(\d+)\s*,')
    match = nat_regex.search(SYSTEM)
    if match is not None:
        nat = match.group(1)
        return int(nat)
    else:
        return None


def get_ntyp(SYSTEM):
    """Gets the number of elements from the ntyp line. For extra security, this function parses only the &SYSTEM namelist, where the parameter should be, instead of the entire file.
    
    Arguments:
    -SYSTEM: string. The &SYSTEM namelist.
    
    Returns:
    -ntyp: The number of elements in the file.
    """
    
    ntyp_regex = re.compile(r'\n\s*ntyp\s*=\s*(\d+)\s*,')
    match = ntyp_regex.search(SYSTEM)
    if match is not None:
        ntyp = match.group(1)
        return int(ntyp)
    else:
        return None


def get_calculation(CONTROL):
    """Extracts the type of calculation from the calculation parameter in a PWscf file. For extra security, this function parses only the &CONTROL namelist, where the parameter should be, instead of the entire file.
    
    Arguments:
    -CONTROL: String with the contents of the &CONTROL namelist.
    
    Returns:
    calculation: string. The value of the file's calculation parameter.
    """
    
    calculation_regex = re.compile(r"\n\s*calculation\s*=\s*'(.+)'\s*,")
    match = calculation_regex.search(CONTROL)
    if match is not None:
        calculation = match.group(1)
    else:
        calculation = None

    return calculation

def get_ibrav(SYSTEM):
    """Extracts the ibrav parameter from a PWscf file. For extra security, this function parses only the &SYSTEM namelist, where the parameter should be, instead of the entire file.
    
    Arguments:
    -SYSTEM: string. The &SYSTEM namelist.

    Returns:
    -ibrav: int. The value of the ibrav parameter.
    """

    ibrav_regex = re.compile(r'\n\s*ibrav\s*=\s*(\d)\s*,')
    match = ibrav_regex.search(SYSTEM)
    if match is not None:
        ibrav = int(match.group(1))
        return ibrav
    else:
        return None


def make_PWfile_obj(filename):
    """Creates a PWfile object.

    Arguments:
    -filename: string. The name of the file from which to create a PWfile object.

    Returns:
    -pwfile: A PWfile object with the attributes of the input.
    """

    ##First: get the contents.
    with open(filename, 'r') as archivo:
        contents = archivo.read()

    ##Second: get the required and optional namelists.
    namelists = get_namelists(contents)
    CONTROL, SYSTEM, ELECTRONS = namelists[:3] #Required
    CELL, IONS = get_opt_namelists(namelists) #Optional    

    ##Third: get the obligatory cards.
    ATOMIC_SPECIES = get_SPECIES_card(contents)
    ATOMIC_POSITIONS = get_POS_card(contents)
    K_POINTS = get_K_POINTS_card(contents)

    ##Fourth: get the optional cards.
    CELL_PARAMETERS = get_CELL_PARAMS_card(contents)

    ##Fifth: get the effective blocks.
    CONTROL = get_effective_block(CONTROL)
    SYSTEM = get_effective_block(SYSTEM)
    ELECTRONS = get_effective_block(ELECTRONS)
    CELL = get_effective_block(CELL)
    IONS = get_effective_block(IONS)
    ATOMIC_SPECIES = get_effective_block(ATOMIC_SPECIES)
    ATOMIC_POSITIONS = get_effective_block(ATOMIC_POSITIONS)
    K_POINTS = get_effective_block(K_POINTS)
    CELL_PARAMETERS = get_effective_block(CELL_PARAMETERS)

    ##Sixth: get the calculation parameter.
    calculation = get_calculation(CONTROL)

    ##Seventh: get the ibrav parameter.
    ibrav = get_ibrav(contents)

    ##Fourth: get the elem parameter, with the file's chemical elements. This parameter can only be unambiguously calculated if the chemical elements of both the ATOMIC_POSITIONS and ATOMIC_SPECIES blocks are the same, which in turn requires them to be written correctly. 

    if ATOMIC_POSITIONS is not None: #Only in this case does it make sense to obtain the elements.
        POS_lines = extract_POS_lines(ATOMIC_POSITIONS) #List with the individual lines of the block
        if POS_lines != []: #If the block wasn't malformed
            POS_lines = [list(line) for line in POS_lines]#Turn from a list of tuples to a list of lists
            POS_elem = [line[0] for line in POS_lines] #The chemical elements
            POS_elem_flag = True #To decide whether or not to compare POS_elem and SPECIES_elem to set the value of elem. Avoids an exception.
        else:
            POS_elem_flag = False
    else:
        POS_elem_flag = False

    if ATOMIC_SPECIES is not None:
        SPECIES_elem = extract_SPECIES_elem(ATOMIC_SPECIES)
        if SPECIES_elem != []:
            SPECIES_elem_flag = True 
        else:
            SPECIES_elem_flag = False
    else:
        SPECIES_elem_flag = False

    if POS_elem_flag and SPECIES_elem_flag:
        if set(POS_elem) == set(SPECIES_elem): #If the two blocks define the same elements
            elem = set(SPECIES_elem)
        else:
            elem = {}
    else:
        elem = {}

    ##Fifth: Get the nat parameter.
    nat = get_nat(SYSTEM)

    ##Sixth: Get the ntyp parameter
    ntyp = get_ntyp(SYSTEM)

    ##Final: Make the object.
    pwfile = PWfile(filename, CONTROL, 
            SYSTEM, ELECTRONS, ATOMIC_SPECIES,
            ATOMIC_POSITIONS, K_POINTS,
            ibrav, nat, ntyp, elem,
            calculation, CELL = CELL, IONS = IONS,
            CELL_PARAMETERS = CELL_PARAMETERS)

    return pwfile


if __name__ == '__main__':
    from sys import argv
    if len(argv) < 2:
        print('Usage: check_PWfile [Path to file]')
        exit()
    else:
        file = make_PWfile_obj(argv[1])