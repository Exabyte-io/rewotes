import re

def formula_has_repeating_unit(formula):
    """
    Decides if a formula has a repeating unit in it.

    Arguments:
        formula (str)

    Returns:
        boolean
    """
    if "(" in formula:
        return True
    else:
        return False

def formula_has_number(formula):
    """
    Decides if a formula has a numeric value in it.

    Arguments:
        formula (str)

    Returns:
        boolean
    """
    for ii in formula:
        if ii.isnumeric():
            return True
    return False

def simple_formula_split(formula):
    """
    Splits a formula into substrings if it doesn't have a numeric value it.

    Arguments:
        formula (str)

    Returns:
        (list of str)
    """
    return re.sub( r"([A-Z])", r" \1", formula).split()

def numbered_formula_split(formula):
    """
    Splits a formula into substrings if it has numeric value it.

    Arguments:
        formula (str)

    Returns:
        (list of str)
    """
    number_split_strings = re.sub( r"([0-9])", r" \1", formula).split()
    split_formula = []
    for string in number_split_strings:
        tmp_splits = re.sub( r"([A-Z])", r" \1", string).split()
        for tmp_split in tmp_splits:
            split_formula.append(tmp_split)
    return split_formula

def get_stoichiometry_for_simple_formula(split_formula, stoichiometry):
    """
    Populates input stoichiometry dictionary for a given formula with no numeric values
    which has been split into its constituent elements.

    Arguments:
        split_formula (list of str): output from simple_formula_split
            Example: ["Cu", "S"] 
        stoichiometry (dict): empty dictionary 
            Example after function: { "Cu": 1, "S": 1 }
    """
    for string in split_formula:
        stoichiometry[string] = 1

def get_stoichiometry_for_numbered_formula(split_formula,stoichiometry):
    """
    Populates input stoichiometry dictionary for a given formula with numeric values
    which has been split into its constituent elements and subscripts.

    Arguments:
        split_formula (list of str): output from simple_formula_split
            Example: ["Cu", 2, "S", 3] 
        stoichiometry (dict): empty dictionary 
            Example after function: { "Cu": 2, "S": 3 }
    """
    for idx, string in enumerate(split_formula):
        if string[0].isupper() and not string.isnumeric():
            if idx+1 != len(split_formula) and not split_formula[idx+1][0].isupper():
                stoichiometry[string] = int(split_formula[idx+1])
            else:
                stoichiometry[string] = 1
        elif not string.isnumeric():
            stoichiometry[string] = 1

def get_stoichiomertry(formula):
    """
    Creates a stoichiometry dictionary which contains a key value pair of the 
    constituent elements and their number of atoms.

    Arguments:
        formula (str): a formula in string format with no spaces or formatting 
            Example: "Cu2S3"

    Returns:
        stoichiometry (dict): map from element to number of atoms 
            Example: { "Cu": 2, "S": 3 }
    """
    stoichiometry = dict()

    if formula_has_repeating_unit(formula):
        split_formula = formula.split("(")

        # before repeating unit
        tmp_stoichiometry = dict()
        if formula_has_number(split_formula[0]): 
            tmp_split_formula = numbered_formula_split(split_formula[0])
            get_stoichiometry_for_numbered_formula(tmp_split_formula,tmp_stoichiometry)

        else: 
            tmp_split_formula = simple_formula_split(split_formula[0])
            get_stoichiometry_for_numbered_formula(tmp_split_formula,tmp_stoichiometry)

        for key, value in tmp_stoichiometry.items():
            stoichiometry[key] = value 

        # handle repeating unit
        repeat_unit, repeat_unit_amount = split_formula[1].split(")")

        tmp_stoichiometry = dict()
        if formula_has_number(repeat_unit): 
            tmp_split_formula = numbered_formula_split(repeat_unit)
            get_stoichiometry_for_numbered_formula(tmp_split_formula,tmp_stoichiometry)

        else: 
            tmp_split_formula = simple_formula_split(repeat)
            get_stoichiometry_for_numbered_formula(tmp_split_formula,tmp_stoichiometry)

        for key, value in tmp_stoichiometry.items():
            stoichiometry[key] = int(value) * int(repeat_unit_amount)
    
    elif formula_has_number(formula):
        split_formula = numbered_formula_split(formula)
        get_stoichiometry_for_numbered_formula(split_formula,stoichiometry)

    else:
        split_formula = simple_formula_split(formula)
        get_stoichiometry_for_simple_formula(split_formula,stoichiometry)

    return stoichiometry

def get_norm_stoichiomertry(formula):
    """
    Creates a normalized stoichiometry dictionary which contains a key value pair of the 
    constituent elements and their elemental ratio.

    Arguments:
        formula (str): a formula in string format with no spaces or formatting 
            Example: "Cu2S3"

    Returns:
        norm_stoichiometry (dict): map from element to number of atoms 
            Example: { "Cu": 0.4, "S": 0.6 }
    """
    stoichiometry = get_stoichiomertry(formula)
    num_atoms = [ atoms for (element, atoms) in stoichiometry.items() ]
    num_atoms_total = sum(num_atoms)
    norm_stoichiometry = { element:atoms/num_atoms_total for (element, atoms) in stoichiometry.items()} 
    return norm_stoichiometry
