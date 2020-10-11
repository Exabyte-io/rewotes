import re

def formula_has_repeating_unit(formula):
    if "(" in formula:
        return True
    else:
        return False

def formula_has_number(formula):
    for ii in formula:
        if ii.isnumeric():
            return True
    return False

def simple_formula_split(formula):
    return re.sub( r"([A-Z])", r" \1", formula).split()

def numbered_formula_split(formula):
    number_split_strings = re.sub( r"([0-9])", r" \1", formula).split()
    split_formula = []
    for string in number_split_strings:
        tmp_splits = re.sub( r"([A-Z])", r" \1", string).split()
        for tmp_split in tmp_splits:
            split_formula.append(tmp_split)
    return split_formula

def repeating_unit_formula_split(formula):
    pass

def get_stoichiometry_for_simple_formula(split_formula,stoichiometry):
    for string in split_formula:
        stoichiometry[string] = 1

def get_stoichiometry_for_numbered_formula(split_formula,stoichiometry):
    for idx, string in enumerate(split_formula):
        if string[0].isupper() and not string.isnumeric():
            if idx+1 != len(split_formula) and not split_formula[idx+1][0].isupper():
                stoichiometry[string] = int(split_formula[idx+1])
            else:
                stoichiometry[string] = 1
        elif not string.isnumeric():
            stoichiometry[string] = 1

def get_stoichiomertry(formula):
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
    stoichiometry = get_stoichiomertry(formula)
    num_atoms = [ atoms for (element, atoms) in stoichiometry.items() ]
    num_atoms_total = sum(num_atoms)
    norm_stoichiometry = { element:atoms/num_atoms_total for (element, atoms) in stoichiometry.items()} 
    return norm_stoichiometry
