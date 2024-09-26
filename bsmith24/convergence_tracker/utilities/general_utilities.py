    
def read_file(filename):
    """
    """
    with open(filename) as file:
        filelines = file.readlines()
    return filelines


def is_converged(values, tolerance):
    prev_value = 0 
    for count, curr_value in enumerate(values):
        if abs(curr_value - prev_value) <= tolerance:
            return [True, values[count-1], count-1]
        else:
            prev_value = curr_value
    return [False, None, None]

