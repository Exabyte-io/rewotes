def encode_composition(composition_dict):
    if type(composition_dict) != dict:
        raise TypeError("Expected dict. Found: " + str(type(composition_dict)))
    for key in composition_dict:
        if type(key) != int:
            raise TypeError("Expected key to be int. Found: " + str(type(key)))
        if key < 1 or key > 118:
            raise ValueError("Expected key to be between 1 and 118. Found: " + str(key))
    byte_values = []
    for key in sorted(list(composition_dict)):
        value = composition_dict[key]
        if type(value) != int:
            raise TypeError("Expected value to be int. Found: " + str(type(value)))
        if value < 1:
            continue
        if value > 255:
            raise ValueError("Expected value to be no greater than 255. Found: " + str(value))
        byte_values.append(key)
        byte_values.append(value)
    return bytes(byte_values)

def decode_composition(composition_bytes):
    if len(composition_bytes) % 2:
        raise ValueError("Expected length to be even. Found: " + str(len(composition_bytes)))
    composition_dict = {}
    for i in range(0,len(composition_bytes),2):
        key = composition_bytes[i]
        value = composition_bytes[i + 1]
        if key < 1 or key > 118:
            raise ValueError("Expected key to be between 1 and 118. Found: " + str(key))
        if value < 1:
            continue
        composition_dict[key] = value
    return composition_dict
