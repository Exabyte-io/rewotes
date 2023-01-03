from material_pb2 import Material as protobuf_material

class Material(object):
    def __init__(self, initializer = None):
        if initializer is None:
            self._value = protobuf_material()
        elif type(initializer) == bytes:
            self._value = protobuf_material()
            self._value.ParseFromString(initializer)
        else:
            raise TypeError("Expected bytes. Found: " + str(type(initializer)))
    def serialize(self):
        return self._value.SerializeToString()
    def _get_material_id(self):
        return self._value.material_id
    def _set_material_id(self, value):
        if value is None:
            self._value.material_id = ""
        elif type(value) == str:
            self._value.material_id = value
        else:
            raise TypeError("Expected str. Found: " + str(type(value)))
    def _del_material_id(self):
        self._set_material_id(None)
    material_id = property(_get_material_id, _set_material_id, _del_material_id)
    def _get_optional_double(self, property_name):
        optional_double = getattr(self._value, property_name)
        if optional_double.is_defined:
            return optional_double.value
        else:
            return None
    def _set_optional_double(self, property_name, value):
        optional_double = getattr(self._value, property_name)
        if value is None:
            optional_double.value = 0
            optional_double.is_defined = False
        else:
            try:
                optional_double.value = float(value)
                optional_double.is_defined = True
            except:
                raise TypeError("Expected type which could be cast to float. Found: " + str(type(value)))
    def _del_optional_double(self, property_name):
        self._set_optional_double(property_name, None)

def apply_optional_double_properties():
    for property_name in ["band_gap"]:
        setattr(Material, property_name, property(
            lambda self: self._get_optional_double(property_name),
            lambda self, value: self._set_optional_double(property_name, value),
            lambda self: self._del_optional_double(property_name)
        ))

apply_optional_double_properties()
