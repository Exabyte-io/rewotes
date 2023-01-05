from .material_pb2 import Material as protobuf_material
from .composition_encoder import encode_composition, decode_composition

class Material(object):
    def __init__(self, initializer = None):
        self.__dict__["_values"] = protobuf_material()
        if initializer is None:
            pass
        elif type(initializer) == bytes:
            self.__dict__["_values"] .ParseFromString(initializer)
        else:
            raise TypeError("Expected bytes. Found: " + str(type(initializer)))
    def __getattr__(self, property_name):
        if property_name == "composition":
            return decode_composition(self.__dict__["_values"].composition)
        elif property_name == "composition_reduced":
            return decode_composition(self.__dict__["_values"].composition_reduced)
        else:
            return getattr(self.__dict__["_values"] , property_name)
    def __setattr__(self, property_name, property_value):
        if property_name == "composition":
            return setattr(self.__dict__["_values"], property_name, encode_composition(property_value))
        elif property_name == "composition_reduced":
            return setattr(self.__dict__["_values"], property_name, encode_composition(property_value))
        else:
            return setattr(self.__dict__["_values"] , property_name, property_value)
    def __delattr__(self, property_name):
        return delattr(self.__dict__["_values"] , property_name)
    def serialize(self):
        return self.__dict__["_values"] .SerializeToString()
    def save_to_file(self, filename):
        if type(filename) != str:
            raise TypeError("Expected str. Found: " + str(type(filename)))
        else:
            with open(filename, "wb") as file:
                file.write(self.serialize())
    @staticmethod
    def load_from_file(filename):
        if type(filename) != str:
            raise TypeError("Expected str. Found: " + str(type(filename)))
        else:
            with open(filename, "rb") as file:
                return Material(file.read())
    # def decode_composition(self):
    #     return decode_composition(self.composition)
    # def encode_composition(self, value):
    #     self.composition = encode_composition(value)
    # def decode_composition_reduced(self):
    #     return decode_composition(self.composition_reduced)
    # def encode_composition_reduced(self, value):
    #     self.composition_reduced = encode_composition(value)
    # composition = property(
    #     lambda self: decode_composition(getattr(self.__dict__["_values"], "composition")),
    #     lambda self, value: setattr(self.__dict__["_values"], "composition", encode_composition(value)),
    #     lambda self: setattr(self.__dict__["_values"], "composition", b'')
    # )
    # composition_reduced = property(
    #     lambda self: decode_composition(getattr(self.__dict__["_values"], "composition_reduced")),
    #     lambda self, value: setattr(self.__dict__["_values"], "composition_reduced", encode_composition(value)),
    #     lambda self: setattr(self.__dict__["_values"], "composition_reduced", b'')
    # )

# def add_composition_properties():
#     for property_name in ["composition", "composition_reduced"]:
#         setattr(Material, property_name, property(
#             lambda self: decode_composition(getattr(self.__dict__["_values"], property_name)),
#             lambda self, value: setattr(self.__dict__["_values"], property_name, encode_composition(value)),
#             lambda self: setattr(self.__dict__["_values"], property_name, b''),
#         ))

# add_composition_properties()
