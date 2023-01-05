from .material_pb2 import Material as protobuf_material

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
        return getattr(self.__dict__["_values"] , property_name)
    def __setattr__(self, property_name, property_value):
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
