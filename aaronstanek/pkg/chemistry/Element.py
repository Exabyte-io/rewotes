class Element(object):
    def __init__(self, english_name, symbol, atomic_number):
        self._english_name = english_name
        self._symbol = symbol
        self._atomic_number = atomic_number
    english_name = property(
        lambda self: self._english_name,
        lambda self, value: None,
        lambda self: None
    )
    symbol = property(
        lambda self: self._symbol,
        lambda self, value: None,
        lambda self: None
    )
    atomic_number = property(
        lambda self: self._atomic_number,
        lambda self, value: None,
        lambda self: None
    )

class ElementLookupTable(object):
    def _put_element_in_table(self, english_name, symbol, atomic_number):
        element = Element(english_name, symbol, atomic_number)
        self._by_english_name[english_name] = element
        self._by_symbol[symbol] = element
        self._by_atomic_number[atomic_number] = element
    def __init__(self):
        self._by_english_name = {}
        self._by_symbol = {}
        self._by_atomic_number = {}
        self._put_element_in_table("Hydrogen", "H", 1)
        self._put_element_in_table("Helium", "He", 2)
        self._put_element_in_table("Lithium", "Li", 3)
        self._put_element_in_table("Beryllium", "Be", 4)
        self._put_element_in_table("Boron", "B", 5)
        self._put_element_in_table("Carbon", "C", 6)
        self._put_element_in_table("Nitrogen", "N", 7)
        self._put_element_in_table("Oxygen", "O", 8)
        self._put_element_in_table("Fluorine", "F", 9)
        self._put_element_in_table("Neon", "Ne", 10)
        self._put_element_in_table("Sodium", "Na", 11)
        self._put_element_in_table("Magnesium", "Mg", 12)
        self._put_element_in_table("Aluminium", "Al", 13)
        self._put_element_in_table("Silicon", "Si", 14)
        self._put_element_in_table("Phosphorus", "P", 15)
        self._put_element_in_table("Sulfur", "S", 16)
        self._put_element_in_table("Chlorine", "Cl", 17)
        self._put_element_in_table("Argon", "Ar", 18)
        self._put_element_in_table("Potassium", "K", 19)
        self._put_element_in_table("Calcium", "Ca", 20)
        self._put_element_in_table("Scandium", "Sc", 21)
        self._put_element_in_table("Titanium", "Ti", 22)
        self._put_element_in_table("Vanadium", "V", 23)
        self._put_element_in_table("Chromium", "Cr", 24)
        self._put_element_in_table("Manganese", "Mn", 25)
        self._put_element_in_table("Iron", "Fe", 26)
        self._put_element_in_table("Cobalt", "Co", 27)
        self._put_element_in_table("Nickel", "Ni", 28)
        self._put_element_in_table("Copper", "Cu", 29)
        self._put_element_in_table("Zinc", "Zn", 30)
        self._put_element_in_table("Gallium", "Ga", 31)
        self._put_element_in_table("Germanium", "Ge", 32)
        self._put_element_in_table("Arsenic", "As", 33)
        self._put_element_in_table("Selenium", "Se", 34)
        self._put_element_in_table("Bromine", "Br", 35)
        self._put_element_in_table("Krypton", "Kr", 36)
        self._put_element_in_table("Rubidium", "Rb", 37)
        self._put_element_in_table("Strontium", "Sr", 38)
        self._put_element_in_table("Yttrium", "Y", 39)
        self._put_element_in_table("Zirconium", "Zr", 40)
        self._put_element_in_table("Niobium", "Nb", 41)
        self._put_element_in_table("Molybdenum", "Mo", 42)
        self._put_element_in_table("Technetium", "Tc", 43)
        self._put_element_in_table("Ruthenium", "Ru", 44)
        self._put_element_in_table("Rhodium", "Rh", 45)
        self._put_element_in_table("Palladium", "Pd", 46)
        self._put_element_in_table("Silver", "Ag", 47)
        self._put_element_in_table("Cadmium", "Cd", 48)
        self._put_element_in_table("Indium", "In", 49)
        self._put_element_in_table("Tin", "Sn", 50)
        self._put_element_in_table("Antimony", "Sb", 51)
        self._put_element_in_table("Tellurium", "Te", 52)
        self._put_element_in_table("Iodine", "I", 53)
        self._put_element_in_table("Xenon", "Xe", 54)
        self._put_element_in_table("Cesium", "Cs", 55)
        self._put_element_in_table("Barium", "Ba", 56)
        self._put_element_in_table("Lanthanum", "La", 57)
        self._put_element_in_table("Cerium", "Ce", 58)
        self._put_element_in_table("Praseodymium", "Pr", 59)
        self._put_element_in_table("Neodymium", "Nd", 60)
        self._put_element_in_table("Promethium", "Pm", 61)
        self._put_element_in_table("Samarium", "Sm", 62)
        self._put_element_in_table("Europium", "Eu", 63)
        self._put_element_in_table("Gadolinium", "Gd", 64)
        self._put_element_in_table("Terbium", "Tb", 65)
        self._put_element_in_table("Dysprosium", "Dy", 66)
        self._put_element_in_table("Holmium", "Ho", 67)
        self._put_element_in_table("Erbium", "Er", 68)
        self._put_element_in_table("Thulium", "Tm", 69)
        self._put_element_in_table("Ytterbium", "Yb", 70)
        self._put_element_in_table("Lutetium", "Lu", 71)
        self._put_element_in_table("Hafnium", "Hf", 72)
        self._put_element_in_table("Tantalum", "Ta", 73)
        self._put_element_in_table("Tungsten", "W", 74)
        self._put_element_in_table("Rhenium", "Re", 75)
        self._put_element_in_table("Osmium", "Os", 76)
        self._put_element_in_table("Iridium", "Ir", 77)
        self._put_element_in_table("Platinum", "Pt", 78)
        self._put_element_in_table("Gold", "Au", 79)
        self._put_element_in_table("Mercury", "Hg", 80)
        self._put_element_in_table("Thallium", "Tl", 81)
        self._put_element_in_table("Lead", "Pb", 82)
        self._put_element_in_table("Bismuth", "Bi", 83)
        self._put_element_in_table("Polonium", "Po", 84)
        self._put_element_in_table("Astatine", "At", 85)
        self._put_element_in_table("Radon", "Rn", 86)
        self._put_element_in_table("Francium", "Fr", 87)
        self._put_element_in_table("Radium", "Ra", 88)
        self._put_element_in_table("Actinium", "Ac", 89)
        self._put_element_in_table("Thorium", "Th", 90)
        self._put_element_in_table("Protactinium", "Pa", 91)
        self._put_element_in_table("Uranium", "U", 92)
        self._put_element_in_table("Neptunium", "Np", 93)
        self._put_element_in_table("Plutonium", "Pu", 94)
        self._put_element_in_table("Americium", "Am", 95)
        self._put_element_in_table("Curium", "Cm", 96)
        self._put_element_in_table("Berkelium", "Bk", 97)
        self._put_element_in_table("Californium", "Cf", 98)
        self._put_element_in_table("Einsteinium", "Es", 99)
        self._put_element_in_table("Fermium", "Fm", 100)
        self._put_element_in_table("Mendelevium", "Md", 101)
        self._put_element_in_table("Nobelium", "No", 102)
        self._put_element_in_table("Lawrencium", "Lr", 103)
        self._put_element_in_table("Rutherfordium", "Rf", 104)
        self._put_element_in_table("Dubnium", "Db", 105)
        self._put_element_in_table("Seaborgium", "Sg", 106)
        self._put_element_in_table("Bohrium", "Bh", 107)
        self._put_element_in_table("Hassium", "Hs", 108)
        self._put_element_in_table("Meitnerium", "Mt", 109)
        self._put_element_in_table("Darmstadtium", "Ds", 110)
        self._put_element_in_table("Roentgenium", "Rg", 111)
        self._put_element_in_table("Copernicium", "Cn", 112)
        self._put_element_in_table("Nihonium", "Nh", 113)
        self._put_element_in_table("Flerovium", "Fl", 114)
        self._put_element_in_table("Moscovium", "Mc", 115)
        self._put_element_in_table("Livermorium", "Lv", 116)
        self._put_element_in_table("Tennessine", "Ts", 117)
        self._put_element_in_table("Oganesson", "Og", 118)
    def lookup_by_english_name(self, english_name):
        if type(english_name) != str:
            raise TypeError("Expected str. Found: " + str(type(english_name)))
        else:
            return self._by_english_name.get(english_name)
    def lookup_by_symbol(self, symbol):
        if type(symbol) != str:
            raise TypeError("Expected str. Found: " + str(type(symbol)))
        else:
            return self._by_symbol.get(symbol)
    def lookup_by_atomic_number(self, atomic_number):
        try:
            return self._by_atomic_number.get(int(atomic_number))
        except:
            raise TypeError("Expected type that could be cast to int. Found: " + str(type(atomic_number)))

element_lookup_table = ElementLookupTable()