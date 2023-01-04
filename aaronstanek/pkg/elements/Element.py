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
