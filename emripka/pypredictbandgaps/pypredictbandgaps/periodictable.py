from urllib.request import urlopen
import json

class PeriodicTable:
    """
    Class which houses the Exabyte.io periodic table json file, as well as a list
    of elemental symbols and a map from symbol to element name.
    """
    def __init__(self):
        # use period table from exabyte.io github:
        self.periodic_table_url = urlopen(
            "https://raw.githubusercontent.com/Exabyte-io/periodic-table.js/master/periodic-table.json"
        )
        self.periodic_table = json.loads(self.periodic_table_url.read())
        self.symbol_to_element_map = { details["symbol"]:element for (element, details) in self.periodic_table.items()}
        self.symbols = list(self.symbol_to_element_map.keys())
