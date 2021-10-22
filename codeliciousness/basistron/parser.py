# -*- coding: utf-8 -*-
"""
A CCCBDB-specific table parser. Borrowed heavily
from https://github.com/marcelo-mason/cccbdb-calculation-parser
"""
from html.parser import HTMLParser

from typing import List, Optional

class TableParser(HTMLParser):
    """Single table parser."""

    def __init__(self):
        super().__init__()
        self.td = False
        self.th = False
        self.tr = False
        self.current_cell = []
        self.current_row = []
        self.table = []

    def handle_starttag(self, tag: str, attrs: List[List[str]]) -> None:
        for t in ["td", "th", "tr"]:
            if tag == t:
                setattr(self, tag, True)
    
    def handle_data(self, data: str) -> None:
        if self.td or self.th:
            self.current_cell.append(data.strip())

    def handle_endtag(self, tag: str) -> None:
        for t in ["td", "th", "tr"]:
            if tag == t:
                setattr(self, tag, False)
        if tag in ["td", "th"]:
            self.current_row.append(self.current_cell[::-1])
            self.current_cell = []
        if tag == "tr":
            self.table.append(self.current_row)
            self.current_row = []

    def pad_table(self, table: List[List[Optional[str]]]) -> List[List[str]]:
        padlen = max((len(row) for row in table))
        for i, row in enumerate(table):
            while len(row) < padlen:
                row.insert(0, "")
            flat = []
            for cell in row:
                if isinstance(cell, str):
                    flat.append(cell)
                    continue
                flat.append("") if not cell else flat.extend(cell)
            table[i] = flat
        return table