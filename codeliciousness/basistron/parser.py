# -*- coding: utf-8 -*-
"""
A CCCBDB-specific table parser. Borrowed heavily
from https://github.com/marcelo-mason/cccbdb-calculation-parser
"""
from html.parser import HTMLParser
from collections import defaultdict

from typing import List, Optional
import pandas as pd
import numpy as np

from basistron import utils

log = utils.get_logger(__name__)

class TableParser(HTMLParser):
    """Single table parser."""

    def __init__(self):
        super().__init__()
        self.td = False
        self.th = False
        self.tr = False
        self.caption = False
        self.title = ""
        self.current_cell = []
        self.current_row = []
        self.table = []

    def handle_starttag(self, tag: str, attrs: List[List[str]]) -> None:
        for t in ["td", "th", "tr", "caption"]:
            if tag == t:
                setattr(self, tag, True)
    
    def handle_data(self, data: str) -> None:
        if self.td or self.th:
            self.current_cell.append(data.strip())
        if self.caption:
            self.title += data.strip()

    def handle_endtag(self, tag: str) -> None:
        for t in ["td", "th", "tr", "caption"]:
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

    def to_df(self):
        """This is where it gets messy."""
        padded = self.pad_table(self.table)
        try:
            df = pd.DataFrame(padded[1:], columns=padded[0])
        except (ValueError, TypeError):
            return None

        def clean_values(df: pd.DataFrame) -> pd.DataFrame:
            return df.replace(r'^\s*$', np.nan, regex=True)

        def clean_columns(df: pd.DataFrame) -> pd.DataFrame:
            index = None
            if df.columns[:2].tolist() == ["", ""]:
                index = ["theory", "implementation"]
                df.columns = index + df.columns[2:].tolist()
            unique_columns = []
            seen = defaultdict(int)
            for column in df.columns:
                seen[column] += 1
                if seen[column] > 1:
                    unique_columns.append(f"{column}{seen[column]}")
                else:
                    unique_columns.append(column)
            df.columns = unique_columns
            return df, index

        def clean_index(df, index):
            if index is not None:
                df.set_index(index, inplace=True)
            try:
                df.drop(("", ""), inplace=True)
            except Exception as e:
                log.error(f"cleaning index failed: {repr(e)}")
            return df.droplevel(0) if index else df

        df = clean_values(df)
        df, index = clean_columns(df)
        df = clean_index(df, index)
        if df.index.duplicated().sum():
            log.warning("found duplicated index entries")

        warn_threshold = len(df.index) // 2
        for column in df.columns:
            df[column] = pd.to_numeric(df[column], errors="coerce")
            nulls = df[column].isnull().sum()
            if nulls > warn_threshold:
                log.warning(f"column {column} has {nulls} nulls")
        df.name = self.title
        return df
