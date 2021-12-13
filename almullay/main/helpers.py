'''
    Here we include functions to read input data from .csv,
    collect data from MAPI http(s) requests, and clean prior
    to pre-processing and machine learning in core.py
'''

from numpy import zeros
from pandas import read_csv, to_numeric, DataFrame

import os, glob, sys, importlib, urllib.request
from utils.generic import display_JSON
import settings; importlib.reload(settings)
from settings import ENDPOINT_ARGS, MATERIALS_PROJECT_API_KEY

from exabyte_api_client.endpoints.materials import MaterialEndpoints


def sheetData():
    '''
    Read and clean data from EXABYTE.IO-HT-BG-2018-PIII-FINAL - Phase III.csv
    (from "Electronic properties of binary compounds with high fidelity and
    high throughput" by Das and Bazhirov)
    '''

    Y = read_csv('../data/EXABYTE.IO-HT-BG-2018-PIII-FINAL - Phase III.csv')

    Y.rename(columns={ Y.columns[5]: "GGAin", Y.columns[7]: "HSEin" }, inplace = True) 

    # Remove extra info and labels
    Y.drop(Y.columns[[2,3,8,9,10]], axis=1, inplace=True)
    Y.drop(Y.index[:2], inplace=True) 

    # If GGA or HSE calculations are missing (or are 'zero'), remove row
    DataFrame.dropna(Y, axis=0, how='any', thresh=None, subset=None, inplace=True)
    
    Y = Y.apply(to_numeric, errors='ignore', downcast='float')
    Y.drop(Y[(Y.GGA == 0.0) | (Y.GGAin == 0.0) | (Y.HSE == 0.0) | \
        (Y.HSEin == 0.0)].index, inplace=True)

    # Target function will be the correction (difference) between GGA and HSE columns

    return Y

def mapiData(MATERIALS_PROJECT_IDS):
    '''
        Given list of mapi IDs, obtain material properties and construct
        feature matrix
    '''

    endpoint = MaterialEndpoints(*ENDPOINT_ARGS)
    materials = endpoint.import_from_materialsproject(MATERIALS_PROJECT_API_KEY, \
        MATERIALS_PROJECT_IDS)

    # Pssible element encodings:
    # 1. Two components specify proton number and valence shell occupancy of 
    # element to uniquely specify it (to assist in generalizability to elements
    # not seen by learning algorithm). Three components specify absolute position
    # of element in lattice.
    # 2. Enumerate elements and specify by index, with no extra chemical info (no
    # generalizability to unseen elements).

    # We choose approach (1.) out of interest in generalizability.

    # Six components specify elements in binary structure: atomic number of each
    # element, number of valence electrons of each element (periodic table grouping),
    # and the quantity of each in the composition.
    # Six numbers specify the lattice structure (three lattice separations and three 
    # angles).

    from re import findall, split
    from pymatgen.core.periodic_table import Element

    X = zeros((len(MATERIALS_PROJECT_IDS), 6 + 6))
    row = 0
    for m in materials:
        X[row, 6:] = [
                        m['lattice']['a'],
                        m['lattice']['b'],
                        m['lattice']['c'],
                        m['lattice']['alpha'],
                        m['lattice']['beta'],
                        m['lattice']['gamma'],
                    ]
        composition = findall('[A-Z][^A-Z]*', m['formula'])
        parse_comp = [split(r'(^[^\d]+)', comp)[1:] for comp 
                        in composition]
        elements = [Element(e[0]) for e in parse_comp]

        # Easy to generalize beyond binary compositions with a loop
        X[row, 0:2] = [elements[0].group, elements[0].Z]
        X[row, 2] = int(parse_comp[0][1])
        X[row, 3:7] = [elements[1].group, elements[1].Z]
        X[row, 7] = int(parse_comp[1][1]) 
  
    return X


