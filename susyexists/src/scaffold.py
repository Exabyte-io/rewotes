import numpy as np
import pandas as pd
from . import reads
from . import writes
from . import utils
import json

def pw(parameters):
    with open(parameters['file_path'], 'w') as file:
        for i,j in parameters.items():
            try:
                if(type(j)==dict):
                    file.write(f"&{i.upper()} \n")
                for k,l in j.items():
                    try:
                        float(l)
                        file.write(f"{k} = {l} \n")
                    except:
                        try:
                            json.loads(l)
                            file.write(f"{k} = .{l}. \n")
                        except:
                            if(bool(l)):
                                file.write(f"{k} = '{l}' \n")
                file.write("/ \n")
            except:
                pass
    writes.write_atom_species(parameters['file_path'], parameters['atomic_species'])
    writes.write_atom_positions(parameters['file_path'],parameters['atomic_positions'])
    writes.write_cell_parameters(parameters['file_path'], parameters['cell_parameters'])
    if parameters["control"]['calculation']=='bands':
        writes.write_k_points_bands(parameters['file_path'], parameters['k_points_bands'])
    else:
        writes.write_k_points(parameters['file_path'], parameters['k_points'])
    return


def bands_pp(parameters):
    filedata = f"""
&bands
outdir = '{parameters["control"]['outdir']}'
prefix = '{parameters["control"]['prefix']}'
filband = '{parameters["control"]['outdir']}bands.dat'
/
"""
    with open(parameters['file_path'], 'w') as file:
        file.write(filedata)
    return


def ph(parameters):
    with open(parameters['file_path'], 'w') as file:
        for i,j in parameters.items():
            try:
                if(type(j)==dict):
                    file.write(f"&{i.upper()} \n")
                for k,l in j.items():
                    try:
                        float(l)
                        file.write(f"{k} = {l} \n")
                    except:
                        try:
                            json.loads(l)
                            file.write(f"{k} = .{l}. \n")
                        except:
                            if(bool(l)):
                                file.write(f"{k} = '{l}' \n")
                file.write("/ \n")
            except:
                pass
    return

def q2r(parameters):
    with open(parameters['file_path'], 'w') as file:
        for i,j in parameters.items():
            try:
                if(type(j)==dict):
                    file.write(f"&{i.upper()} \n")
                for k,l in j.items():
                    try:
                        float(l)
                        file.write(f"{k} = {l} \n")
                    except:
                        try:
                            json.loads(l)
                            file.write(f"{k} = .{l}. \n")
                        except:
                            if(bool(l)):
                                file.write(f"{k} = '{l}' \n")
                file.write("/ \n")
            except:
                pass
    return

def matdyn(parameters):
    with open(parameters['file_path'], 'w') as file:
        for i,j in parameters.items():
            try:
                if(type(j)==dict):
                    file.write(f"&{i.upper()} \n")
                for k,l in j.items():
                    try:
                        float(l)
                        file.write(f"{k} = {l} \n")
                    except:
                        try:
                            json.loads(l)
                            file.write(f"{k} = .{l}. \n")
                        except:
                            if(bool(l)):
                                file.write(f"{k} = '{l}' \n")
                file.write("/ \n")
            except:
                pass
    writes.write_k_points_matdyn(parameters['file_path'], parameters['k_points_bands'])
    return

def plotband(parameters):
    filedata = f"""
{parameters['prefix']}.freq
0 550
{parameters['prefix']}.xmgr
"""
    with open(f"./dyn/{parameters['prefix']}/{parameters['prefix']}-plotband.in", 'w') as file:
        file.write(filedata)
    return