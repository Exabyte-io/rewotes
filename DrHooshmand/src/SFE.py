import numpy as np
import copy
import math
import argparse
from os import listdir
from os.path import isfile, join
import os
from pos import POSCAR
import shutil




def gen_tilt(file,range_shift=([0,1.5],[0,1.5]), y_disp=False, last_del=(False,False), neg= (False, False), no_shift=None, d_shift=None, chain_fold= "POSCARs_O", no_x=1,no_y=1  ):

    '''
    Routine for generating input files for generalized stacking fault energy (GSFE) surface calculation:
    Fault is applied along XY plane, and Z in the input file should be directed normal to the fault plane

    :param file: input file
    :param range_shift: range of shifts to be applied with respect to the planar XY defect cell dimension normalized between (0,1)
    :param y_disp: If True, displacement is applied on the direction prependicular to the Burgers Vector
    :param last_del: If True, last file (which is the periodic image of initial input) will be deleted
    :param neg: if True, displacement is applied in the negative direction
    :param no_shift: if not None, this number of shifts are automatically applied
    :param d_shift: If not None, customized shift displacements are applied
    :param chain_fold: Name of the folder to output the generated files
    :param no_x: Prediocity of the supercell along X provided in the input file (wrt to the conventional cell)
    :param no_y: Prediocity of the supercell along Y provided in the input file (wrt to the conventional cell)

    :return: Generated POSCARs with the corresponding GSFE shifts can be found in the "chain_fold"
    '''

    sh = POSCAR(file)  # Read the input file

    llat_x = list(sh.pos[0][:] / no_x) # single periodicity of the conventional cell along X (direction of Burgers vector)
    llat_y = list(sh.pos[1][:] / no_y) # single periodicity of the conventional cell along Y (direction prependicular to the Burgers vector)

    shift = []
    shift_dir = []


    range_shift_x, range_shift_y = range_shift
    last_del_x , last_del_y = last_del
    neg_x , neg_y = neg


    # Generate the list of shifts to be applied to the cell

    if no_shift is not None:
        no_shift_x, no_shift_y = no_shift
        ar_x = np.arange(range_shift_x[0], range_shift_x[1] + range_shift_x[1] / (no_shift_x - 1),
                         range_shift_x[1] / (no_shift_x - 1))
        ar_x = ar_x.tolist()
        ar_x = [x for x in ar_x if x <= range_shift_x[1]]
        if neg_x:
            ar_x = [-x for x in ar_x]

        if last_del_x:
            del ar_x[-1]

        if not y_disp:

            for i in ar_x:
                shift_dir.append((i, 0, 0))
                s = i * np.array(llat_x)
                shift.append((s[0], s[1], s[2]))

        elif y_disp:
            ar_y = np.arange(range_shift_y[0], range_shift_y[1] + range_shift_y[1] / (no_shift_y - 1),
                             range_shift_y[1] / (no_shift_y - 1))
            ar_y = ar_y.tolist()
            ar_y = [y for y in ar_y if y <= range_shift_y[1]]
            if neg_y:
                ar_y=[-y for y in ar_y]

            if last_del_y:
                del ar_y[-1]

            for i in ar_x:
                for j in ar_y:
                    shift_dir.append((i, j, 0))
                    s = i * np.array(llat_x) + j * np.array(llat_y)
                    shift.append((s[0], s[1], s[2]))

    elif d_shift is not None:
        d_shift_x, d_shift_y = d_shift
        ar_x = np.arange(range_shift_x[0], range_shift_x[1] + d_shift_x, d_shift_x)
        ar_x = ar_x.tolist()
        ar_x = [x for x in ar_x if x <= range_shift_x[1]]
        if neg_x:
            ar_x = [-x for y in ar_x]

        if last_del_x:
            del ar_x[-1]

        if not y_disp:

            for i in ar_x:
                shift_dir.append((i, 0, 0))
                s = i * np.array(llat_x)
                shift.append((s[0], s[1], s[2]))

        elif y_disp:
            ar_y = np.arange(range_shift_y[0], range_shift_y[1] + d_shift_y, d_shift_y)
            ar_y = ar_y.tolist()
            ar_y = [x for x in ar_y if x <= range_shift_y[1]]
            if neg_y:
                ar_y=[-y for y in ar_y]

            if last_del_y:
                del ar_y[-1]

            for i in ar_x:
                for j in ar_y:
                    shift_dir.append((i, j, 0))
                    s = i * np.array(llat_x) + j * np.array(llat_y)
                    shift.append((s[0], s[1], s[2]))
    elif no_shift is None and d_shift is None:
        raise NameError('Define shift parameters')



    # Shift amounts are stored in "shift" list. Now we itereate over the list and create the corresponding POSCARs
    shift = list(shift)

    ### Translate
    for i,(disp, disp_dir) in enumerate(zip(shift,shift_dir)):
        disp_x, disp_y, disp_z = disp



        sh = POSCAR(file)
        sh.pos[2][0] += disp_x
        sh.pos[2][1] += disp_y

        if not os.path.exists(chain_fold):
            os.mkdir(chain_fold)
        out_name = chain_fold+"/"+"POSCAR_"+str(round(i,3))


        sh.wPOSCAR(out_name, str(disp_dir[0])+"_"+str(disp_dir[1]))

    return shift_dir


def clear_fold(fold):

    cdir = os.getcwd()
    if os.path.isdir(os.path.join(cdir, fold)):
        shutil.rmtree(join(cdir,fold))
        os.mkdir(os.path.join(cdir,fold))
    else:
        os.mkdir(os.path.join(cdir,fold))

if __name__ =="__main__":


    ch_fold = "POSCARs_tilt"  #Chain folder to write the files
    cdir = os.getcwd()
    os.chdir(join(cdir,"../example") )
    cdir = os.getcwd()
    print(cdir)
    clear_fold(ch_fold)

    file = "POSCAR_UC_cartesian_basal_3X3_cart" #input file (should be provided in Cartesian coordinates)

    ss = gen_tilt(file, no_shift=(4,4), range_shift=([0, 1], [0, 1]), neg=(False, True),
                  y_disp=True, chain_fold=ch_fold, no_x=3, no_y=3)  # SF

    print("done")




