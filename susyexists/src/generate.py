from . import scaffold
from . import reads
from . import utils
import os
import sys


def pw_input(project_id, calculation,config=False, degauss=None, job_id=None, initial_guess=None, k_points=None, poscar=None, layer=None, lattice_constant=False,atomic_positions=False,pseudo=False):
    
    #Check file name, if it is None, uses degauss instead
    if (not job_id):
        if (degauss):
            job_id = degauss
        else:
            # sys.stdout.writelines("No file name \n")
            # print("No file name")
             raise Exception("Define a job_id")
    
    #Default config
    if config==False:
        config = utils.configure_pw()

    #Create directory for input files
    try:
        os.makedirs(f'./Projects/')
    except:
        pass
    try:
        os.makedirs(f'./Projects/{project_id}')
        os.makedirs(f'./Projects/{project_id}/{job_id}')
    except:
        try:
            os.makedirs(f'./Projects/{project_id}/{job_id}')
        except:
            pass

    
    # Set parameters from input
    config["file_path"] = f"./Projects/{project_id}/{job_id}/{calculation}.in"

    #If degauss is given explicitly use it instead
    if(degauss!=None):
        config["system"]["degauss"] = degauss
    config["control"]["outdir"] = f"./Projects/{project_id}/{job_id}/"
    # nat = int(config["system"]['nat'])
    

    

    if calculation == 'vc-relax':
        # Import initial cell and atom parameters
        if poscar != None:
            cell, atoms = reads.read_poscar(f'{poscar}')
        elif initial_guess != None:
            cell, atoms = reads.read_vc_relax(f"{initial_guess}")
        else:
            try:
                cell, atoms = config['cell_parameters'], config['atomic_positions']
            except:
                raise Exception("PLease enter atomic position and lattice constants")

    elif calculation == 'relax':
        if poscar != None:
            cell, atoms = reads.read_poscar(f'{poscar}')
        try:
            cell, atoms = reads.read_vc_relax(f"./Projects/{project_id}/{job_id}/vc-relax.out")
        except:
            try:
                cell, atoms = config['cell_parameters'], config['atomic_positions']
            except:
                raise Exception("PLease enter atomic position and lattice constants")


    else:
        if poscar != None:
            cell, atoms = reads.read_poscar(f'{poscar}')
        try:
            cell, atoms = reads.read_vc_relax(f"./Projects/{project_id}/{job_id}/vc-relax.out")
            try:
                atoms = reads.read_relax(f"./Projects/{project_id}/{job_id}/relax.out")
            except:
                pass
        except:
            try:
                cell, atoms = config['cell_parameters'], config['atomic_positions']
            except:
                raise Exception("PLease enter atomic position and lattice constants")

    if(layer=='mono'):
            atoms = utils.make_monolayer(atoms)

    if k_points != None:
        config['k_points'] = k_points

    if(lattice_constant):
        cell=lattice_constant
    if(atomic_positions):
        atoms=atomic_positions




    #Check for pseudopotentials
    if pseudo==False:
        config['atomic_species']=utils.default_pseudo(atoms)

    config["control"]['prefix'] = job_id
    config["system"]['nat'] = len(atoms)
    #check types of atoms
    try:
        config["system"]['ntyp']
    except:
        config["system"]['ntyp'] = utils.atom_type(atoms)
    
    #check ibrav
    try:
        config["system"]['ibrav']
    except:
        config["system"]['ibrav'] = 0
    config['atomic_positions'] = atoms
    config['cell_parameters'] = cell
    config["control"]['calculation'] = calculation

    if calculation == 'bands-pp':
        scaffold.bands_pp(config)
    else:
        scaffold.pw(config)
    return


def ph_input(project_id,calculation,job_id='default',config=False):
    #Default config
    if config==False:
        config = utils.configure_ph()
    # Set parameters from input
    config["file_path"] = f"./Projects/{project_id}/{job_id}/{calculation}.in"
    if calculation == 'ph':
        config["inputph"]["outdir"] = f"./Projects/{project_id}/{job_id}/"
        config["inputph"]['prefix'] = job_id
        config["inputph"]["fildyn"]= f'./Projects/{project_id}/{job_id}/{job_id}.dyn'
        scaffold.ph(config)
    elif calculation == 'q2r':
        config['input']['fildyn'] = f'./Projects/{project_id}/{job_id}/{job_id}.dyn'
        config['input']['flfrc']  = f'./Projects/{project_id}/{job_id}/{job_id}.fc'
        scaffold.q2r(config)
    elif calculation == 'matdyn':
        config['input']['flfrq'] = f'./Projects/{project_id}/{job_id}/{job_id}.freq'
        config['input']['flfrc']  = f'./Projects/{project_id}/{job_id}/{job_id}.fc'
        scaffold.matdyn(config)
    elif calculation == 'plotband':
        scaffold.plotband(config)
    elif calculation == "ph_plot":
        scaffold.ph_plot(config)


def runner(project_name,iteration,file_names,calculation,qe_path,ncore):
    pre = f'''QE={qe_path}
SLURM_NTASKS={ncore}
project_name={project_name}
    '''
    # print(pre)
    post =f'''
echo Starting {iteration} calculation
for file_name in {" ".join(file_names)}
do
echo $file_name is started
mpirun -n $SLURM_NTASKS $QE/pw.x -inp ./$project_name/$file_name/{calculation}.in > ./$project_name/$file_name/{calculation}.out
echo $file_name is done
done
echo All {iteration} calculations are completed
    '''
    # print(post)
    with open(f"./{iteration}.sh", 'w') as file:
        file.write(pre)
        file.write(post)
    return