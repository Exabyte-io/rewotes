import sys
import os

directory = sys.argv[1]
for root, dirs, files in os.walk(directory, topdown=False):
    for f in files:
        filepath = os.path.join(root, f)
        if os.path.exists(filepath):
            if f == 'POSCAR':
                pass
            elif f == 'INCAR': 
                pass
            elif f == 'POTCAR':
                pass
            elif f == 'KPOINTS': 
                pass
            elif f == 'submit.sh':
                pass
            elif f == 'run.py': 
                pass
            else:
                os.remove(filepath) 
