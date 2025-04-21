import time
import os

def wait_qe(output_file: str) -> None:
    """Waits for QE calculation to finish"
 
        output_file: output file of point calc
    """

    finished = False
    while not finished:
        if os.path.exists(output_file):
            with open(output_file, "r") as file:
                lines = file.readlines()
                
            for line in reversed(lines): 
                if ("JOB DONE." in line):
                    finished = True
                    break
        time.sleep(5)