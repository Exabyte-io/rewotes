import pytest



def test_loading():
    from ConvergenceTracker import kpoint_scheduler_uniform, mae
    from ConvergenceTracker import DriverQE, ConvergenceTrackerQE
    from ConvergenceTracker import calculator_qe, calculator_qe_par

from ConvergenceTracker import kpoint_scheduler_uniform, mae
from ConvergenceTracker import DriverQE, ConvergenceTrackerQE
from ConvergenceTracker import calculator_qe, calculator_qe_par



def test_ksh():
    """Test of sequential kpoint scheduler"""
    k_sch = kpoint_scheduler_uniform(2, 30, 2)
    for el in range(2, 30, 2):
        assert k_sch.get_next() == el



def test_mae_float():
    """Test of mae metrix"""
    assert mae()(5, 10) == 5


def test_mae_vec():
    """Test of mae metrix"""
    assert int(mae()([2, 4], [4, 8])) == 3


    
path_to_exec = "/pscratch/sd/v/vladygin/tools/q-e_new/q-e/bin"


def test_driver():
    """Test of sequential sequential driver on Si"""
    calculator = calculator_qe(path_to_exec)
    driver = DriverQE(workdir = './reference', 
                      calculator = calculator,
                      input_file_name = "pw.in", 
                      encut = 40)

    driver.gen_input(4)
    driver.calc()
    total_energy = driver.extract_target('total_energy')
    assert int(total_energy // 100) == -3


def test_convergence_tracker_Si_seq():
    """Test of sequential convergence tracker on Si"""
    calculator = calculator_qe(path_to_exec)
    
    driver = DriverQE(workdir = './reference', 
                      calculator = calculator,
                      input_file_name = "pw.in", 
                      encut = 40)
    
    tracker = ConvergenceTrackerQE('./reference', 'total_energy', 
                                    driver = driver)
    
    k_list, target_list, errors = tracker.find_opt()

    assert errors[-1] <= 1e-2 or k_list[-1] == 28 

def test_convergence_tracker_Si_par():
    """Test of parallel convergence tracker on Si"""
    calculator = calculator_qe_par(path_to_exec, 1, 1)

    driver = DriverQE(workdir = './reference',
                      calculator = calculator,
                      input_file_name = "pw.in", 
                      encut = 40)
    
    tracker = ConvergenceTrackerQE('./reference', 'total_energy', 
                                    driver = driver)
    
    k_list, target_list, errors = tracker.find_opt()

    assert errors[-1] <= 1e-2 or k_list[-1] == 28 


def test_convergence_tracker_Si_par_job():
    """Test of parallel convergence tracker on Si"""
    calculator = calculator_qe_par(path_to_exec, 1, 1, job = True, workdir = './reference')

    driver = DriverQE(workdir = './reference',
                      calculator = calculator,
                      input_file_name = "pw.in", 
                      encut = 40)
    
    tracker = ConvergenceTrackerQE('./reference', 'total_energy', 
                                    driver = driver)
    
    k_list, target_list, errors = tracker.find_opt()

    assert errors[-1] <= 1e-2 or k_list[-1] == 28 
