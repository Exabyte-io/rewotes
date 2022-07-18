import re
import urllib.request


class ConvTracker:
    """
    Args:
        config (dict): Exabyte API config.
        job_endpoints (JobEndpoints): Exabyte API endpoint.
        kwargs:
            cutoff (float): Desired energy cutoff in eV.
            energy (list): Total energy values. Can be used as a pseudo-restart to convergence.

    Attributes:
        config (dict): Exabyte API config.
        job_endpoints (JobEndpoints): Exabyte API endpoint.
        cutoff (float): Desired energy cutoff in eV.
        energy (list): List of energy values used to check for convergence.
    """

    def __init__(self, config, job_endpoints, cutoff=1e-5, energy=[]):
        self.config = config
        self.job_endpoints = job_endpoints
        self.cutoff = cutoff            # Units = eV
        self.energy = energy            # Array of energies can be passed in to continue a job set.

    def create_submit_job(self, kp, jobs_set=None, job_name_prefix="kpoint"):
        """
        Creates and submits a given job.

        Args:
            kp (int): Value of kpoints. Also used to generate job name.

        kwargs:
            jobs_set (str): ID of job set.
            job_name_prefix (str): Name of job prepended to kpoint value.
        """
        job_name = {"name": f"{job_name_prefix}_{kp}"}
        self.config.update(job_name)
        job = self.job_endpoints.create(self.config)

        if jobs_set is not None:
            self.job_endpoints.move_to_set(job["_id"], "", jobs_set["_id"])

        # Update K-Point Values
        # This is not an ideal way to set kpoints, but the built in convergence tool did npt work as expected, and adjusting the workflow did not update render.
        job["workflow"]["subworkflows"][0]["units"][0]["input"][0]["rendered"] = job["workflow"]["subworkflows"][0]["units"][0]["input"][0]["rendered"].replace("K_POINTS automatic\n10 10 10 0 0 0", f"K_POINTS automatic\n{kp} {kp} {kp} 0 0 0")
        self.job_endpoints.update(job["_id"], job)
        self.job_endpoints.submit(job['_id'])

        return job["_id"]

    def parse_output(self, job_id):
        """
        Read energy from results file

        Args:
            job_id (str): ID of job to get results from.
        """
        files = self.job_endpoints.list_files(job_id)
        output_file = [file for file in files if file["name"] ==  'pw_scf.out'][0]
        server_response = urllib.request.urlopen(output_file['signedUrl'])
        output_file_bytes = server_response.read()
        output_file = output_file_bytes.decode(encoding="UTF-8")
        output_as_array = output_file.split("\n")
        total_energy_ry = float(re.split(" +", [row for row in output_as_array if "!    total energy" in row][0])[-2])
        total_energy_ev = total_energy_ry * 13.6056980659

        return total_energy_ev

    def check_convergence(self):
        """
        Check if energy convergence reached.
        """
        if len(self.energy) < 2:
            return False
        else:
            return abs(self.energy[-1] - self.energy[-2]) <= self.cutoff

