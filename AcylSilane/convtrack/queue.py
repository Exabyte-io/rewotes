import subprocess
import re


class Queue(object):
    def __init__(self):
        pass

    def qstat(self):
        raw_qstat = subprocess.check_output("qstat").split("\n")

        # Format is Headers, followed by a series of dashed lines, followed by actual data
        keys = re.split("\s{2,}", raw_qstat[0])
        vals = []
        for index, line in enumerate(raw_qstat):
            if index <= 1:
                # Skip header and decorative dashes
                continue
            elif re.match("^\s*$", line):
                # Skip blank lines
                continue
            else:
                vals.append(re.split("\s{2,}", line.strip()))
        return [dict(zip(keys, val)) for val in vals]
