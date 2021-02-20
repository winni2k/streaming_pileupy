import timeit
import os
import sys
from pathlib import Path

fixture_dir = Path(os.environ["FIXTURES"])

class Suite:
    timer = timeit.default_timer

    def time_head5k(self):
        os.system(f"spileup {fixture_dir}/gawad/gawad_chr20.head5k.sam {fixture_dir}/gawad/sample_names.txt > /dev/null")
