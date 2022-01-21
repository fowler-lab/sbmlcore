# Tests for the '.AminoAcidHydropathyChange' defs for different hydropathy scales
import pathlib
import shutil
import subprocess

def test_stride_ok(tmp_path):

    if pathlib.Path("./stride").exists():
        stride = pathlib.Path("./stride").resolve()

    # or if there is one in the $PATH use that one
    elif shutil.which('stride') is not None:
         stride = pathlib.Path(shutil.which('stride'))

    process = subprocess.Popen(
        [
            stride,
            'tests/3pl1.pdb'
        ],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )

    stdout, stderr = process.communicate()

    # insist that the above command did not fail
    assert process.returncode == 0
