# Tests for DSSP secondary structure code
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

    process.wait()

    # insist that the above command did not fail
    assert process.returncode == 0



def test_dssp_ok(tmp_path):

    if pathlib.Path("./mkdssp").exists():
        mkdssp = pathlib.Path("./mkdssp").resolve()

    # or if there is one in the $PATH use that one
    elif shutil.which('mkdssp') is not None:
         mkdssp = pathlib.Path(shutil.which('mkdssp'))

    process = subprocess.Popen(
        [
            mkdssp,
            'tests/3pl1.pdb'
        ],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )

    process.wait()

    # insist that the above command did not fail
    assert process.returncode == 0
