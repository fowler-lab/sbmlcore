# Tests for StructuralDistances.py
import pandas
import numpy
import pytest
import sbmlcore

def test_freesasa_working():
    b = {'segid': ['A', 'A', 'A'], 'mutation': ['M1D','R2K', 'A3V']}
    df = pandas.DataFrame(b)

    a = sbmlcore.FreeSASA('tests/3pl1.pdb')
    df = a.add_feature(df)

    # These should fail!
    # Loading the wrong file
    with pytest.raises(IOError):
        a = sbmlcore.FreeSASA('tests/3pl2.pdb')
