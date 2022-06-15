# Tests for FreeSASA class in ExternalCode.py
import pandas
import numpy
import pytest
import sbmlcore

def test_freesasa_rnap():
#    b = {'segid': ['A', 'A', 'A'], 'mutation': ['M1D','R2K', 'A3V']}
    b = {'segid': ['A', 'A', 'A', 'B', 'C', 'C'], 'mutation': ['I3D','S4K', 'Q5V', 'R6D', 'S450F', 'D435F']}
    df = pandas.DataFrame(b)

    a = sbmlcore.FreeSASA('tests/5uh6.pdb', offsets = {'A': 0, 'B': 0, 'C': -6})
    df = a._add_feature(df)

    # These should fail!
    # Loading the wrong file
    with pytest.raises(IOError):
        a = sbmlcore.FreeSASA('tests/5uh7.pdb')

    # Offsets is not a dictionary
    with pytest.raises(AssertionError):
        a = sbmlcore.FreeSASA('tests/5uh6.pdb', offsets = [0, 0, -6])

    # Incorrectly naming a segid offset with an incorrect letter
    with pytest.raises(AssertionError):
        a = sbmlcore.FreeSASA('tests/5uh6.pdb', offsets = {'A': 0, 'B': 0, 'Z': -6})
        df = a._add_feature(df)

    # Incorrectly naming a segid offset with a lower case letter
    with pytest.raises(AssertionError):
        a = sbmlcore.FreeSASA('tests/5uh6.pdb', offsets = {'A': 0, 'B': 0, 'c': -6})
        df = a._add_feature(df)

    # Specifying an offset value as a non-integer e.g. as a string
    with pytest.raises(AssertionError):
        a = sbmlcore.FreeSASA('tests/5uh6.pdb', offsets = {'A': 0, 'B': 0, 'C': 'c'})
        df = a._add_feature(df)

    # Specifying the offset value as a non-integer e.g. as a float
    with pytest.raises(AssertionError):
        a = sbmlcore.FreeSASA('tests/5uh6.pdb', offsets = {'A': 0, 'B': 0, 'C': 2.3})
        df = a._add_feature(df)

#Tests here were for when the testing columns for chain_offsets and pdb_resid had not been dropped
#def test_freesasa_pnca():
#    file = sbmlcore.FreeSASA('tests/3pl1.pdb')
#    b = {'segid': ['A', 'A', 'A'], 'mutation': ['M1D','R2K', 'A3V']}
#    df = pandas.DataFrame.from_dict(b)
#    sasa_df=file._add_feature(df)

    #Check that all entries in the chain_offsets column are zeros because no offset has been specified
    #assert (sasa_df["chain_offsets"] == 0).all(), "All entries in column 'chain_offsets' should be zero as no offset has been specified!"

    #Check that all entries in the resid and pdb_resid columns are the same
    #assert (sasa_df["resid"] == sasa_df["pdb_resid"]).all(), "Columns 'resid' and 'pdb_resid' should be equal as no offset has been specified!"
