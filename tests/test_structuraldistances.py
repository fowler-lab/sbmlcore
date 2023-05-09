# Tests for StructuralDistances.py
import pandas
import numpy
import pytest
import sbmlcore

def test_chain_selection():

    a = {'segid': ['A', 'A', 'A', 'B', 'C', 'C'], 'mutation': ['I3D','S4K', 'Q5V', 'R6D', 'S450F', 'D435F']}
    df = pandas.DataFrame(a)

    a = sbmlcore.StructuralDistances('tests/5uh6.pdb','resname ZN', 'Zn_distance', offsets = {'A': 0, 'B': 0, 'C': -6})
    df =  a._add_feature(df)

    # these should all fail!

    # Loading the wrong file
    with pytest.raises(AssertionError):
        a = sbmlcore.StructuralDistances('tests/5uh7.pdb','resname ZN', 'Zn_distance', offsets = {'A': 0, 'B': 0, 'C': -6})

    # Distance selection is not a string
    with pytest.raises(AssertionError):
        a = sbmlcore.StructuralDistances('tests/5uh6.pdb', 5, 'Zn_distance', offsets = {'A': 0, 'B': 0, 'C': -6})

    # Offsets is not a dictionary
    with pytest.raises(AssertionError):
        a = sbmlcore.StructuralDistances('tests/5uh6.pdb','resname ZN', 'Zn_distance', offsets = [0, 0, -6])

    # Incorrectly naming a segid offset with an incorrect letter
    with pytest.raises(AssertionError):
        a = sbmlcore.StructuralDistances('tests/5uh6.pdb','resname ZN', 'Zn_distance', offsets = {'A': 0, 'B': 0, 'Z': -6})

    # Incorrectly naming a segid offset with a lower case letter
    with pytest.raises(AssertionError):
        a = sbmlcore.StructuralDistances('tests/5uh6.pdb','resname ZN', 'Zn_distance', offsets = {'A': 0, 'B': 0, 'c': -6})

    # Specifying an offset value as a non-integer e.g. as a string
    with pytest.raises(AssertionError):
        a = sbmlcore.StructuralDistances('tests/5uh6.pdb','resname ZN', 'Zn_distance', offsets = {'A': 0, 'B': 0, 'C': 'c'})

    # Specifying the offset value as a non-integer e.g. as a float
    with pytest.raises(AssertionError):
        a = sbmlcore.StructuralDistances('tests/5uh6.pdb','resname ZN', 'Zn_distance', offsets = {'A': 0, 'B': 0, 'C': 2.3})

    # Joining distances to the wrong residue in the existing mutation dataframe because the offsets for the distance df have been specified incorrectly and therefore do not correspond to those in the mutation df --> will return NaNs
    # Tests to check that code will fail if more than half or half of the column is NaNs
    # Testing fail if half of the distance column is NaNs
    with pytest.raises(AssertionError):
        a = sbmlcore.StructuralDistances('tests/5uh6.pdb','resname MG', 'Mg_distance', offsets = {'A': 1, 'B': 0, 'C': -6})
        df = a._add_feature(df)

    #Testing fail if more than half of the distance column is NaNs
    with pytest.raises(AssertionError):
        a = sbmlcore.StructuralDistances('tests/5uh6.pdb','resname MG', 'Mg_distance', offsets = {'A': 1, 'B': -1, 'C': -6})
        df = a._add_feature(df)
