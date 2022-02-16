# Tests for the '.AminoAcidHydropathyChange' defs for different hydropathy scales
import pandas
import numpy
import pytest
import sbmlcore

def test_chain_selection():

    a = {'segid': ['A', 'A', 'A', 'B', 'C', 'C'], 'mutation': ['I3D','S4K', 'Q5V', 'R6D', 'S450F', 'D435F']}
    df = pandas.DataFrame(a)

    a = sbmlcore.StructuralDistances('tests/5uh6.pdb','resname ZN', 'Zn_distance', offsets = {'A': 0, 'B': 0, 'C': -6})
    df =  a.add_feature(df)

    # these should all fail!
    with pytest.raises(AssertionError):
        a = sbmlcore.StructuralDistances('tests/5uh6.pdb','resname ZN', 'Zn_distance', offsets = {'A': 0, 'B': 0, 'Z': -6})

    with pytest.raises(AssertionError):
        a = sbmlcore.StructuralDistances('tests/5uh6.pdb','resname ZN', 'Zn_distance', offsets = {'A': 0, 'B': 0, 'c': -6})

    with pytest.raises(AssertionError):
        a = sbmlcore.StructuralDistances('tests/5uh7.pdb','resname ZN', 'Zn_distance', offsets = {'A': 0, 'B': 0, 'C': -6})

    with pytest.raises(AssertionError):
        a = sbmlcore.StructuralDistances('tests/5uh6.pdb','resname ZP', 'Zn_distance', offsets = {'A': 0, 'B': 0, 'C': -6})
