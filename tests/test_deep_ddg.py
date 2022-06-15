# Tests for DeepDDG.py

import pandas
import pytest
import sbmlcore

# let's be tough and use the RNAP structure, 5uh6
# note that B:S4V and F:A3W are not resolved in the PDB and
# the PDB has an offset for chain C where all the reported resids are +6

def test_missing_file():

    a = {'segid': ['A', 'A', 'B', 'B', 'C', 'C', 'C', 'D', 'D', 'E', 'F', 'F'], 'mutation': ['I3F','Q5T', 'S10Q', 'S4V', 'S450L', 'H445D', 'I491F', 'I491V', 'G594E', 'A101E', 'E245D', 'A3W']}
    df = pandas.DataFrame(a)

    with pytest.raises(AssertionError):
        a = sbmlcore.DeepDDG('tests/5uh7-protein-chains.ddg', offsets = {'A': 0, 'B': 0, 'C':-6, 'D':0, 'E':0, 'F':0})


def test_bad_file_1():

    a = {'segid': ['A', 'A', 'B', 'B', 'C', 'C', 'C', 'D', 'D', 'E', 'F', 'F'], 'mutation': ['I3F','Q5T', 'S10Q', 'S4V', 'S450L', 'H445D', 'I491F', 'I491V', 'G594E', 'A101E', 'E245D', 'A3W']}
    df = pandas.DataFrame(a)

    with pytest.raises(OSError):
        a = sbmlcore.DeepDDG('tests/3pl1-bad.ddg', offsets = {'A': 0})


def test_bad_offsets_1():

    a = {'segid': ['A', 'A', 'B', 'B', 'C', 'C', 'C', 'D', 'D', 'E', 'F', 'F'], 'mutation': ['I3F','Q5T', 'S10Q', 'S4V', 'S450L', 'H445D', 'I491F', 'I491V', 'G594E', 'A101E', 'E245D', 'A3W']}
    df = pandas.DataFrame(a)

    with pytest.raises(AssertionError):
        a = sbmlcore.DeepDDG('tests/5uh6-protein-chains.ddg', offsets = {'A': 0.5, 'B': 0, 'C':56.2, 'D':0, 'E':0, 'F':0})


def test_bad_offsets_2():

    a = {'segid': ['A', 'A', 'B', 'B', 'C', 'C', 'C', 'D', 'D', 'E', 'F', 'F'], 'mutation': ['I3F','Q5T', 'S10Q', 'S4V', 'S450L', 'H445D', 'I491F', 'I491V', 'G594E', 'A101E', 'E245D', 'A3W']}
    df = pandas.DataFrame(a)

    with pytest.raises(AssertionError):
        a = sbmlcore.DeepDDG('tests/5uh6-protein-chains.ddg', offsets = {'A': 0.0, 'B': 0, 'C':-6.0, 'D':0, 'E':0, 'F':0})


def test_bad_offsets_3():

    a = {'segid': ['A', 'A', 'B', 'B', 'C', 'C', 'C', 'D', 'D', 'E', 'F', 'F'], 'mutation': ['I3F','Q5T', 'S10Q', 'S4V', 'S450L', 'H445D', 'I491F', 'I491V', 'G594E', 'A101E', 'E245D', 'A3W']}
    df = pandas.DataFrame(a)

    with pytest.raises(AssertionError):
        a = sbmlcore.DeepDDG('tests/5uh6-protein-chains.ddg', offsets = [1,2,1,3,2,1])

def test_correct_offsets():

    a = {'segid': ['A', 'A', 'B', 'B', 'C', 'C', 'C', 'D', 'D', 'E', 'F', 'F'], 'mutation': ['I3F','Q5T', 'S10Q', 'S4V', 'S450L', 'H445D', 'I491F', 'I491V', 'G594E', 'A101E', 'E245D', 'A3W']}
    df = pandas.DataFrame(a)

    a = sbmlcore.DeepDDG('tests/5uh6-protein-chains.ddg', offsets = {'A': 0, 'B': 0, 'C':-6, 'D':0, 'E':0, 'F':0})

    df2 = a._add_feature(df)

    assert len(df2) == 12

    # only two mutations should be missing a ddG since they are not resolved in the PDB
    # this means the join was correct so the offset is working correctly
    assert df2.deep_ddG.isna().sum() == 2

    # let's explicitly check rpoB S450L
    mask = (df2.segid=='C') & (df2.mutation=='S450L')

    assert (df2[mask].deep_ddG < 0).all()


def test_incorrect_offsets():

    a = {'segid': ['A', 'A', 'B', 'B', 'C', 'C', 'C', 'D', 'D', 'E', 'F', 'F'], 'mutation': ['I3F','Q5T', 'S10Q', 'S4V', 'S450L', 'H445D', 'I491F', 'I491V', 'G594E', 'A101E', 'E245D', 'A3W']}
    df = pandas.DataFrame(a)

    a = sbmlcore.DeepDDG('tests/5uh6-protein-chains.ddg', offsets = {'A': 1, 'B': 1, 'C':0, 'D':1, 'E':1, 'F':1})

    df2 = a._add_feature(df)

    assert len(df2) == 12

    # only two mutations should be missing a ddG since they are not resolved in the PDB
    # this means the join was correct so the offset is working correctly
    assert df2.deep_ddG.isna().sum() == 11

    # let's explicitly the mutation in chain F that does match since both 244 and 245 are Glu
    mask = (df2.segid=='F') & (df2.mutation=='E245D')

    assert (df2[mask].deep_ddG < 0).all()
