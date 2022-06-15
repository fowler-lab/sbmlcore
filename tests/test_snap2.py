# Tests for SNAP2 class in ExternalCode.py
import pandas
import numpy
import pytest
import sbmlcore

#Loading the wrong .csv file
def snap2_loaded():
    with pytest.raises(IOError):
        file = sbmlcore.SNAP2('tests/3pl2-snap2.csv')

def csv_interpret():
    file = sbmlcore.SNAP2('tests/3pl1-snap2.csv')
    b = {'segid': ['A', 'A', 'A'], 'mutation': ['M1D','R2K', 'A3V']}
    df = pandas.DataFrame(b)
    snap2_df=file._add_feature(df)

    #Check that mutation has been correctly assigned and aligned

def snap2_offsets():
    file = sbmlcore.SNAP2('tests/5uh6-snap2.csv')
    b = {'segid': ['A', 'A', 'A', 'B', 'C', 'C'], 'mutation': ['I3D','S4K', 'Q5V', 'R6D', 'S450F', 'D435F']}
    df = pandas.DataFrame(b)
    df = a._add_feature(df)

#Check that dictionary specification is correct
    #Should all fail!
    # Offsets is not a dictionary
    with pytest.raises(AssertionError):
        a = sbmlcore.SNAP2('tests/5uh6-snap2.csv', offsets = [0, 0, -6])

    # Incorrectly naming a segid offset with an incorrect letter
    with pytest.raises(AssertionError):
        a = sbmlcore.SNAP2('tests/5uh6-snap2.csv', offsets = {'A': 0, 'B': 0, 'Z': -6})
        df = a._add_feature(df)

    # Incorrectly naming a segid offset with a lower case letter
    with pytest.raises(AssertionError):
        a = sbmlcore.SNAP2('tests/5uh6-snap2.csv', offsets = {'A': 0, 'B': 0, 'c': -6})
        df = a._add_feature(df)

    # Specifying an offset value as a non-integer e.g. as a string
    with pytest.raises(AssertionError):
        a = sbmlcore.SNAP2('tests/5uh6-snap2.csv', offsets = {'A': 0, 'B': 0, 'C': 'c'})
        df = a._add_feature(df)

    # Specifying the offset value as a non-integer e.g. as a float
    with pytest.raises(AssertionError):
        a = sbmlcore.SNAP2('tests/5uh6-snap2.csv', offsets = {'A': 0, 'B': 0, 'C': 2.3})
        df = a._add_feature(df)


    #Check that offsets correctly align
