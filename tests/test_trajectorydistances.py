import pandas
import pytest
import sbmlcore

def test_missing_file():

    with pytest.raises(AssertionError):
        a = sbmlcore.TrajectoryDistances('tests/3pl1.pdb', 'tests/missing.xtc', 'resname FE', 'dist_FE')

def test_wrong_type():

    with pytest.raises(AssertionError):
        a = sbmlcore.TrajectoryDistances('tests/rpob-5uh6-3-warm.gro.gz', ['tests/rpob-5uh6-3-md-1-50ns-dt10ns-nojump.xtc'], 'resname RFP', 'dist_RIF', distance_type='average')

def test_not_list():

    with pytest.raises(AssertionError):
        a = sbmlcore.TrajectoryDistances('tests/3pl1.pdb', 'tests/3pl1.xtc', 'resname FE', 'dist_FE')

def test_not_list_1():

    with pytest.raises(AssertionError):
        a = sbmlcore.TrajectoryDistances('tests/3pl1.pdb', ['tests/3pl1.xtc'], 'resname FE', 'dist_FE', start_time=-1, end_time=12.2)

def test_not_list_2():

    with pytest.raises(AssertionError):
        a = sbmlcore.TrajectoryDistances('tests/3pl1.pdb', ['tests/3pl1.xtc'], 'resname FE', 'dist_FE', start_time=1, end_time='b')

def test_not_list():

    with pytest.raises(AssertionError):
        a = sbmlcore.TrajectoryDistances('tests/3pl1.pdb', ['tests/3pl1.xtc'], 'resname FE', 'dist_FE', start_time=20, end_time=10.0)


def test_working():

    a = sbmlcore.TrajectoryDistances('tests/rpob-5uh6-3-warm.gro.gz', ['tests/rpob-5uh6-3-md-1-50ns-dt10ns-nojump.xtc'], 'resname RFP', 'dist_RIF')
