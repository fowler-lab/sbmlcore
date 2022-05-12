import pandas
import numpy
import pytest
import sbmlcore


def test_missing_file():

    with pytest.raises(AssertionError):
        a = sbmlcore.TrajectoryDistances(
            "tests/3pl1.pdb",
            "tests/missing.xtc",
            "tests/3pl1.pdb",
            "resname FE",
            "dist_FE",
        )


def test_wrong_type():

    with pytest.raises(AssertionError):
        a = sbmlcore.TrajectoryDistances(
            "tests/rpob-5uh6-3-warm.gro.gz",
            ["tests/rpob-5uh6-3-md-1-50ns-dt10ns-nojump.xtc"],
            "tests/3pl1.pdb",
            "resname RFP",
            "dist_RIF",
            distance_type="average",
        )


def test_not_list():

    with pytest.raises(AssertionError):
        a = sbmlcore.TrajectoryDistances(
            "tests/3pl1.pdb",
            "tests/3pl1.xtc",
            "tests/3pl1.pdb",
            "resname FE",
            "dist_FE",
        )


def test_start_time():

    with pytest.raises(AssertionError):
        a = sbmlcore.TrajectoryDistances(
            "tests/3pl1.pdb",
            ["tests/3pl1.xtc"],
            "tests/3pl1.pdb",
            "resname FE",
            "dist_FE",
            start_time=-1,
            end_time=12.2,
        )


def test_end_time():

    with pytest.raises(AssertionError):
        a = sbmlcore.TrajectoryDistances(
            "tests/3pl1.pdb",
            ["tests/3pl1.xtc"],
            "tests/3pl1.pdb",
            "resname FE",
            "dist_FE",
            start_time=5,
            end_time=2,
        )

    with pytest.raises(AssertionError):
        a = sbmlcore.TrajectoryDistances(
            "tests/3pl1.pdb",
            ["tests/3pl1.xtc"],
            "tests/3pl1.pdb",
            "resname FE",
            "dist_FE",
            start_time=20,
            end_time=10,
        )


def test_working():

    a = sbmlcore.TrajectoryDistances(
        "tests/rpob-5uh6-3-warm.gro.gz",
        ["tests/rpob-5uh6-3-md-1-50ns-dt10ns-nojump.xtc"],
        "tests/3pl1.pdb",
        "resname RFP",
        "dist_RIF",
    )

def test_exclude_percentiles():
    input = numpy.array([1,2,3,4,5,6,7,8,9,10,11,12,13,14,15])
    expected_output = numpy.array([3,4,5,6,7,8,9,10,11,12,13,14])
    output = sbmlcore.TrajectoryDistances._exclude_percentiles(input)
    assert output == expected_output


def test_init():
    a = sbmlcore.TrajectoryDistances(
        "tests/rpob-5uh6-3-warm.gro",
        [
            "tests/rpob-5uh6-3-md-1-50ns-dt100ps-nojump.xtc",
            "tests/rpob-5uh6-3-md-2-50ns-dt100ps-nojump.xtc",
        ],
        "tests/5uh6.pdb",
        "resname RFP",
        "max",
        percentile_exclusion=True,
    ).return_dist_df()
    b = pandas.read_csv("tests/5uh6_traj_distances.csv", index_col=0)
    pandas.testing.assert_frame_equal(a, b)


    
