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
            "tests/5uh6.pdb",
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


def test_not_string():

    with pytest.raises(AssertionError):
        a = sbmlcore.TrajectoryDistances(
            "tests/3pl1.pdb",
            ["tests/3pl1.xtc"],
            "tests/3pl1.pdb",
            5,
            "dist_FE",
        )

    with pytest.raises(AssertionError):
        b = sbmlcore.TrajectoryDistances(
            "tests/3pl1.pdb",
            ["tests/3pl1.xtc"],
            "tests/3pl1.pdb",
            "resname FE",
            ["dist_FE"],
        )


def test_boolean():

    with pytest.raises(AssertionError):
        a = sbmlcore.TrajectoryDistances(
            "tests/3pl1.pdb",
            ["tests/3pl1.xtc"],
            "tests/3pl1.pdb",
            "resname FE",
            "dist_FE",
            percentile_exclusion="true",
        )


def test_start_time():

    with pytest.raises(AssertionError):
        a = sbmlcore.TrajectoryDistances(
            "tests/3pl1.pdb",
            ["tests/3pl1.xtc"],
            "tests/3pl1.pdb",
            "resname FE",
            "dist_FE",
            start_time=1,
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
            start_time=5.0,
            end_time=7,
        )

    with pytest.raises(AssertionError):
        a = sbmlcore.TrajectoryDistances(
            "tests/3pl1.pdb",
            ["tests/3pl1.xtc"],
            "tests/3pl1.pdb",
            "resname FE",
            "dist_FE",
            start_time=20.0,
            end_time=10.0,
        )


def test_runs():
    a = sbmlcore.TrajectoryDistances(
        "tests/rpob-5uh6-3-warm.gro.gz",
        ["tests/rpob-5uh6-3-md-1-50ns-dt10ns-nojump.xtc"],
        "tests/5uh6.pdb",
        "resname RFP",
        "dist_RIF",
    )


def test_exclude_percentiles():
    input = numpy.array([[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]])
    expected_output = numpy.array([[2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14]])
    output = sbmlcore.TrajectoryDistances._exclude_percentiles(input)
    numpy.testing.assert_array_equal(expected_output, output, verbose=True)


def test_apply_offsets():
    a = {
        "segid": ["A", "A", "A", "A", "A", "B", "B", "C", "C", "C", "C", "C", "C"],
        "resid": [1, 2, 3, 4, 5, 1, 2, 12, 13, 14, 15, 16, 17],
    }
    b = {
        "segid": ["A", "A", "A", "A", "A", "B", "B", "C", "C", "C", "C", "C", "C"],
        "resid": [-2, -1, 0, 1, 2, 1, 2, 19, 20, 21, 22, 23, 24],
    }
    input_df = pandas.DataFrame.from_dict(a)
    expected_output = pandas.DataFrame.from_dict(b)
    input_offsets = {"A": -3, "B": 0, "C": 7}
    output = sbmlcore.TrajectoryDistances._apply_offsets(input_df, input_offsets)
    pandas.testing.assert_frame_equal(output, expected_output)


def test_init():
    a = sbmlcore.TrajectoryDistances(
        "tests/rpob-5uh6-3-warm.gro.gz",
        [
            "tests/rpob-5uh6-3-md-1-50ns-dt10ns-nojump.xtc",
            "tests/rpob-5uh6-3-md-2-50ns-dt10ns-nojump.xtc",
            "tests/rpob-5uh6-3-md-3-50ns-dt10ns-nojump.xtc",
        ],
        "tests/5uh6.pdb",
        "resname RFP",
        "max RFP",
        distance_type="max",
        offsets={"A": 0, "B": 0, "C": -6},
        percentile_exclusion=True,
    ).return_dist_df()
    b = pandas.read_csv("tests/5uh6_traj_distances.csv", index_col=0)
    pandas.testing.assert_frame_equal(a, b)


def test._add_feature():
    # although this requires running the entire class to merge the dfs,
    # all other components of the class should have been individually tested by now
    a = {
        "segid": ["A", "A", "A", "B", "C", "C"],
        "mutation": ["I3D", "S4K", "Q5V", "R6D", "S450F", "D435F"],
    }
    df = pandas.DataFrame.from_dict(a)
    b = sbmlcore.TrajectoryDistances(
        "tests/rpob-5uh6-3-warm.gro.gz",
        [
            "tests/rpob-5uh6-3-md-1-50ns-dt10ns-nojump.xtc",
            "tests/rpob-5uh6-3-md-2-50ns-dt10ns-nojump.xtc",
            "tests/rpob-5uh6-3-md-3-50ns-dt10ns-nojump.xtc",
        ],
        "tests/5uh6.pdb",
        "resname RFP",
        "max RFP",
        distance_type="max",
        offsets={"A": 0, "B": 0, "C": -6},
        percentile_exclusion=True,
    )
    c = sbmlcore.TrajectoryDistances(
        "tests/rpob-5uh6-3-warm.gro.gz",
        [
            "tests/rpob-5uh6-3-md-1-50ns-dt10ns-nojump.xtc",
            "tests/rpob-5uh6-3-md-2-50ns-dt10ns-nojump.xtc",
            "tests/rpob-5uh6-3-md-3-50ns-dt10ns-nojump.xtc",
        ],
        "tests/5uh6.pdb",
        "resname RFP",
        "mean RFP",
        distance_type="mean",
        offsets={"A": 0, "B": 0, "C": -6},
        percentile_exclusion=True,
    )
    test_df = pandas.read_csv("tests/5uh6_added_traj_distances.csv", index_col=0)
    df = b._add_feature(df)
    df = c._add_feature(df)
    pandas.testing.assert_frame_equal(test_df, df)
