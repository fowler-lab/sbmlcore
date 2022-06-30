import pandas
import numpy
import pytest
import sbmlcore


def test_missing_file():

    with pytest.raises(AssertionError):
        a = sbmlcore.TrajectoryDistances(
            "tests/dhfr-3fre-tmp-1-1.gro",
            "tests/missing.xtc",
            "tests/dhfr-3fre-tmp-1-1.pdb",
            "resname NDP",
            "dist_NDP",
        )


def test_wrong_type():

    with pytest.raises(AssertionError):
        a = sbmlcore.TrajectoryDistances(
            "tests/dhfr-3fre-tmp-1-1.gro",
            ["tests/dhfr-3fre-tmp-1-2-nojump-skip100.xtc"],
            "tests/dhfr-3fre-tmp-1-1.pdb",
            "resname RFP",
            "dist_RIF",
            distance_type="average",
        )


def test_not_list():

    with pytest.raises(AssertionError):
        a = sbmlcore.TrajectoryDistances(
            "tests/dhfr-3fre-tmp-1-1.gro",
            "tests/dhfr-3fre-tmp-1-2-nojump-skip100.xtc",
            "tests/dhfr-3fre-tmp-1-1.pdb",
            "resname NDP",
            "dist_NDP",
        )


def test_not_string():

    with pytest.raises(AssertionError):
        a = sbmlcore.TrajectoryDistances(
            "tests/dhfr-3fre-tmp-1-1.gro",
            ["tests/dhfr-3fre-tmp-1-2-nojump-skip100.xtc"],
            "tests/dhfr-3fre-tmp-1-1.pdb",
            5,
            "dist_NDP",
        )

    with pytest.raises(AssertionError):
        b = sbmlcore.TrajectoryDistances(
            "tests/dhfr-3fre-tmp-1-1.gro",
            ["tests/dhfr-3fre-tmp-1-2-nojump-skip100.xtc"],
            "tests/dhfr-3fre-tmp-1-1.pdb",
            "resname NDP",
            ["dist_NDP"],
        )


def test_boolean():

    with pytest.raises(AssertionError):
        a = sbmlcore.TrajectoryDistances(
            "tests/dhfr-3fre-tmp-1-1.gro",
            ["tests/dhfr-3fre-tmp-1-2-nojump-skip100.xtc"],
            "tests/dhfr-3fre-tmp-1-1.pdb",
            "resname NDP",
            "dist_NDP",
            percentile_exclusion="true",
        )


def test_start_time():

    with pytest.raises(AssertionError):
        a = sbmlcore.TrajectoryDistances(
            "tests/dhfr-3fre-tmp-1-1.gro",
            ["tests/dhfr-3fre-tmp-1-2-nojump-skip100.xtc"],
            "tests/dhfr-3fre-tmp-1-1.pdb",
            "resname NDP",
            "dist_NDP",
            start_time=1,
            end_time=12.2,
        )


def test_end_time():

    with pytest.raises(AssertionError):
        a = sbmlcore.TrajectoryDistances(
            "tests/dhfr-3fre-tmp-1-1.gro",
            ["tests/dhfr-3fre-tmp-1-2-nojump-skip100.xtc"],
            "tests/dhfr-3fre-tmp-1-1.pdb",
            "resname NDP",
            "dist_NDP",
            start_time=5.0,
            end_time=7,
        )

    with pytest.raises(AssertionError):
        a = sbmlcore.TrajectoryDistances(
            "tests/dhfr-3fre-tmp-1-1.gro",
            ["tests/dhfr-3fre-tmp-1-2-nojump-skip100.xtc"],
            "tests/dhfr-3fre-tmp-1-1.pdb",
            "resname NDP",
            "dist_NDP",
            start_time=20.0,
            end_time=10.0,
        )


def test_runs():
    a = sbmlcore.TrajectoryDistances(
        "tests/dhfr-3fre-tmp-1-1.pdb",
        ["tests/dhfr-3fre-tmp-1-2-nojump-skip100.xtc"],
        "tests/dhfr-3fre-tmp-1-1.pdb",
        "resname TMP",
        "dist_TMP",
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
        "tests/dhfr-3fre-tmp-1-1.gro",
        [
            "tests/dhfr-3fre-tmp-1-2-nojump-skip100.xtc",
            "tests/dhfr-3fre-tmp-1-2-nojump-skip100.xtc",
            "tests/dhfr-3fre-tmp-1-2-nojump-skip100.xtc",
        ],
        "tests/dhfr-3fre-tmp-1-1.pdb",
        "resname TMP",
        "max_TMP",
        distance_type="max",
        offsets={"A": 0},
        percentile_exclusion=True,
    ).return_dist_df()

    assert a.max_TMP.sum() == pytest.approx(5120.556)


def test_add_feature():
    # although this requires running the entire class to merge the dfs,
    # all other components of the class should have been individually tested by now
    a = {
        "segid": ["A", "A", "A", "A"],
        "mutation": ["T1D", "L2K", "S3V", "S3F"],
    }
    df = pandas.DataFrame.from_dict(a)
    features = sbmlcore.FeatureDataset(
        df, species="S. aureus", gene="folA", protein="DHFR"
    )
    b = sbmlcore.TrajectoryDistances(
        "tests/dhfr-3fre-tmp-1-1.gro",
        [
            "tests/dhfr-3fre-tmp-1-2-nojump-skip100.xtc",
        ],
        "tests/dhfr-3fre-tmp-1-1.pdb",
        "resname NDP",
        "max NDP",
        distance_type="max",
        offsets={"A": 0},
        percentile_exclusion=True,
    )
    c = sbmlcore.TrajectoryDistances(
        "tests/dhfr-3fre-tmp-1-1.gro",
        [
            "tests/dhfr-3fre-tmp-1-2-nojump-skip100.xtc",
        ],
        "tests/dhfr-3fre-tmp-1-1.pdb",
        "resname TMP",
        "min TMP",
        distance_type="min",
        offsets={"A": 0},
        percentile_exclusion=True,
    )
    d = sbmlcore.TrajectoryDistances(
        "tests/dhfr-3fre-tmp-1-1.gro",
        [
            "tests/dhfr-3fre-tmp-1-2-nojump-skip100.xtc",
        ],
        "tests/dhfr-3fre-tmp-1-1.pdb",
        "resname NDP",
        "mean NDP",
        distance_type="mean",
        offsets={"A": 0},
        percentile_exclusion=True,
    )

    # method 1:
    features.add_feature([b])

    # method 2:
    features_df = (features + c).df

    # method 3:
    features_df = d._add_feature(features_df)

    test_df = pandas.read_csv('tests/3fre_added_traj_distances.csv', index_col=0)
    pandas.testing.assert_frame_equal(test_df, features_df)
