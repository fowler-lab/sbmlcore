import pandas
import numpy
import pytest
import sbmlcore
import MDAnalysis


def test_missing_file():

    with pytest.raises(AssertionError):
        a = sbmlcore.TrajectoryDihedrals(
            "tests/dhfr-3fre-tmp-1-1.gro",
            "tests/missing.xtc",
            "tests/dhfr-3fre-tmp-1-1.pdb",
            "psi",
            "mean_psi",
        )

    with pytest.raises(AssertionError):
        a = sbmlcore.TrajectoryDihedrals(
            "tests/dhfr-3fre-tmp-1-1.gro",
            ["tests/dhfr-3fre-tmp-1-2-nojump-skip100.xtc"],
            "tests/missing.pdb",
            "psi",
            "mean_psi",
        )

    with pytest.raises(AssertionError):
        a = sbmlcore.TrajectoryDihedrals(
            "tests/missing.gro",
            ["tests/dhfr-3fre-tmp-1-2-nojump-skip100.xtc"],
            "tests/dhfr-3fre-tmp-1-1.pdb",
            "psi",
            "mean_psi",
        )


def test_wrong_type():

    with pytest.raises(AssertionError):
        a = sbmlcore.TrajectoryDihedrals(
            "tests/dhfr-3fre-tmp-1-1.gro",
            ["tests/dhfr-3fre-tmp-1-2-nojump-skip100.xtc"],
            "tests/dhfr-3fre-tmp-1-1.pdb",
            "psi",
            "mean_psi",
            angle_type="average",
        )


def test_illegal_dihedral():

    with pytest.raises(AssertionError):
        a = sbmlcore.TrajectoryDihedrals(
            "tests/dhfr-3fre-tmp-1-1.gro",
            ["tests/dhfr-3fre-tmp-1-2-nojump-skip100.xtc"],
            "tests/dhfr-3fre-tmp-1-1.pdb",
            "Ramachandron",
            "mean_Ramachandron",
        )


def test_boolean():

    with pytest.raises(AssertionError):
        a = sbmlcore.TrajectoryDihedrals(
            "tests/dhfr-3fre-tmp-1-1.gro",
            ["tests/dhfr-3fre-tmp-1-2-nojump-skip100.xtc"],
            "tests/dhfr-3fre-tmp-1-1.pdb",
            "psi",
            "mean_psi",
            percentile_exclusion="Yes",
        )


def test_not_list():
    with pytest.raises(AssertionError):
        a = sbmlcore.TrajectoryDihedrals(
            "tests/dhfr-3fre-tmp-1-1.gro",
            "tests/dhfr-3fre-tmp-1-2-nojump-skip100.xtc",
            "tests/dhfr-3fre-tmp-1-1.pdb",
            "psi",
            "mean_psi",
        )


def test_start_time():
    with pytest.raises(AssertionError):
        a = sbmlcore.TrajectoryDihedrals(
            "tests/dhfr-3fre-tmp-1-1.gro",
            ["tests/dhfr-3fre-tmp-1-2-nojump-skip100.xtc"],
            "tests/dhfr-3fre-tmp-1-1.pdb",
            "psi",
            "mean_psi",
            start_time=5,
        )


def test_end_time():
    with pytest.raises(AssertionError):
        a = sbmlcore.TrajectoryDihedrals(
            "tests/dhfr-3fre-tmp-1-1.gro",
            ["tests/dhfr-3fre-tmp-1-2-nojump-skip100.xtc"],
            "tests/dhfr-3fre-tmp-1-1.pdb",
            "psi",
            "mean_psi",
            start_time=12.0,
            end_time=2.0,
        )


def test_runs():
    a = sbmlcore.TrajectoryDihedrals(
        "tests/dhfr-3fre-tmp-1-1.gro",
        ["tests/dhfr-3fre-tmp-1-2-nojump-skip100.xtc"],
        "tests/dhfr-3fre-tmp-1-1.pdb",
        "phi",
        "mean_phi",
    )


def test_exclude_percentiles():
    input = numpy.array([[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]])
    expected_output = numpy.array([[2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14]])
    output = sbmlcore.TrajectoryDihedrals._exclude_percentiles(input)
    numpy.testing.assert_array_equal(expected_output, output, verbose=True)

    input = numpy.array([[0, 0, 0, 0, 0, 0, 0, 0, 0, 0]])
    expected_output = numpy.array([[0, 0, 0, 0, 0, 0, 0, 0]])
    output = sbmlcore.TrajectoryDihedrals._exclude_percentiles(input)
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


def test_filter_frames():

    u = MDAnalysis.Universe(
        "tests/dhfr-3fre-tmp-1-1.gro", "tests/dhfr-3fre-tmp-1-2-nojump-skip100.xtc"
    )
    dt = u.trajectory[1].time - u.trajectory[0].time
    ut = sbmlcore.TrajectoryDihedrals._filter_frames(
        "tests/dhfr-3fre-tmp-1-1.gro", u, "start", 2000, dt
    )
    test_coordinates = (
        MDAnalysis.analysis.base.AnalysisFromFunction(
            lambda ag: ag.positions.copy(), ut.atoms
        )
        .run()
        .results
    )
    test_coordinates = pandas.DataFrame.from_dict(test_coordinates.timeseries[0])
    expected_coordinates = pandas.read_csv(
        "tests/test_start_coordinates_frame0.csv", index_col=0
    )
    expected_coordinates.rename(columns={"0": 0, "1": 1, "2": 2}, inplace=True)
    pandas.testing.assert_frame_equal(
        test_coordinates, expected_coordinates.astype("float32")
    )

    ut = sbmlcore.TrajectoryDihedrals._filter_frames(
        "tests/dhfr-3fre-tmp-1-1.gro", u, "end", 8000, dt
    )
    test_coordinates = (
        MDAnalysis.analysis.base.AnalysisFromFunction(
            lambda ag: ag.positions.copy(), ut.atoms
        )
        .run()
        .results
    )
    assert len(test_coordinates.timeseries) == 8


def test_apply_angle_type():

    instantiated_class = sbmlcore.TrajectoryDihedrals(
        "tests/dhfr-3fre-tmp-1-1.gro",
        ["tests/dhfr-3fre-tmp-1-2-nojump-skip100.xtc"],
        "tests/dhfr-3fre-tmp-1-1.pdb",
        "phi",
        "phi",
        angle_type="mean",
    )
    input = numpy.array([[2, 4, 5, 6, 7, 8, 9], [6, 7, 8, 9, 0, 11, 12]])
    output = instantiated_class.apply_angle_type(input)
    expected_output = numpy.array([41 / 7, 53 / 7])
    numpy.testing.assert_array_equal(expected_output, output, verbose=True)


def test_search_nonetypes():

    u = MDAnalysis.Universe(
        "tests/dhfr-3fre-tmp-1-1.gro", "tests/dhfr-3fre-tmp-1-2-nojump-skip100.xtc"
    )
    instantiated_class = sbmlcore.TrajectoryDihedrals(
        "tests/dhfr-3fre-tmp-1-1.gro",
        ["tests/dhfr-3fre-tmp-1-2-nojump-skip100.xtc"],
        "tests/dhfr-3fre-tmp-1-1.pdb",
        "psi",
        "psi",
    )
    assert instantiated_class.search_nonetypes(u) == [156]


def test_calculate_dihedrals():

    u = MDAnalysis.Universe(
        "tests/dhfr-3fre-tmp-1-1.gro", "tests/dhfr-3fre-tmp-1-2-nojump-skip100.xtc"
    )
    instantiated_class = sbmlcore.TrajectoryDihedrals(
        "tests/dhfr-3fre-tmp-1-1.gro",
        ["tests/dhfr-3fre-tmp-1-2-nojump-skip100.xtc"],
        "tests/dhfr-3fre-tmp-1-1.pdb",
        "phi",
        "phi",
    )
    protein_res = u.select_atoms("protein")
    output = instantiated_class.calculate_dihedrals(u, protein_res)
    output = pandas.DataFrame(output)
    expected_output = pandas.read_csv("tests/test_dihedrals.csv", index_col=0)
    expected_output.columns = expected_output.columns.astype(int)

    pandas.testing.assert_frame_equal(output, expected_output.astype("float64"))


def test_init():
    a = sbmlcore.TrajectoryDihedrals(
        "tests/dhfr-3fre-tmp-1-1.gro",
        [
            "tests/dhfr-3fre-tmp-1-2-nojump-skip100.xtc",
        ],
        "tests/dhfr-3fre-tmp-1-1.pdb",
        "phi",
        "phi",
        angle_type="max",
        offsets={"A": 0},
        percentile_exclusion=True,
        start_time=2000.0,
        end_time=8000.0,
    ).return_angle_df()
    b = pandas.read_csv("tests/3fre_traj_angles.csv", index_col=0)
    pandas.testing.assert_frame_equal(a, b)


def test_add_feature():

    a = {
        "segid": ["A", "A", "A", "A"],
        "mutation": ["T1D", "L2K", "S3V", "S3F"],
    }
    df = pandas.DataFrame.from_dict(a)
    features = sbmlcore.FeatureDataset(
        df, species="S. aureus", gene="folA", protein="DHFR"
    )
    b = sbmlcore.TrajectoryDihedrals(
        "tests/dhfr-3fre-tmp-1-1.gro",
        [
            "tests/dhfr-3fre-tmp-1-2-nojump-skip100.xtc",
        ],
        "tests/dhfr-3fre-tmp-1-1.pdb",
        "psi",
        "max psi",
        angle_type="max",
        offsets={"A": 0},
        percentile_exclusion=True,
        end_time=8000.0,
    )
    c = sbmlcore.TrajectoryDihedrals(
        "tests/dhfr-3fre-tmp-1-1.gro",
        [
            "tests/dhfr-3fre-tmp-1-2-nojump-skip100.xtc",
        ],
        "tests/dhfr-3fre-tmp-1-1.pdb",
        "omega",
        "min omega",
        angle_type="min",
        offsets={"A": 0},
        percentile_exclusion=True,
        end_time=8000.0,
    )

    d = sbmlcore.TrajectoryDihedrals(
        "tests/dhfr-3fre-tmp-1-1.gro",
        [
            "tests/dhfr-3fre-tmp-1-2-nojump-skip100.xtc",
        ],
        "tests/dhfr-3fre-tmp-1-1.pdb",
        "phi",
        "mean phi",
        angle_type="mean",
        offsets={"A": 0},
        percentile_exclusion=True,
        end_time=8000.0,
    )

    # method 1:
    features.add_feature([b])

    # method 2:
    features_df = (features + c).df

    # method 3:
    features_df = d._add_feature(features_df)

    test_df = pandas.read_csv("tests/3fre_added_traj_angles.csv", index_col=0)
    pandas.testing.assert_frame_equal(test_df, features_df)
