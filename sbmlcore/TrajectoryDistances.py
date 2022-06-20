import pathlib

import pandas
import numpy
import MDAnalysis

amino_acid_3to1letter = {
    "CYS": "C",
    "ASP": "D",
    "SER": "S",
    "GLN": "Q",
    "LYS": "K",
    "ILE": "I",
    "PRO": "P",
    "THR": "T",
    "PHE": "F",
    "ASN": "N",
    "GLY": "G",
    "HIS": "H",
    "LEU": "L",
    "ARG": "R",
    "TRP": "W",
    "ALA": "A",
    "VAL": "V",
    "GLU": "E",
    "TYR": "Y",
    "MET": "M",
}


class TrajectoryDistances(object):
    """
    Calculate and add distances from a molecular dynamics trajectory.

    Parameters
    ----------
    1st - path to structure file
    2nd - list containing paths to trajectory files
    3rd - path to pdb_file of static structure - used to pull segment ids as these are often lost in trajectory structure files
    4th - group of atoms you want to calculate the distances to - uses MDAnalysis syntax, and distances are calculated from the centre of mass of this whole selection to each Ca in the structure
    5th - your choice of name for the resulting distance column in the dataframe
    6th - type of distance metric that is being calculated (mean, median, max or min) (default=mean)
    7th - resid offsets for the different chains - must be a dictionary in the form {'segid': int, ...}.
    8th - desired starting time of the trajectory
    9th - desired end time of the trajectory
    10th - if percentile_exlcusion is set to True, only data between the 5th and 9th percentile is considered (default=False)

    E.g. a = sbmlcore.StructuralDistances('tests/5uh6.gro',['tests/5uh6.xtc'], '5hu6.pdb', 'resname RFP', 'RFP_distance', distance_type='median', offsets = {'A': 0, 'B': 0, 'C': -6}, 'start_time=1000', 'end_time=49000', percentile_exclusion=True)

    Returns
    -------
    Instantiating the class instantiates a df with segment ids, residue ids, residue names, and associated distances to the specified reference selection

    add_feature adds distances to existing mutation dataframe, and returns new joined dataframe

    """

    def __init__(
        self,
        pdb_file,
        trajectory_list,
        static_pdb,
        distance_selection,
        distance_name,
        distance_type="mean",
        offsets=None,
        start_time=None,
        end_time=None,
        percentile_exclusion=False,
    ):

        # check files exists
        assert pathlib.Path(pdb_file).is_file(), "File does not exist!"
        assert pathlib.Path(static_pdb).is_file(), "File does not exist!"

        # ensure trajectories are provided as a list
        assert isinstance(trajectory_list, list), "trajectories not provided as a list"
        for i in trajectory_list:
            assert pathlib.Path(i).is_file(), "File does not exist!"

        # ensure the distance type is legal, and the distance selection and name are strings
        assert distance_type in [
            "mean",
            "min",
            "max",
            "median",
        ], "provided distance_type not recognised!"
        assert isinstance(
            distance_selection, str
        ), "Distance selection must be a string!"
        assert isinstance(distance_name, str), "Distance name must be a string!"
        # checks whether percentile_exlcusion is True/False
        assert isinstance(
            percentile_exclusion, bool
        ), "Confidence interval must be True or False"

        if offsets is not None:
            # ensure the offsets are specified in a dictionary
            assert isinstance(
                offsets, dict
            ), "Offsets should be specified as a dictionary e.g. offsets = {'A': 3, 'B': -4}"

        # ensure start and end times are legal floats
        if start_time is not None:
            assert start_time > 0
            assert isinstance(start_time, float)
        if end_time is not None:
            assert end_time > 0
            assert isinstance(end_time, float)
            if start_time is not None:
                assert start_time < end_time

        first_pass = True

        for trajectory in trajectory_list:

            u = MDAnalysis.Universe(pdb_file, trajectory)
            u_static = MDAnalysis.Universe(static_pdb)

            reference_com = u.select_atoms(distance_selection).center_of_mass()

            # check atom selection exists
            assert (
                u.select_atoms(distance_selection).n_atoms > 0
            ), "Atom selection does not exist! Is your selection using the correct MDAnalysis syntax?"

            Ca_all = u.select_atoms("name CA")

            for ts in u.trajectory:

                if ts.frame == 0:
                    continue

                if start_time is not None and ts.time < start_time:
                    continue

                if end_time is not None and ts.time > end_time:
                    continue

                distances = MDAnalysis.lib.distances.distance_array(
                    reference_com, Ca_all.positions
                )

                if first_pass:
                    distance_array = distances
                    first_pass = False
                else:
                    distance_array = numpy.concatenate([distance_array, distances])

        # inverts the array to make subseqeunt calculations more intuitive
        distance_array = distance_array.T

        if percentile_exclusion == True:
            distance_array = TrajectoryDistances._exclude_percentiles(distance_array)

        # calculates distances as per the specified distance type
        if distance_type == "mean":
            distances = numpy.mean(distance_array, axis=1)
        elif distance_type == "min":
            distances = numpy.min(distance_array, axis=1)
        elif distance_type == "max":
            distances = numpy.max(distance_array, axis=1)
        elif distance_type == "median":
            distances = numpy.median(distance_array, axis=1)

        # Pulls segment ids from the static pdb file
        segids = [i for i in u_static.select_atoms("name CA").residues.segids]

        # constructs the dictionary containihg the distances and assocaited residue labels
        Ca_data = {
            "segid": segids,
            "resid": Ca_all.resids,
            "resname": Ca_all.resnames,
            distance_name: distances,
        }

        def one_letter(row):
            return amino_acid_3to1letter[row.resname]

        results = pandas.DataFrame(Ca_data)
        results["amino_acid"] = results.apply(one_letter, axis=1)
        results.drop(columns=["resname"], inplace=True)
        # otherwise column names post merge in add_feature are messy

        if offsets is not None:
            results = TrajectoryDistances._apply_offsets(results, offsets)

        self.results = results
        self.distance_name = distance_name

    def _add_feature(self, existing_df):
        """
        Adds distances to existing mutation dataframe, and returns new joined dataframe.

        Parameters
        ---------
        1st - existing dataframe

        E.g. if a = sbmlcore.TrajectoryDistances(...),
        use new_df = a.add_feature(existing_df)

        Returns
        --------
        New, joined dataframe

        """

        assert isinstance(
            existing_df, pandas.DataFrame
        ), "You must be adding the extra feature to an existing dataframe!"

        assert (
            self.distance_name not in existing_df.columns
        ), "You've already added that feature!"

        assert (
            "mutation" in existing_df.columns
        ), "Passed dataframe must contain a column called mutation"

        assert (
            "segid" in existing_df.columns
        ), "Passed dataframe must contain a column called segid containing chain information e.g. A"

        # Identifies the original one letter resname and resid (consistent with offset) so that the new feature can be subsequently linked to these
        def split_mutation(row):
            return pandas.Series([row.mutation[0], int(row.mutation[1:-1])])

        existing_df[["amino_acid", "resid"]] = existing_df.apply(split_mutation, axis=1)

        # Create MultiIndex using segid, resid and amino_acid
        existing_df.set_index(["segid", "resid", "amino_acid"], inplace=True)
        self.results.set_index(["segid", "resid", "amino_acid"], inplace=True)

        new_df = pandas.merge(
            existing_df, self.results, how="left", left_index=True, right_index=True
        )

        half_data = len(new_df[self.distance_name]) // 2

        total_nans = new_df[self.distance_name].isna().sum()

        assert (
            total_nans < half_data
        ), "Too many NaNs! Have you defined your offsets correctly?"

        for i in new_df[self.distance_name]:
            assert isinstance(i, float), "Distances must be floats!"

        new_df.reset_index(inplace=True)
        self.results.reset_index(inplace=True)

        new_df.drop(
            columns=[
                "amino_acid",
            ],
            inplace=True,
        )

        return new_df

    def return_dist_df(self):
        """returns the df generated by init method (useful for tests)"""
        return self.results

    @staticmethod
    def _exclude_percentiles(data):
        """returns array without 5% tails"""
        data_list = []
        for resnum in range(len(data)):
            arr = data[resnum]
            p5 = numpy.percentile(arr, 5)
            p95 = numpy.percentile(arr, 95)
            data_list.append(arr[(arr >= p5) & (arr <= p95)])
        data_arr = numpy.array(data_list)
        return data_arr

    @staticmethod
    def _apply_offsets(df, offsets):
        """returns dataframe with numbering offsets applied to the resids"""
        for chain in offsets:
            assert chain in set(
                df["segid"]
            ), "Need to specifify a segid that exists in the static pdb!"
            assert isinstance(
                offsets[chain], int
            ), "Offsets for each segid must be an integer!"

            resid_list = []
            for i in df.index:
                if df["segid"][i] == chain:
                    resid_list.append(df["resid"][i] + offsets[chain])
                else:
                    resid_list.append(df["resid"][i])

            df["resid"] = resid_list
        return df
