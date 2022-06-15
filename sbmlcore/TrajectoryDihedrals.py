import pathlib

import pandas
import numpy
import MDAnalysis
from MDAnalysis.analysis.dihedrals import Dihedral

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


class TrajectoryDihedrals(object):
    """Calculate and add dihedral angles from a molecular dynamics trajectory.

    Parameters
    ----------
    1st - path to structure file
    2nd - list containing paths to trajectory files
    3rd - path to pdb_file of static structure - used to pull segment ids as these are often lost in trajectory structure files
    4th - the target angle (phi, psi, omega)
    5th - your choice of name for the resulting distance column in the dataframe
    6th - type of angle metic that is being calculated (mean, median, max, or min)(default=mean)
    7th - if the structure file doesn't contain all bonds, one can set add_bonds=True to fire an MDAnalysis bond prediction algorithm
    8th - resid offsets for the different chains - must be a dictionary in the form {'segid':int, ...}.
    9th - desired starting time of the trajectory
    10th - deisred end time of the trajectory
    11th - if percentile_exclusion is set to True, only data between the 5th and 9th percentile is considered (default=False)

    E.g. a = a = sbmlcore.TrajectoryDihedrals(
        "./tests/rpob-5uh6-3-warm.gro.gz",
        [
            "./tests/rpob-5uh6-3-md-1-50ns-dt10ns-nojump.xtc",
        ],
        "./tests/5uh6.pdb",
        "phi",
        "mean_phi",
        angle_type="mean",
        add_bonds=True,
        offsets = {'A': 0, 'B': 0, 'C': -6},
        percentile_exclusion=True
        )

    Returns
    ----------
    Instantiating the class instantiates a df with segment ids, residue ids, residue names, and associated distances to the specified reference selection

    add_feature adds distances to existing mutation dataframe, and returns new joined dataframed

    """

    def __init__(
        self,
        pdb_file,
        trajectory_list,
        static_pdb,
        dihedral,
        angle_name,
        angle_type="mean",
        add_bonds=False,
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

        # ensure the dihedral is legal
        assert dihedral in ["phi", "psi", "omega"]

        # ensure the angle type is legal, and the distance selection name are strings
        assert angle_type in [
            "mean",
            "min",
            "max",
            "median",
        ], "provided angle_type not recognised!"
        assert isinstance(angle_name, str), "Angle name must be a string!"

        # checks whether add_bonds is True/False
        assert isinstance(add_bonds, bool), "add_bonds must be True or False"

        # checks whether percentile_exclusion is True/False
        assert isinstance(
            percentile_exclusion, bool
        ), "percentile_exclusion must be True or False"

        if offsets is not None:
            # ensure the offsets are specificed in a dictionary
            assert isinstance(
                offsets, dict
            ), "Offsets should be specific as a dictionary eg. offsets = {'A':3, 'B':-4}"

        # ensure start and end times are legal floats
        if start_time is not None:
            assert start_time > 0
            assert isinstance(start_time, float), "start_time must be a float"
        if end_time is not None:
            assert end_time > 0
            assert isinstance(end_time, float), "end_time must be a float"
            if start_time is not None:
                assert start_time < end_time, "start_ime should be less than end_time"

        self.dihedral = dihedral
        self.angle_type = angle_type

        first_pass = True

        for trajectory in trajectory_list:

            u = MDAnalysis.Universe(pdb_file, trajectory)
            u_static = MDAnalysis.Universe(static_pdb)

            # calculate time step for regeneration of trajectory during frame filtering
            dt = u.trajectory[1].time - u.trajectory[0].time

            # remove the first frame (often the pdb structure)
            u = TrajectoryDihedrals._filter_frames(
                pdb_file, u, "start", u.trajectory[1].time, dt
            )

            # remove all frames that exist before the specified start time
            if start_time is not None:
                u = TrajectoryDihedrals._filter_frames(
                    pdb_file, u, "start", start_time, dt
                )
                # rewriting the trajectory shifts time back to zero - therefore have to recalibrate end_time
                end_time = end_time - start_time
            # remove all frames that exist after the specified start time
            if end_time is not None:
                u = TrajectoryDihedrals._filter_frames(pdb_file, u, "end", end_time, dt)

            residues_all = u.select_atoms("name CA")

            if add_bonds:
                # Runs the MDAnalysis bond prediction algorithm
                protein_res = TrajectoryDihedrals._add_bonds(u).select_atoms("protein")
            else:
                protein_res = u.select_atoms("protein")

            # calculate all angles of the specified bond - returns an array of shape (timesteps, residues)
            dihedrals = self.calculate_dihedrals(u, protein_res)

            if first_pass:
                dihedral_array = dihedrals
                first_pass = False
            else:
                dihedral_array = numpy.concatenate([dihedral_array, dihedrals])

        # inverts the array to make subseqeunt calculations more intuitive
        dihedral_array = dihedral_array.T

        if percentile_exclusion == True:
            # exclude 5% tails
            dihedral_array = TrajectoryDihedrals._exclude_percentiles(dihedral_array)

        # calculate the specified metric for each residue (mean, max, min, median)
        angles = self.apply_angle_type(dihedral_array)

        # pull segment ids from the static pdb file
        segids = [i for i in u_static.select_atoms("protein").residues.segids]

        # constructs the dictionary containing the distances and assocaited residue labels
        data = {
            "segid": segids,
            "resid": residues_all.resids,
            "resname": residues_all.resnames,
            angle_name: angles,
        }

        def one_letter(row):
            return amino_acid_3to1letter[row.resname]

        results = pandas.DataFrame(data)
        results["amino_acid"] = results.apply(one_letter, axis=1)
        results.drop(columns=["resname"], inplace=True)
        # otherwise column names post merge in add_feature are messy

        if offsets is not None:
            results = TrajectoryDihedrals._apply_offsets(results, offsets)

        self.results = results
        self.angle_name = angle_name

    def calculate_dihedrals(self, trajectory, protein_res):
        """
        Calculates specified dihedral angles for all residues in the protein
        and returns array of shape (timesteps, residues)
        """

        selection_call = "res." + self.dihedral + "_selection()"

        # generate list of nonetype dihedral indexes (residue index)
        # (MDAnalysis can't calculate angles for a small number of specific residues depending on the angle)
        nonetypes = self.search_nonetypes(trajectory)

        # Create an array of zeros to use in place of nonetypes
        NaNs = numpy.empty((len(trajectory.trajectory), 1))
        NaNs[:] = 0

        # if the first residue has nonetype angles
        if nonetypes[0] == 0:
            x = NaNs.copy()
        else:
            # calculate dihedrals for residues up until the first nonetype residue
            x = Dihedral(
                [eval(selection_call) for res in protein_res.residues[0 : nonetypes[0]]]
            ).run()
            # add a zero array to x to represent the first nonetype residue
            x = numpy.concatenate((x.results["angles"], NaNs), axis=1)

        for i in range(0, len(nonetypes)):
            if i + 1 < len(nonetypes):
                # calculate angles from residue after the current nonetype to the next nonetype = a 'segment of angles'
                y = Dihedral(
                    [
                        eval(selection_call)
                        for res in protein_res.residues[
                            nonetypes[i] + 1 : nonetypes[i + 1]
                        ]
                    ]
                ).run()
                # add the segment of angles to the dihedral array
                x = numpy.concatenate((x, y.results["angles"], NaNs), axis=1)
            else:
                # if the final residue is a nonetype
                if nonetypes[-1] == len(protein_res.residues) - 1:
                    continue
                else:
                    # if not, calculate angles right through to the last residue
                    y = Dihedral(
                        [
                            eval(selection_call)
                            for res in protein_res.residues[nonetypes[i] + 1 :]
                        ]
                    ).run()
                    # add the segment of angles to the dihedral array
                    x = numpy.concatenate((x, y.results["angles"]), axis=1)

        return x

    def search_nonetypes(self, traj):
        """searches for nonetype dihedral angles and returns a list of their indexes"""

        selection_call = "i." + self.dihedral + "_selection()"

        index = 0
        nonetype_list = []
        residues = traj.select_atoms("protein")
        for i in residues.residues:
            if not eval(selection_call):
                nonetype_list.append(index)
            index += 1
        return nonetype_list

    def apply_angle_type(self, dihedral_array):
        """calculates the specified angle for each residue across all frames"""

        if self.angle_type == "mean":
            angles = numpy.mean(dihedral_array, axis=1)
        elif self.angle_type == "min":
            angles = numpy.min(dihedral_array, axis=1)
        elif self.angle_type == "max":
            angles = numpy.max(dihedral_array, axis=1)
        elif self.angle_type == "median":
            angles = numpy.median(dihedral_array, axis=1)

        return angles

    def _add_feature(self, existing_df):
        """
        Adds angles to existing mutation dataframe, and returns new joined dataframe.
        Arguments: existing dataframe
        e.g. if a = sbmlcore.TrajectoryDihedrals(...),
        use new_df = a.add_feature(existing_df)
        """

        assert isinstance(
            existing_df, pandas.DataFrame
        ), "You must be adding the extra feature to an existing dataframe!"

        assert (
            self.angle_name not in existing_df.columns
        ), "You've already added that feature!"

        assert (
            "mutation" in existing_df.columns
        ), "Passed dataframe must contain a column called mutation"

        assert (
            "segid" in existing_df.columns
        ), "Passed dataframe must contain a column called segid containing chain information e.g. A "

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

        # To check that the chain offsets have been correctly set
        # and/or that structural data is of adequate quality -
        # defining this as the no NaNs should be less than half no rows

        half_data = (
            len(new_df[self.angle_name]) // 2
        )  # // divides and rounds DOWN to nearest int

        total_nans = new_df[self.angle_name].isna().sum()

        assert (
            total_nans < half_data
        ), "Too many NaNs! Have you defined your offsets correctly?"

        for i in new_df[self.angle_name]:
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

    def return_angle_df(self):
        """returns the df generated by init method (useful for tests)"""
        return self.results

    @staticmethod
    def _filter_frames(pdb, traj, boundary, spec_time, dt):
        """returns universe with frames greater and
        less than the specified start and end times"""

        # Becuase a new universe is essentially being created, every coordinate in the original is needed
        coordinates = (
            MDAnalysis.analysis.base.AnalysisFromFunction(
                lambda ag: ag.positions.copy(), traj.atoms
            )
            .run()
            .results
        )

        # Create arrays to signifiy whether time condition has passed or failed
        if boundary == "start":
            bools = coordinates.times < spec_time
        elif boundary == "end":
            bools = coordinates.times >= spec_time
        # delete times, ts, and frames for which the above conditions have passed (=outside of desired range)
        filtered_times = numpy.delete(coordinates.times, bools)
        filtered_timeseries = numpy.delete(coordinates.timeseries, bools, axis=0)
        filtered_frames = numpy.delete(coordinates.frames, bools)

        coordinates = {
            "timeseries": filtered_timeseries,
            "frames": filtered_frames,
            "times": filtered_times,
        }

        # recreate universe with new specifications (reads into memory -
        # //if running this class bombs the machine, its likely due to insufficent memory )
        u = MDAnalysis.Universe(
            pdb,
            coordinates["timeseries"],
            dimensions=traj.dimensions,
            dt=dt,
            format=MDAnalysis.coordinates.memory.MemoryReader,
        )
        return u

    @staticmethod
    def _add_bonds(traj):
        """add bonds to protein (only) in a trajectory"""

        protein_res = traj.select_atoms("protein")
        # run the bond guessing algorithm
        bonds = MDAnalysis.topology.guessers.guess_bonds(
            protein_res, protein_res.positions
        )
        # add the bonds to the trajectory
        traj.add_TopologyAttr("bonds", bonds)
        return traj

    @staticmethod
    def _exclude_percentiles(data):
        """return array without 5% tails"""
        data_list = []
        for resnum in range(len(data)):
            # calculate the expected length of an array after the tailes have been removed
            len_no_tails = len(data[resnum]) - (0.1 * len(data[resnum]))
            if data[resnum].all() == numpy.zeros(len(data[resnum])).all():
                # if the arrays is just an array of zeros, numpy can't calculate p5 nor p95 -
                # //so we do it manually
                arr = numpy.zeros(int(len_no_tails) - 1)
                data_list.append(arr)
            else:
                arr = data[resnum]
                p5 = numpy.percentile(arr, 5)
                p95 = numpy.percentile(arr, 95)
                data_list.append(arr[(arr > p5) & (arr <= p95)])

        data_arr = numpy.array(data_list)

        return data_arr

    @staticmethod
    def _apply_offsets(df, offsets):
        """returns dataframe with number offsets applied to the resids"""
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
