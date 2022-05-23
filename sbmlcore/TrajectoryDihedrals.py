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

            #calculate time step for regeneration of trajectory during frame filtering
            dt = u.trajectory[1].time - u.trajectory[0].time

            if start_time is not None:
                u = TrajectoryDihedrals._filter_frames(
                    pdb_file, u, "start", start_time, dt
                )
                #rewriting the trajectory shifts time back to zero - therefore have to recalibrate end_time
                end_time = end_time - start_time
            if end_time is not None:
                u = TrajectoryDihedrals._filter_frames(
                    pdb_file, u, "end", end_time, dt
                )

            residues_all = u.select_atoms("name CA")

            if add_bonds:
                protein_res = TrajectoryDihedrals._add_bonds(u).select_atoms("protein")
            else:
                protein_res = u.select_atoms("protein")

            dihedrals = self.calculate_dihedrals(u, protein_res)

            if first_pass:
                dihedral_array = dihedrals
                first_pass = False
            else:
                dihedral_array = numpy.concatenate([dihedral_array, dihedrals])

        # inverts the array to make subseqeunt calculations more intuitive
        dihedral_array = dihedral_array.T

        if percentile_exclusion == True:
            dihedral_array = TrajectoryDihedrals._exclude_percentiles(dihedral_array)

        angles = self.apply_angle_type(dihedral_array)

        # pull segment ids from the static pdb file
        segids = [i for i in u_static.select_atoms("protein").residues.segids]

        # constructs the dictionary containihg the distances and assocaited residue labels
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

        nonetypes = self.search_nonetypes(trajectory)

        NaNs = numpy.empty((len(trajectory.trajectory), 1))
        NaNs[:] = 0

        if nonetypes[0] == 0:
            x = NaNs.copy()
        else:
            x = Dihedral(
                [eval(selection_call) for res in protein_res.residues[0 : nonetypes[0]]]
            ).run()
            x = numpy.concatenate((x.results["angles"], NaNs), axis=1)

        for i in range(0, len(nonetypes)):
            if i + 1 < len(nonetypes):
                y = Dihedral(
                    [
                        eval(selection_call)
                        for res in protein_res.residues[
                            nonetypes[i] + 1 : nonetypes[i + 1]
                        ]
                    ]
                ).run()
                x = numpy.concatenate((x, y.results["angles"], NaNs), axis=1)
            else:
                if nonetypes[-1] == len(protein_res.residues) - 1:
                    continue
                else:
                    y = Dihedral(
                        [
                            eval(selection_call)
                            for res in protein_res.residues[nonetypes[i] + 1 :]
                        ]
                    ).run()
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
        """calculates the specified angle for each residue"""

        if self.angle_type == "mean":
            angles = numpy.mean(dihedral_array, axis=1)
        elif self.angle_type == "min":
            angles = numpy.min(dihedral_array, axis=1)
        elif self.angle_type == "max":
            angles = numpy.max(dihedral_array, axis=1)
        elif self.angle_type == "median":
            angles = numpy.median(dihedral_array, axis=1)

        return angles

    def add_feature(self, existing_df):

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

        half_data = len(new_df[self.angle_name]) // 2

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
        """returns trajectory with frames greater and
        less than the specified start and end times"""

        coordinates = (
            MDAnalysis.analysis.base.AnalysisFromFunction(
                lambda ag: ag.positions.copy(), traj.atoms
            )
            .run()
            .results
        )

        if boundary == "start":
            bools = coordinates.times < spec_time
        elif boundary == "end":
            bools = coordinates.times >= spec_time
        filtered_times = numpy.delete(coordinates.times, bools)
        filtered_timeseries = numpy.delete(coordinates.timeseries, bools, axis=0)
        filtered_frames = numpy.delete(coordinates.frames, bools)

        coordinates = {
            "timeseries": filtered_timeseries,
            "frames": filtered_frames,
            "times": filtered_times,
        }

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
        bonds = MDAnalysis.topology.guessers.guess_bonds(
            protein_res, protein_res.positions
        )
        traj.add_TopologyAttr("bonds", bonds)
        return traj

    @staticmethod
    def _exclude_percentiles(data):
        """return array without 5% tails"""
        data_list = []
        for resnum in range(len(data)):
            if data[resnum].all() == numpy.zeros(len(data[resnum])).all():
                len_no_tails = len(data[resnum]) - (0.1 * len(data[resnum]))
                arr = numpy.zeros(int(len_no_tails)-1)
                data_list.append(arr)
            else:
                arr = data[resnum]
                p5 = numpy.percentile(arr, 5)
                p95 = numpy.percentile(arr, 95)
                data_list.append(arr[(arr >= p5) & (arr <= p95)])

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


# need to figure out way to exclude first frame + start/end times
