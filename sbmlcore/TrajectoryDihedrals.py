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
        ), "Confidence interval must be True or False"

        if offsets is not None:
            # ensure the offsets are specificed in a dictionary
            assert isinstance(
                offsets, dict
            ), "Offsets should be specific as a dictionary eg. offsets = {'A':3, 'B':-4}"

        # ensure start and end times are legal floats
        if start_time is not None:
            assert start_time > 0
            assert isinstance(start_time, float)
        if end_time is not None:
            assert end_time > 0
            assert isinstance(end_time, float)
            if start_time is not None:
                assert start_time < end_time

        self.dihedral = dihedral

        for trajectory in trajectory_list:

            u = MDAnalysis.Universe(pdb_file, trajectory)
            u_static = MDAnalysis.Universe(static_pdb)

            if add_bonds:
                protein_res = add_bonds(trajectory).select_atoms("protein")
            else:
                protein_res = trajectory.select_atoms("protein")

            self.calculate_dihedrals(trajectory, protein_res)
            #STOPPED HERE            

            

    def calculate_dihedrals(self, trajectory, protein_res):
        """
        Calculates specified dihedral angles for all residues in the protein 
        and returns array of shape (timesteps, residues)
        """

        selection_call = 'res.' + self.dihedral + '_selection()'

        nonetypes = TrajectoryDihedrals._search_nonetypes(trajectory)

        NaNs = numpy.empty((len(trajectory.trajectory), 1))
        NaNs[:] = 0

        if nonetypes[0] == 0:
            x = NaNs.copy()
        else:
            x = Dihedral([eval(selection_call) for res in protein_res.residues[0:nonetypes[0]]]).run()
            x = numpy.concatenate((x.results['angles'], NaNs), axis=1)
        
        for i in range(0, len(nonetypes)):
            if i+1 < len(nonetypes):
                y = Dihedral([eval(selection_call) for res in protein_res.residues[nonetypes[i]+1:nonetypes[i+1]]]).run()
                x = numpy.concatenate((x, y.results['angles'], NaNs), axis=1)
            else:
                if nonetypes[-1] == len(protein_res.residues)-1:
                    continue
                else:
                    y = Dihedral([eval(selection_call) for res in protein_res.residues[nonetypes[i]+1:]]).run()
                    x = numpy.concatenate((x, y.results['angles']), axis=1)

        return (x)

    def _search_nonetypes(self, traj):
        """searches for nonetype dihedral angles and returns a list of their indexes"""

        selection_call = 'i.' + self.dihedral + '_selection()'

        index = 0
        nonetype_list = []
        residues = traj.select_atoms("protein")
        for i in residues.residues:
            if not eval(selection_call):
                nonetype_list.append(index)
            index += 1
        return nonetype_list

    @staticmethod
    def _add_bonds(traj):
        """add bonds to protein (only) in a trajectory"""

        protein_res = traj.select_atoms("protein")
        bonds = MDAnalysis.topology.guessers.guess_bonds(
            protein_res, protein_res.positions
        )
        traj.add_TopologyAttr("bonds", bonds)
        return traj

