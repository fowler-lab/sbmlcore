import pathlib

import pandas
import numpy
import MDAnalysis


class TrajectoryDistances(object):
    """Calculate and add distances from a molecular dynamics trajectory.

    Parameters
    ----------

    Returns
    -------

    """

    def __init__(self, pdb_file, trajectory_list, distance_selection, distance_name, distance_type='mean', offsets=None, start_time=None, end_time=None):

        # check files exists
        assert pathlib.Path(pdb_file).is_file(), "File does not exist!"

        assert isinstance(trajectory_list, list), 'trajectories not provided as a list'
        for i in trajectory_list:
            assert pathlib.Path(i).is_file(), "File does not exist!"

        assert distance_type in ['mean'], 'provided distance_type not recognised!'
        assert isinstance(distance_selection, str), "Distance selection must be a string!"
        assert isinstance(distance_name, str), "Distance name must be a string!"

        if start_time is not None:
            assert start_time>0
            assert isinstance(start_time, float)
        if end_time is not None:
            assert end_time>0
            assert isinstance(end_time, float)
            if start_time is not None:
                assert start_time < end_time

        first_pass = True

        for trajectory in trajectory_list:

            u = MDAnalysis.Universe(pdb_file, trajectory)

            reference_com = u.select_atoms(distance_selection).center_of_mass()

            # check atom selection exists
            assert u.select_atoms(distance_selection).n_atoms > 0, "Atom selection does not exist! Is your selection using the correct MDAnalysis syntax?"

            Ca_all = u.select_atoms("name CA")

            for ts in u.trajectory:

                if ts.frame == 0:
                    continue

                if start_time is not None and ts.time < start_time:
                    continue

                if end_time is not None and ts.time > end_time:
                    continue

                distances = MDAnalysis.lib.distances.distance_array(reference_com, Ca_all.positions)

                if first_pass:
                    distance_array = distances
                    first_pass = False
                else:
                    distance_array = numpy.concatenate([distance_array, distances])

            if distance_type == 'mean':
                distances = numpy.mean(distance_array, axis=0)

            print(distances)
