import pandas
import pathlib
import MDAnalysis
import sbmlcore


class StructuralDistances(object):
    """
    Distances between a specified region (i.e. origin) and all amino acids in a protein.

    Args:
        pdb_file (file): path to the Protein DataBank file
        distance_selection (str): the MDAnalysis style selection text that defines the origin
                                  e.g. "resname MG"
        distance_name (str): what you would like to call this distance
        infer_masses (bool): whether to allow MDAnalysis to infer the masses of the atoms in the distance_selection,
                             allows the centre of mass to be calculated. This will generally work for residues, but
                             for small molecules the masses may not be able to be inferred, in which case set this to False
                             and ensure only heavy atoms are in the supplied PDB and it will find the centre of geometry.
        offsets (dict): dictionary of form {segid (str): value (int)} where value is
                        the numerical offset between the genetic sequence and the PDB

    Examples:
        The below will calculate the distance from the Magnesium ion in the M. tuberculosis
        RNA polymerase structure to all amino acids in the protein. Since chain D (rpoB) in the
        PDB file is numbered differently to the gene an offset dictionary is also supplied.

        >>> a = sbmlcore.StructuralDistances('tests/5uh6.pdb','resname MG', 'Mg_distance', offsets = {'A': 0, 'B': 0, 'C': -6})
    """

    def __init__(
        self,
        pdb_file,
        distance_selection,
        distance_name,
        dataset_type="mutation",
        infer_masses=True,
        offsets=None,
    ):

        # check file exists
        assert pathlib.Path(pdb_file).is_file(), "File does not exist!"

        u = MDAnalysis.Universe(pdb_file)

        assert dataset_type in [
            "mutation",
            "amino_acid",
        ], 'dataset_type must be either "mutation" or "amino_acid"!'
        self.dataset_type = dataset_type

        # ensure distance selection is a string
        assert isinstance(
            distance_selection, str
        ), "Distance selection must be a string!"

        # prefer to calculate the centre of mass, but the masses cannot always be inferred from the PDB file
        if infer_masses:
            reference_com = u.select_atoms(distance_selection).center_of_mass()
        else:
            reference_com = u.select_atoms(distance_selection).center_of_geometry()

        # check atom selection exists
        assert (
            u.select_atoms(distance_selection).n_atoms > 0
        ), "Atom selection does not exist! Is your selection using the correct MDAnalysis syntax?"

        # apply any offsets to the residue numbering
        # as specified in the supplied offsets dict e.g. {'A': 3, 'B': -4}
        # Chain is segid i.e. A, B, C etc.
        if offsets is not None:
            assert isinstance(
                offsets, dict
            ), "Offsets should be specified as a dictionary e.g. offsets = {'A': 3, 'B': -4}"
            for chain in offsets:
                assert chain in set(
                    u.residues.segids
                ), "Need to specify a segid that exists in pdb!"
                assert isinstance(
                    offsets[chain], int
                ), "Offsets for each segid must be an integer!"
                chainGroup = u.select_atoms("segid " + chain)
                chainGroup.residues.resids = chainGroup.residues.resids + offsets[chain]

        Ca_all = u.select_atoms("name CA")

        distances = MDAnalysis.lib.distances.distance_array(
            reference_com, Ca_all.positions
        )

        Ca_data = {
            "segid": Ca_all.segids,
            "resid": Ca_all.resids,
            "resname": Ca_all.resnames,
            distance_name: distances[0],
        }

        def one_letter(row):
            return sbmlcore.amino_acid_3to1letter[row.resname]

        results = pandas.DataFrame(Ca_data)
        results["amino_acid"] = results.apply(one_letter, axis=1)

        self.results = results
        self.distance_name = distance_name

    def _add_feature(self, other):
        """
        Private method to add distances to existing mutation dataframe
        """

        assert isinstance(
            other, pandas.DataFrame
        ), "You must be adding the extra feature to an existing dataframe!"

        assert (
            self.distance_name not in other.columns
        ), "You've already added that feature!"

        assert (
            "segid" in other.columns
        ), "Passed dataframe must contain a column called segid containing chain information e.g. A"

        if self.dataset_type == "mutation":

            assert (
                "mutation" in other.columns
            ), "Passed dataframe must contain a column called mutation"

            # Identifies the original one letter resname and resid (consistent with offset) so that the new feature can be subsequently linked to these
            def split_mutation(row):
                return pandas.Series([row.mutation[0], int(row.mutation[1:-1])])

            other[["amino_acid", "resid"]] = other.apply(split_mutation, axis=1)

        elif self.dataset_type == "amino_acid":

            assert (
                "amino_acid" in other.columns
            ), "Passed dataframe must contain a column called amino_acid"

            assert (
                "resid" in other.columns
            ), "Passed dataframe must contain a column called resid"

        # create MultiIndex using segid, resid and amino_acid
        other.set_index(["segid", "resid", "amino_acid"], inplace=True)
        self.results.set_index(["segid", "resid", "amino_acid"], inplace=True)

        other = other.join(self.results, how="left")

        # to check that the chain offsets have been correctly set
        # and/or that structural data is of adequate quality -
        # defining this as the no NaNs should be less than half no rows

        # // divides and rounds DOWN to nearest integer
        half_data = len(other[self.distance_name]) // 2

        total_nans = other[self.distance_name].isna().sum()

        assert (
            total_nans < half_data
        ), "Too many NaNs! Have you defined your offsets correctly?"

        for i in other[self.distance_name]:
            assert isinstance(i, float), "Distances must be floats!"

        other.reset_index(inplace=True)
        self.results.reset_index(inplace=True)

        if self.dataset_type == "mutation":
            other.drop(columns=["amino_acid", "resid", "resname"], inplace=True)
        elif self.dataset_type == "amino_acid":
            other.drop(columns=["resname"], inplace=True)

        return other
