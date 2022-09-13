import pandas
import pathlib
import MDAnalysis
import sbmlcore
 
class TempFactors(object):
    """
    Distances between a specified region (i.e. origin) and all amino acids in a protein.

    Args:
        pdb_file (file): path to the Protein DataBank file
        distance_selection (str): the MDAnalysis style selection text that defines the origin
                                  e.g. "resname MG"
        distance_name (str): what you would like to call this distance
        offsets (dict): dictionary of form {segid (str): value (int)} where value is 
                        the numerical offset between the genetic sequence and the PDB

    Examples:
        The below will calculate the distance from the Magnesium ion in the M. tuberculosis
        RNA polymerase structure to all amino acids in the protein. Since chain D (rpoB) in the
        PDB file is numbered differently to the gene an offset dictionary is also supplied.

        >>> a = sbmlcore.StructuralDistances('tests/5uh6.pdb','resname MG', 'Mg_distance', offsets = {'A': 0, 'B': 0, 'C': -6})
    """

    def __init__(self, pdb_file, offsets=None):

        # check file exists
        assert pathlib.Path(pdb_file).is_file(), "File does not exist!"

        u = MDAnalysis.Universe(pdb_file)

        # apply any offsets to the residue numbering
        # as specified in the supplied offsets dict e.g. {'A': 3, 'B': -4}
        # Chain is segid i.e. A, B, C etc.
        if offsets is not None:
            assert isinstance(offsets, dict), "Offsets should be specified as a dictionary e.g. offsets = {'A': 3, 'B': -4}"
            for chain in offsets:
                assert chain in set(u.residues.segids), "Need to specify a segid that exists in pdb!"
                assert isinstance(offsets[chain], int), "Offsets for each segid must be an integer!"
                chainGroup = u.select_atoms('segid ' + chain)
                chainGroup.residues.resids = chainGroup.residues.resids + offsets[chain]

        Ca_all = u.select_atoms("name CA")

        Ca_data = {'segid': Ca_all.segids, 'resid': Ca_all.resids,
                   'resname': Ca_all.resnames, 'temp_factor': Ca_all.tempfactors}

        def one_letter(row):
            return(sbmlcore.amino_acid_3to1letter[row.resname])

        results = pandas.DataFrame(Ca_data)
        results['amino_acid'] = results.apply(one_letter, axis=1)

        self.results = results


    def _add_feature(self, other):
        """
        Private method to add temp_factors to existing mutation dataframe
        """

        assert isinstance(other, pandas.DataFrame), "You must be adding the extra feature to an existing dataframe!"

        assert 'temp_factor' not in other.columns, "You've already added that feature!"

        assert 'mutation' in other.columns, "Passed dataframe must contain a column called mutation"

        assert 'segid' in other.columns, "Passed dataframe must contain a column called segid containing chain information e.g. A"

        # Identifies the original one letter resname and resid (consistent with offset) so that the new feature can be subsequently linked to these
        def split_mutation(row):
            return pandas.Series([row.mutation[0], int(row.mutation[1:-1])])

        other[['amino_acid', 'resid']] = other.apply(split_mutation, axis=1)

        # create MultiIndex using segid, resid and amino_acid
        other.set_index(['segid', 'resid', 'amino_acid'], inplace=True)
        self.results.set_index(['segid', 'resid', 'amino_acid'], inplace=True)

        other = other.join(self.results, how='left')

        other.reset_index(inplace=True)
        self.results.reset_index(inplace=True)

        other.drop(columns=['amino_acid', 'resid', 'resname'], inplace=True)

        return(other)
