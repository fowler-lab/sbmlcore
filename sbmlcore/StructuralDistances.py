import pandas
import pathlib
import MDAnalysis

amino_acid_3to1letter = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

class StructuralDistances(object):

    def __init__(self, pdb_file, distance_selection, distance_name, offsets=None):

        #Check file exists
        assert pathlib.Path(pdb_file).is_file(), "File does not exist!"

        u = MDAnalysis.Universe(pdb_file)
        #print(set(u.residues.segids))

        #Ensure distance selection is a string
        assert isinstance(distance_selection, str), "Distance selection must be a string!"

        reference_com = u.select_atoms(distance_selection).center_of_mass()
        #Check atom selection exists
        assert u.select_atoms(distance_selection).n_atoms > 0, "Atom selection does not exist! Is your selection using the correct MDAnalysis syntax?"

        # apply any offsets to the residue numbering
        # as specified in the supplied offsets dict e.g. {'A': 3, 'B': -4}
        if offsets is not None:
            assert isinstance(offsets, dict), "Offsets should be specified as a dictionary e.g. offsets = {'A': 3, 'B': -4}"
            for chain in offsets:
                assert chain in set(u.residues.segids), "Need to specify a segid that exists in pdb!"
                assert isinstance(offsets[chain], int), "Offsets for each segid must be an integer!"
                chainGroup = u.select_atoms('segid ' + chain)
                chainGroup.residues.resids = chainGroup.residues.resids + offsets[chain]

        Ca_all = u.select_atoms("name CA")

        distances = MDAnalysis.lib.distances.distance_array(reference_com, Ca_all.positions)

        Ca_data = {'segid': Ca_all.segids, 'resid': Ca_all.resids,
                   'resname': Ca_all.resnames, distance_name: distances[0]}

        def one_letter(row):
            return(amino_acid_3to1letter[row.resname])

        results = pandas.DataFrame(Ca_data)
        results['amino_acid'] = results.apply(one_letter, axis=1)

        self.results = results
        self.distance_name = distance_name


    def add_feature(self, other):

        assert isinstance(other, pandas.DataFrame), "You must be adding the extra feature to an existing dataframe!"

        assert self.distance_name not in other.columns, "You've already added that feature!"

        assert 'mutation' in other.columns, "Passed dataframe must contain a column called mutation"

        assert 'segid' in other.columns, "Passed dataframe must contain a column called segid containing chain information e.g. A"

        # Identifies the original one letter resname and resid (consistent with offset) so that the new feature can be subsequently linked to these
        def split_mutation(row):
            return pandas.Series([row.mutation[0], int(row.mutation[1:-1])])

        other[['amino_acid', 'resid']] = other.apply(split_mutation, axis=1)

        #Create MultiIndex using segid, resid and amino_acid
        other.set_index(['segid', 'resid', 'amino_acid'], inplace=True)
        self.results.set_index(['segid', 'resid', 'amino_acid'], inplace=True)

        other = other.join(self.results, how='left')
        #print(self.distance_name)
        #print(other)
        assert isinstance(other[self.distance_name], float), "Offset has been incorrectly specified such that they do not align with initial mutation resids!"
        #Want to check that all the distance entries in the newly added distances column are all floats - i.e. there are no NaNs (or weird things like ints or strings)

        other.reset_index(inplace=True)
        self.results.reset_index(inplace=True)

        other.drop(columns = ['amino_acid', 'resid', 'resname'], inplace=True)

        return(other)
