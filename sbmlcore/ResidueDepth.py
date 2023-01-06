import pathlib

import pandas

from Bio.PDB.ResidueDepth import residue_depth
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.ResidueDepth import get_surface

class ResidueDepth(object):

    def __init__(self, pdb_file, segids=None, offsets=None):

        # check file exists
        assert pathlib.Path(pdb_file).is_file(), "File does not exist!" + pdb_file

        assert not ((segids is not None) and (offsets is not None)), 'cannot specify both segids and offsets - pick one'
        assert ((segids is not None) or (offsets is not None)), 'must specify one of segids or offsets'

        parser = PDBParser()
        structure = parser.get_structure('structure', pdb_file)
        model = structure[0]
        chain_list = [i.get_id() for i in model.get_chains()]

        surface = get_surface(model)

        if segids is not None:
            assert isinstance(segids, list)
            for chain in segids:
                assert chain in chain_list, 'segid ' + chain + ' is not present in the PDB file.'        
        elif offsets is not None:
            assert isinstance(offsets, dict), "Offsets should be specified as a dictionary e.g. offsets = {'A': 3, 'B': -4}"
            for chain in offsets:
                assert chain in chain_list, "Need to specify a segid that exists in pdb!"
                assert isinstance(offsets[chain], int), "Offsets for each segid must be an integer!"

        rows = {'segid': [],\
                'resid': [],\
                'depth': []}


        for segid in segids:

            chain = model[segid]

            for i in chain.get_residues():

                resid = i.id[1]
                
                residue = chain[resid]

                depth = residue_depth(residue, surface)

                if offsets is not None:
                    resid += offsets[segid]

                rows['segid'].append(segid)
                rows['resid'].append(resid)
                rows['depth'].append(depth)

        self.results = pandas.DataFrame.from_dict(rows)

    def _add_feature(self, other):

        assert isinstance(other, pandas.DataFrame)

        assert 'mutation' in other.columns, 'passed dataframe must contain a column called mutation'

        assert 'segid' in other.columns, 'passed dataframe must contain a column called segid containing chain information e.g. A'

        def split_mutation(row):
            return  int(row.mutation[1:-1])

        other['resid'] = other.apply(split_mutation, axis=1)

        other.set_index(['segid',  'resid'], inplace=True)
        self.results.set_index(['segid', 'resid'], inplace=True)

        other = other.join(self.results, how='left')

        other.reset_index(inplace=True)
        self.results.reset_index(inplace=True)

        other.drop(columns = ['resid'], inplace=True)

        return(other)