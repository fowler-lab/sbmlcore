import pathlib
import shutil
import subprocess

import pandas
# unsure whether to try and do a base class

# class ExternalCode(object):
#
#     def __init__(self, PBBFile):
#

amino_acid_3to1letter = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

def wrap_angle(angle):
    angle =  angle % 360
    if angle > 180:
        angle -= 360;
    return(angle)


class Stride(object):

    def __init__(self, PDBFile):

        assert pathlib.Path(PDBFile).is_file(), "specified PDB file does not exist!"

        self.pdb_file = PDBFile

        if pathlib.Path("./stride").exists():
            stride = pathlib.Path("./stride").resolve()

        # or if there is one in the $PATH use that one
        elif shutil.which('stride') is not None:
             stride = pathlib.Path(shutil.which('stride'))

        output = subprocess.getoutput(str(stride)+' ' + self.pdb_file)

        rows=[]
        for line in output.split('\n'):
            if line[:3] == 'ASG':
                rows.append(line.split()[1:])

        self.results = pandas.DataFrame(rows, columns=['resname', 'segid', 'resid',\
                                                       'ordinal_resid', 'secondary_structure',\
                                                       'secondary_structure_long', 'phi', 'psi',\
                                                       'residue_sasa', 'pdb_code'])

        def short_amino_acid(row):
            return amino_acid_3to1letter[row.resname]

        self.results['amino_acid'] = self.results.apply(short_amino_acid, axis=1)

        self.results = self.results.astype({'resid': 'int',\
                                            'phi': 'float',\
                                            'psi': 'float',\
                                            'residue_sasa': 'float'})
        def correct_torsions(row):
            return pandas.Series([wrap_angle(row.phi), wrap_angle(row.psi)])

        self.results[['phi','psi']] = self.results.apply(correct_torsions, axis=1)

        self.results = self.results[['resid', 'amino_acid', 'resname', 'segid',
                                     'secondary_structure', 'secondary_structure_long',
                                     'phi', 'psi', 'residue_sasa']]

    def add_feature(self, other, feature_name='all'):

        assert isinstance(other, pandas.DataFrame)

        assert 'mutation' in other.columns, 'passed dataframe must contain a column called mutations'

        assert 'segid' in other.columns, 'passed dataframe must contain a column called segid containing chain information e.g. A'

        assert feature_name not in other.columns, 'trying to add a column that is already in the dataset!'

        if feature_name != 'all':
            assert feature_name in self.results.columns, 'supplied feature_name not supplied by STRIDE'

        def split_mutation(row):
            return pandas.Series([row.mutation[0], int(row.mutation[1:-1])])

        other[['amino_acid', 'resid']] = other.apply(split_mutation, axis=1)

        other.set_index(['segid', 'amino_acid', 'resid'], inplace=True)
        self.results.set_index(['segid', 'amino_acid', 'resid'], inplace=True)

        if feature_name == 'all':
            other = other.join(self.results, how='left')
        else:
            other = other.join(self.results[[feature_name]], how='left')

        other.reset_index(inplace=True)
        self.results.reset_index(inplace=True)

        other.drop(columns = ['amino_acid', 'resid'], inplace=True)

        return(other)
