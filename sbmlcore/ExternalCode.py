import pathlib
import shutil
import subprocess

import pandas
# unsure whether to try and do a base class

# class ExternalCode(object):
#
#     def __init__(self, PBBFile):
#
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

        output = subprocess.getoutput(str(stride)+' tests/3pl1.pdb')

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
