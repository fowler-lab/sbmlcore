import pathlib
import shutil
import subprocess

import pandas
import freesasa
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

        else:
            raise IOError("No stride installed!")

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

        tmp = pandas.get_dummies(self.results.secondary_structure)
        self.results = self.results.join(tmp, how='left')

    def add_feature(self, other, feature_name='all'):

        assert isinstance(other, pandas.DataFrame)

        assert 'mutation' in other.columns, 'passed dataframe must contain a column called mutation'

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


class FreeSASA(object):
    """
    Uses external FreeSASA Python module to obtain the solvent accessible
    surface areas of each residue in a user-specified pdb file.

    Takes one argument: path to pdb file
    e.g. a = sbmlcore.FreeSASA('path_pdbfile').

    Functions:
    add_feature - adds the SASA for each residue to the existing mutation dataframe
    """

    def __init__(self, PDBFile):

        if not pathlib.Path(PDBFile).is_file():
            raise IOError("Specified PDB file does not exist!")

        self.pdb_file = PDBFile

#        structure = freesasa.Structure(self.pdb_file)
#        values = freesasa.calc(structure)
#        area_classes = freesasa.classifyResults(values, structure)

#        print("Total : %.2f A2" % values.totalArea())
#        for key in area_classes:
#            print(key, ": %.2f A2" % area_classes[key])

    def add_feature(self, other):
        """
        Calculates and adds the SASA for each residue to the existing mutation dataframe.

        Arguments: existing dataframe
        e.g. if a = sbmlcore.FreeSASA('path_pdbfile'),
        use new_df = a.add_feature(existing_df)
        """
        assert isinstance(other, pandas.DataFrame), "You must be adding the extra feature to an existing dataframe!"

        assert 'mutation' in other.columns, "Passed dataframe must contain a column called mutation"

        assert 'segid' in other.columns, "Passed dataframe must contain a column called segid containing chain information e.g. A"

    #Split mutation df to create new index in form of segid-resid from mutation
        def split_mutation(row):
            m=row.mutation
            return(int(m[1:-1]))

        other['resid'] = other.apply(split_mutation, axis=1)
        other['id'] = other['segid'] + other['resid'].astype(str)
        other.set_index('id', inplace=True)

        #Creates correct text input for FreeSASA
        sele_text = ["%s%i, resi %i and chain %s" % (j,i,i,j) for i,j in zip(other.resid, other.segid)]

        #Obtain SASAs for each residue
        structure = freesasa.Structure(self.pdb_file)
        values = freesasa.calc(structure)
        results = freesasa.selectArea(sele_text, structure, values)
        s = pandas.Series(results)
        b = pandas.DataFrame(s, columns=['surface_area'])
        #print(b)

        #Join SASA df to original mutation df
        other = other.join(b, how='left')

        return(other)
#        self.results = pandas.DataFrame(rows, columns=['resname', 'segid', 'resid',\
#                                                       'ordinal_resid', 'secondary_structure',\
#                                                       'secondary_structure_long', 'phi', 'psi',\
#                                                       'residue_sasa', 'pdb_code'])
