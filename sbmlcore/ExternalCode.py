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

    def __init__(self, PDBFile, offsets=None):

        if not pathlib.Path(PDBFile).is_file():
            raise IOError("Specified PDB file does not exist!")

        self.pdb_file = PDBFile

        # apply any offsets to the residue numbering
        # as specified in the supplied offsets dict e.g. {'A': 3, 'B': -4}
        # Chain is segid i.e. A, B, C etc.
        if offsets is not None:
            assert isinstance(offsets, dict), "Offsets should be specified as a dictionary e.g. offsets = {'A': 3, 'B': -4}"
            self.offsets = offsets
        else:
            print("offsets = None, are you sure?")
            self.offsets = offsets
#        structure = freesasa.Structure(self.pdb_file)
#        values = freesasa.calc(structure)
#        area_classes = freesasa.classifyResults(values, structure)

#        print("Total : %.2f A2" % values.totalArea())
#        for key in area_classes:
#            print(key, ": %.2f A2" % area_classes[key])

    def add_feature(self, other):
        """
        Calculates and adds the SASA for each residue to the existing mutation dataframe (other).

        Arguments: existing dataframe
        e.g. if a = sbmlcore.FreeSASA('path_pdbfile'),
        use new_df = a.add_feature(existing_df)
        """
        assert isinstance(other, pandas.DataFrame), "You must be adding the extra feature to an existing dataframe!"

        assert 'mutation' in other.columns, "Passed dataframe must contain a column called mutation"

        assert 'segid' in other.columns, "Passed dataframe must contain a column called segid containing chain information e.g. A"

        # Checks on offset specification
        if self.offsets is not None:
            for chain in self.offsets:
                assert chain in set(other.segid), "Need to specify a segid that exists in pdb!"
                assert isinstance(self.offsets[chain], int), "Offsets for each segid must be an integer!"
        else:
            pass


    #Split mutation df to create new index in form of segid-resid from mutation
        def split_mutation(row):
            m=row.mutation
            return(int(m[1:-1]))

        other['resid'] = other.apply(split_mutation, axis=1)
        other['id'] = other['segid'] + other['resid'].astype(str)
        other.set_index('id', inplace=True)

    #Create single letter resname column for mutation dataframe
        def resname_1(row):
            m=row.mutation
            return(str(m[0:1]))

        other['resname_1'] = other.apply(resname_1, axis=1)



        #Adds offsets to mutation dataframe
        if self.offsets is not None:
            other["chain_offsets"] = [self.offsets[chain] for chain in other.segid]
        else:
            other["chain_offsets"] = 0

        #Add three letter amino acid to mutation dataframe (needed for FreeSASA input)
        amino_acid_onetothreeletter = {'C': 'CYS', 'D': 'ASP', 'S': 'SER', 'Q': 'GLN', 'K': 'LYS',
        'I': 'ILE', 'P': 'PRO', 'T': 'THR', 'F': 'PHE', 'N': 'ASN',
        'G': 'GLY', 'H': 'HIS', 'L': 'LEU', 'R': 'ARG', 'W': 'TRP',
        'A': 'ALA', 'V': 'VAL', 'E': 'GLU', 'Y': 'TYR', 'M': 'MET'}

        other["resname_3"] = [amino_acid_onetothreeletter[resname] for resname in other.resname_1]

        #Adds column for pdb resids (i.e. the resid as given in the pdb which may not be the same as in the mutation df)
        #N.B. Specify the offsets in the same way as you did for structural features class!
        other["pdb_resid"] = other["resid"] - other["chain_offsets"]

        #Creates correct text input for FreeSASA - N.B. includes offsets in  pdb_resid i.e. these new resids should be the same as in the pdb (not the mutation dataframe)
        sele_text = ["%s%i, resi %i and chain %s and resn %s" % (k,i,j,k,l) for i,j,k,l in zip(other.resid, other.pdb_resid, other.segid, other.resname_3)]

        #sele_text = ["%s%i, resi %i and chain %s" % (j,i,i,j) for i,j in zip(other.resid, other.segid)]

        #Obtain SASAs for each residue
        structure = freesasa.Structure(self.pdb_file)
        values = freesasa.calc(structure)
        results = freesasa.selectArea(sele_text, structure, values)
        s = pandas.Series(results)
        b = pandas.DataFrame(s, columns=['SASA'])
        #print(b)

        #Join SASA df to original mutation df
        other = other.join(b, how='left')

        return(other)

class SNAP2(object):
    """
    Uses the .csv output file from Snap2 https://rostlab.org/services/snap2web/ which predicts the likelihood of each amino acid mutation affecting the function of the protein.

    Arguments: .csv output from the SNAP2 webserver
    """
    def __init__(self, CSVFile, offsets=None):

        if not pathlib.Path(CSVFile).is_file():
            raise IOError("Specified CSV file does not exist!")

        self.csv_file = CSVFile

        # apply any offsets to the residue numbering
        # as specified in the supplied offsets dict e.g. {'A': 3, 'B': -4}
        # Chain is segid i.e. A, B, C etc.
        if offsets is not None:
            assert isinstance(offsets, dict), "Offsets should be specified as a dictionary e.g. offsets = {'A': 3, 'B': -4}"
            self.offsets = offsets
        else:
            print("offsets = None, are you sure?")
            self.offsets = offsets

        #Create dataframe from .csv file
        snap2_df = pandas.read_csv(self.csv_file)
        print(snap2_df)
