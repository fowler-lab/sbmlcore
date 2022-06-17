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
    """
    Secondary structure and SASA prediction

    Parameters
    ----------
    .pdb file
    offsets - as dictionary of form {segid (str): value (int)}

    Returns
    -------
    dataframe with additional columns for:
    secondary_structure
    secondary_structure_long
    phi
    psi
    residue_sasa
    B
    C
    E
    G
    H
    T
    """

    def __init__(self, PDBFile, offsets=None):

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

        self.results = self.results[['resid', 'amino_acid', 'segid',
                                     'secondary_structure', 'secondary_structure_long',
                                     'phi', 'psi', 'residue_sasa']]

        tmp = pandas.get_dummies(self.results.secondary_structure)
        self.results = self.results.join(tmp, how='left')

        # apply any offsets to the residue numbering
        # as specified in the supplied offsets dict e.g. {'A': 3, 'B': -4}
        # Chain is segid i.e. A, B, C etc.

        def update_resid(row, offsets):
            offset = offsets[row['segid']]
            return row['resid'] + offset

        if offsets is not None:

            assert isinstance(offsets, dict), "Offsets should be specified as a dictionary e.g. offsets = {'A': 3, 'B': -4}"

            assert set(offsets.keys()) == set(self.results.segid.unique()), 'Specified segids in offsets do not match the pdb!'

            for chain in offsets:
                assert isinstance(offsets[chain], int), "Offsets for each segid must be an integer!"

            self.results['resid'] = self.results.apply(update_resid, args=(offsets,), axis=1)



    def _add_feature(self, other, feature_name='all'):

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

    Parameters
    ----------
    .pdb file
    offsets - as dictionary of form {segid (str): value (int)}

    Returns
    -------
    dataframe with additional column for SASA

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
#        structure = freesasa.Structure(self.pdb_file)
#        values = freesasa.calc(structure)
#        area_classes = freesasa.classifyResults(values, structure)

#        print("Total : %.2f A2" % values.totalArea())
#        for key in area_classes:
#            print(key, ": %.2f A2" % area_classes[key])

    def _add_feature(self, other):
        """
        Calculates and adds the SASA for each residue to the existing mutation dataframe (other).
        """
        assert isinstance(other, pandas.DataFrame), "You must be adding the extra feature to an existing dataframe!"

        assert 'mutation' in other.columns, "Passed dataframe must contain a column called mutation"

        assert 'segid' in other.columns, "Passed dataframe must contain a column called segid containing chain information e.g. A"

        # Checks on offset specification
        if self.offsets is not None:
            for chain in self.offsets:
                assert chain in set(other.segid), "Must only specify the segids in the mutation dataframe! And these must also be present in the pdb."
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
        other.reset_index(drop=True, inplace=True)

        other.drop(columns = ['resname_1', 'chain_offsets', 'resname_3', 'pdb_resid', 'resid'], inplace=True)

        return(other)

class SNAP2(object):
    """
    Uses the .csv output file from Snap2 https://rostlab.org/services/snap2web/ which predicts the likelihood of each amino acid mutation affecting the function of the protein.

    Parameters
    ----------
     .csv output from the SNAP2 webserver
             N.B. If structure contains multiple chains, each one needs to be loaded into the webserver individually, then the .csv files need to be concatenated using csv_segid_concat.ipynb and the segid for each chain correctly assigned. Only then can the SNAP2 class here be used.
     offsets - as dictionary of form {segid (str): value (int)}

    Returns
    -------
    dataframe with additional columns for:
    Predicted Effect
    Score
    Expected Accuracy
        N.B. Score is probably the only useful column so dropping the other two columns is advised.
    """
    def __init__(self, CSVFile, offsets=None):

        if not pathlib.Path(CSVFile).is_file():
            raise IOError("Specified CSV file does not exist!")

        self.csv_file = CSVFile

        # apply any offsets to the residue numbering
        # as specified in the supplied offsets dict e.g. {'A': 3, 'B': -4}
        # Chain is segid i.e. A, B, C etc.
        if offsets is not None:
            assert isinstance(offsets, dict), "Offsets should be specified as a dictionary e.g. offsets = {'A': 3, 'B': 4}"
            self.offsets = offsets
        else:
            print("Must supply offsets, and .csv must have segids")
            self.offsets = offsets

        #Create dataframe from .csv file
        snap2_df = pandas.read_csv(self.csv_file)
        #print(snap2_df)

        #Check offsets are correctly specified or raise KeyError

        for segid in snap2_df['segid']:
            if segid not in offsets:
                raise KeyError('Need to specify an offset for ALL segids!')

        #Remove entries with less than 80% accuracy
        #N.B. 27/05/22 Removed this feature by changing to remove less than 0%
        #Could make this into a user option instead?
        #First need to change 'Expected Accuracy' from strings to floats
        no_percentage = snap2_df['Expected Accuracy'].replace(to_replace='%', value='', regex=True)

        #Turn str into int64 so that inequalities can be used
        series = pandas.to_numeric(no_percentage)

        #no_percentage[ series < 80 ] #creates series with only those rows for which accuracy < 80%
        #no_percentage[series < 80].index #extracts index for each of the rows for which accuracy < 80%

        #Remove the entries for which the indices are specified above
        snap2_df.drop(no_percentage[series < 0].index, inplace=True)

        #Add offsets column and correct mutation resid column and mutated resname (i.e. the residue change resulting from the mutation)

        def split_mutation_toresname(row):
            return pandas.Series([row.Variant[-1], int(row.Variant[1:-1])])

        snap2_df[["mutated_to_resname", "resid"]] = snap2_df.apply(split_mutation_toresname, axis=1)


        #Adds column for offsets
        snap2_df["chain_offsets"] = [offsets[chain] for chain in snap2_df.segid]

        #Applies offsets - adds them, as is also the case for Structural Features
        snap2_df["mutation_resid"] = snap2_df["resid"] + snap2_df["chain_offsets"]

        self.snap2_df = snap2_df


    def _add_feature(self, other):
        """
        Adds distances to existing mutation dataframe, and returns new joined dataframe.
        """

        assert isinstance(other, pandas.DataFrame), "You must be adding the extra feature to an existing dataframe!"

        assert 'Score' not in other.columns, "You've already added that feature!"

        assert 'mutation' in other.columns, "Passed dataframe must contain a column called mutation"

        assert 'segid' in other.columns, "Passed dataframe must contain a column called segid containing chain information e.g. A"

        # Identifies the original one letter resname and resid (consistent with offset) so that the new feature can be subsequently linked to these
        #def split_mutation(row):
        #    return pandas.Series([row.mutation[0], int(row.mutation[1:-1])])

        #other[['amino_acid', 'mutation_resid']] = other.apply(split_mutation, axis=1)

        #Add column for mutated_to_resname into original mutation df (i.e. the residue change as a result of the mutation) and the mutation resid so that the new feature can subsequently be linked to these
        def split_mutation_toresname(row):
            return pandas.Series([row.mutation[-1], int(row.mutation[1:-1])])

        other[['mutated_to_resname', 'mutation_resid']] = other.apply(split_mutation_toresname, axis=1)

        #Create MultiIndex using segid, resid and amino_acid
        other.set_index(['segid', 'mutation_resid', 'mutated_to_resname'], inplace=True)
        self.snap2_df.set_index(['segid', 'mutation_resid', 'mutated_to_resname'], inplace=True)

        other = other.join(self.snap2_df, how='left')
        other.reset_index(inplace=True) #Removes multi-index

        #Remove superfluous columns
        other.drop(columns = ['mutation_resid', 'mutated_to_resname', 'Variant', 'resid', 'chain_offsets'], inplace=True)
        #print(other)
        return(other)
