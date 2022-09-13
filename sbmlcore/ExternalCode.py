import pathlib
import shutil
import subprocess

import pandas
import freesasa
import sbmlcore

# unsure whether to try and do a base class

# class ExternalCode(object):
#
#     def __init__(self, PBBFile):
#


def wrap_angle(angle):
    angle =  angle % 360
    if angle > 180:
        angle -= 360;
    return(angle)


class Stride(object):
    """
    Prediction of secondary structure, SASA, backbone dihedral angles and numbers of hbond donor and acceptors using STRIDE

    Args:
        PDBFile (file): Protein Databank file containing protein coordinates
        offsets (dict): dictionary of form {segid (str): value (int)} where value is 
                the numerical offset between the genetic sequence and the PDB
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

        output = subprocess.getoutput(str(stride)+' -h ' + self.pdb_file)

        rows=[]
        hbond_acc=[]
        hbond_dnr=[]
        for line in output.split('\n'):
            if line[:3] == 'ASG':
                rows.append(line.split()[1:])
            elif line[:3] == 'ACC':
                cols = line.split()
                hbond_acc.append([cols[1], cols[2], cols[3], cols[10]])
            elif line[:3] == 'DNR':
                cols = line.split()
                hbond_dnr.append([cols[1], cols[2], cols[3], cols[10]])

        self.results = pandas.DataFrame(rows, columns=['resname', 'segid', 'resid',\
                                                       'ordinal_resid', 'secondary_structure',\
                                                       'secondary_structure_long', 'phi', 'psi',\
                                                       'residue_sasa', 'pdb_code'])

        acceptors = pandas.DataFrame(hbond_acc, columns=['resname', 'segid', 'resid', 'n_hbond_acceptors'])
        acceptors = acceptors.groupby(['resname', 'segid', 'resid']).count()

        donors = pandas.DataFrame(hbond_dnr, columns=['resname', 'segid', 'resid', 'n_hbond_donors'])
        donors = donors.groupby(['resname', 'segid', 'resid']).count()

        self.results.set_index(['resname', 'segid', 'resid'], inplace=True)
        self.results = self.results.join(acceptors, how='left')
        self.results = self.results.join(donors, how='left')
        self.results.fillna(0, inplace=True)
        self.results.reset_index(inplace=True)

        def short_amino_acid(row):
            return sbmlcore.amino_acid_3to1letter[row.resname]

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
                                     'phi', 'psi', 'residue_sasa', 'n_hbond_acceptors', 'n_hbond_donors']]

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
    The Solvent Accessible Surface Area (SASA) information for each amino acid in a protein. 

    Notes:
        Uses the external FreeSASA Python module to obtain the solvent accessible
        surface areas of each residue in a user-specified pdb file.

    Args:
        PDBFile (file): path to a Protein Databank file containing protein coordinates
        offsets (dict): dictionary of form {segid (str): value (int)} where value is 
                        the numerical offset between the genetic sequence and the PDB
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


        # Split mutation df to create new index in form of segid-resid from mutation
        def split_mutation(row):
            m=row.mutation
            return(int(m[1:-1]))

        other['resid'] = other.apply(split_mutation, axis=1)
        other['id'] = other['segid'] + other['resid'].astype(str)
        other.set_index('id', inplace=True)

        # Create single letter resname column for mutation dataframe
        def resname_1(row):
            m=row.mutation
            return(str(m[0:1]))

        other['resname_1'] = other.apply(resname_1, axis=1)

        # Adds offsets to mutation dataframe
        if self.offsets is not None:
            other["chain_offsets"] = [self.offsets[chain] for chain in other.segid]
        else:
            other["chain_offsets"] = 0

        # Add three letter amino acid to mutation dataframe (needed for FreeSASA input)

        other["resname_3"] = [sbmlcore.amino_acid_1to3letter[resname] for resname in other.resname_1]

        # Adds column for pdb resids (i.e. the resid as given in the pdb which may not be the same as in the mutation df)
        # N.B. Specify the offsets in the same way as you did for structural features class!
        other["pdb_resid"] = other["resid"] - other["chain_offsets"]

        # Creates correct text input for FreeSASA - N.B. includes offsets in  pdb_resid i.e. these new resids should be the same as in the pdb (not the mutation dataframe)
        sele_text = ["%s%i, resi %i and chain %s and resn %s" % (k,i,j,k,l) for i,j,k,l in zip(other.resid, other.pdb_resid, other.segid, other.resname_3)]

        # sele_text = ["%s%i, resi %i and chain %s" % (j,i,i,j) for i,j in zip(other.resid, other.segid)]

        # Obtain SASAs for each residue
        structure = freesasa.Structure(self.pdb_file)
        values = freesasa.calc(structure)
        results = freesasa.selectArea(sele_text, structure, values)
        s = pandas.Series(results)
        b = pandas.DataFrame(s, columns=['SASA'])
        # print(b)

        # Join SASA df to original mutation df
        other = other.join(b, how='left')
        other.reset_index(drop=True, inplace=True)

        other.drop(columns = ['resname_1', 'chain_offsets', 'resname_3', 'pdb_resid', 'resid'], inplace=True)

        return(other)


class SNAP2(object):
    """
    Likelihood of each amino acid mutation affecting the function of the protein.

    Notes:
        Uses the .csv output file from Snap2 https://rostlab.org/services/snap2web/
        If structure contains multiple chains, each one needs to be loaded into the 
        webserver individually, then the .csv files need to be concatenated using 
        csv_segid_concat.ipynb and the segid for each chain correctly assigned. 
        Only then can the SNAP2 class here be used.

    Args:
        CSVFile (.csv file): file output from the SNAP2 webserver
        offsets (dict): dictionary of form {segid (str): value (int)} where value is 
                        the numerical offset between the genetic sequence and the PDB
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

        # create dataframe from .csv file
        snap2_df = pandas.read_csv(self.csv_file)

        # check offsets are correctly specified or raise KeyError

        for segid in snap2_df['segid']:
            if segid not in offsets:
                raise KeyError('Need to specify an offset for ALL segids!')

        # remove entries with less than 80% accuracy
        # N.B. 27/05/22 Removed this feature by changing to remove less than 0%
        # could make this into a user option instead?
        # first need to change 'Expected Accuracy' from strings to floats
        no_percentage = snap2_df['Expected Accuracy'].replace(to_replace='%', value='', regex=True)

        snap2_df['Expected Accuracy'].replace(to_replace='%', value='', regex=True, inplace=True)
        snap2_df['Expected Accuracy'] = snap2_df['Expected Accuracy'].astype('int')

        # turn str into int64 so that inequalities can be used
        series = pandas.to_numeric(no_percentage)

        # no_percentage[ series < 80 ] #creates series with only those rows for which accuracy < 80%
        # no_percentage[series < 80].index #extracts index for each of the rows for which accuracy < 80%

        # remove the entries for which the indices are specified above
        snap2_df.drop(no_percentage[series < 0].index, inplace=True)

        # add offsets column and correct mutation resid column and mutated resname (i.e. the residue change resulting from the mutation)

        def split_mutation_toresname(row):
            return pandas.Series([row.Variant[-1], int(row.Variant[1:-1])])

        snap2_df[["alt_amino_acid", "resid"]] = snap2_df.apply(split_mutation_toresname, axis=1)

        # adds column for offsets
        snap2_df["chain_offsets"] = [offsets[chain] for chain in snap2_df.segid]

        # applies offsets - adds them, as is also the case for Structural Features
        snap2_df["resid"] = snap2_df["resid"] + snap2_df["chain_offsets"]

        # snap2_df.drop(columns='resid', inplace=True)

        snap2_df.rename(columns={'Predicted Effect':'snap2_effect',\
                                 'Score':'snap2_score',\
                                 'Expected Accuracy':'snap2_accuracy',\
                                 'Variant': 'mutation'
                                }, inplace=True)

        self.results = snap2_df


    def _add_feature(self, other):
        """
        Private method to add distances to existing mutation dataframe
        """

        assert isinstance(other, pandas.DataFrame), "You must be adding the extra feature to an existing dataframe!"

        assert 'snap2_score' not in other.columns, "You've already added that feature!"

        assert 'mutation' in other.columns, "Passed dataframe must contain a column called mutation"

        assert 'segid' in other.columns, "Passed dataframe must contain a column called segid containing chain information e.g. A"

        # identifies the original one letter resname and resid (consistent with offset) so that the new feature can be subsequently linked to these
        #def split_mutation(row):
        #    return pandas.Series([row.mutation[0], int(row.mutation[1:-1])])

        #other[['amino_acid', 'mutation_resid']] = other.apply(split_mutation, axis=1)

        # add column for mutated_to_resname into original mutation df (i.e. the residue change as a result of the mutation) and the mutation resid so that the new feature can subsequently be linked to these
        def split_mutation_toresname(row):
            return pandas.Series([row.mutation[-1], int(row.mutation[1:-1])])

        other[['alt_amino_acid', 'resid']] = other.apply(split_mutation_toresname, axis=1)

        # create MultiIndex using segid, resid and amino_acid
        other.set_index(['segid', 'resid', 'alt_amino_acid'], inplace=True)
        self.results.set_index(['segid', 'resid', 'alt_amino_acid'], inplace=True)
        other = other.join(self.results[['snap2_score', 'snap2_accuracy']], how='left')
        other.reset_index(inplace=True) #Removes multi-index

        # remove superfluous columns
        other.drop(columns = ['resid', 'alt_amino_acid'], inplace=True)

        return(other)
