"""get identifier from CAS number"""

from urllib import request, error, parse
from time import sleep
import random

def get_identifier_from_cas(cas, identifier):
    url = 'https://cactus.nci.nih.gov/chemical/structure/' + parse.quote(cas) + '/' + str(identifier)
    sleep(random.uniform(0, 0.5))
    try:
        response = request.urlopen(url)
        stdinchikey = [line.strip().decode("utf-8") for line in response]
        if len(stdinchikey) > 6:
            stdinchikey = None
    except error.HTTPError as err:
        stdinchikey = None
    except error.URLError as err:
        stdinchikey = None
    except TimeoutError as err:
        stdinchikey = None
    if stdinchikey is not None and len(stdinchikey) == 1:
        stdinchikey = ''.join(stdinchikey)
    return stdinchikey
cactus_smiles = get_identifier_from_cas('50-78-2', 'smiles')
cactus_stdinchikey = get_identifier_from_cas('50-78-2', 'stdinchikey')

"""get cid from smiles"""
from urllib import request, error, parse
def get_identifier_from_smile(smile, identifier):
    url = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/'+smile+'/'+identifier+'/TXT'
    try:
        response = request.urlopen(url)
        cid = [line.strip().decode("utf-8") for line in response]
        cid = int(cid[0])
    except:
        cid = None
        print('error')
    return cid

"""get cid from smiles using pubchempy"""
import pandas as pd
import pubchempy as pcp

pcp.get_cids('CN(C)N=NC1=C(NC=N1)C(=O)N','smiles')[0]

"""get smiles from CID"""
from pubchempy import get_compounds, Compound
comp = Compound.from_cid(1423)
print(comp.canonical_smiles)
# can also get compound inchi, inchikey, molecular weight and formula etc


"""get CID from name"""
from urllib import request, error, parse
from time import sleep
import random

def get_identifier_from_name(name, identifier):
    sleep(random.uniform(0, 0.05))
    url = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/'+name+'/cids/TXT'
    try:
        response = request.urlopen(url)
        cid = [line.strip().decode("utf-8") for line in response]
        cid = int(cid[0])
    except:
        cid = None
        print('error')
    return cid


""" get CID from CAS"""
from urllib import request, error, parse
from time import sleep
import random
import os
import pandas as pd

def ifCas(compound):
    """ Returns List of CIDS
    compound (str): A Cas registery number or common name identifier
    """
    PUBCHEM_BASE = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/'
    sleep(random.uniform(0, 0.05))
    try:
        url = PUBCHEM_BASE+"compound/name/"+compound+"/cids/TXT"
        response = request.urlopen(url)
        for line in response:
            CIDS = line.strip().split()
            CIDS = [x.decode('utf-8') for x in CIDS]
        return CIDS
    except:
        return [None]

"""check whether two compounds are duplicates"""

from rdkit import Chem
m = Chem.MolFromSmiles('OC(=O)CON=Cc1nc2ccccc2n1[C@H]3C[C@H]4CCC[C@@H](C3)N4[C@H]5C[C@@H]6CCC[C@@H](C6)C5')
n = Chem.MolFromSmiles('OC(=O)CON=Cc1nc2ccccc2[n]1C3CC4CCCC(C3)N4C5CC6CCCC(C6)C5')
if (m is not None) and (n is not None):
    if m.HasSubstructMatch(n) and n.HasSubstructMatch(m):
        print("m and n are duplicates")


"""check whether compound is organic"""
from rdkit.Chem import MolFromSmiles
organic_mols = [mol for mol in all_mols if mol.HasSubstructMatch(MolFromSmiles('C'))]
