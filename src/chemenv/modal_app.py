from modal import App, Image
from chemenv.tools.cheminformatics import (
    get_tanimoto_similarity as _get_tanimoto_similarity,
    get_number_of_topologically_distinct_atoms as _get_topologically_distinct_atoms,
    get_element_info as _get_element_info,
)
from chemenv.tools.pubchem import (
    PubChem,
    pubchem_image as _pubchem_image,
)
from chemenv.tools.pubchem import (
    Smiles2Name as _Smiles2Name,
    converters_image as _converters_image,
)
import os

# Define the images
rdkit_image = (
    Image.debian_slim(python_version="3.12").pip_install("rdkit").pip_install("numpy")
)
mendeleev_image = Image.debian_slim().pip_install("mendeleev")

chemenv_name = os.getenv("CHEMENV_NAME", "")
if chemenv_name and not chemenv_name.startswith("-"):
    chemenv_name = f"-{chemenv_name}"

# Create the app
app = App(f"chemenv{chemenv_name}")


@app.function(image=rdkit_image)
def get_tanimoto_similarity(*args, **kwargs):
    return _get_tanimoto_similarity(*args, **kwargs)


@app.function(image=rdkit_image)
def get_number_of_topologically_distinct_atoms(*args, **kwargs):
    return _get_topologically_distinct_atoms(*args, **kwargs)


@app.function(image=mendeleev_image)
def get_element_info(*args, **kwargs):
    return _get_element_info(*args, **kwargs)


@app.function(image=_pubchem_image)
async def get_number_atoms(compound_id: str) -> int:
    pubchem = await PubChem.create(compound_id)
    number_atoms = await pubchem._get_number_atoms()
    if number_atoms is None:
        raise ValueError("No number of atoms found")
    return number_atoms


@app.function(image=_pubchem_image)
async def get_isomeric_smiles(compound_id: str) -> str:
    pubchem = await PubChem.create(compound_id)
    isomeric_smiles = await pubchem._get_isomeric_smiles()
    if isomeric_smiles is None:
        raise ValueError("No isomeric SMILES found")
    return isomeric_smiles


@app.function(image=_pubchem_image)
async def get_canonical_smiles(compound_id: str) -> str:
    pubchem = await PubChem.create(compound_id)
    smiles = await pubchem._get_canonical_smiles()
    if smiles is None:
        raise ValueError("No canonical SMILES found")
    return smiles


@app.function(image=_pubchem_image)
async def get_compound_mass(compound_id: str) -> float:
    pubchem = await PubChem.create(compound_id)
    molecular_weight = await pubchem._get_compound_mass()
    if molecular_weight is None:
        raise ValueError("No molecular weight found")
    return molecular_weight


@app.function(image=_pubchem_image)
async def get_compound_charge(compound_id: str) -> int:
    pubchem = await PubChem.create(compound_id)
    charge = await pubchem._get_compound_charge()
    if charge is None:
        raise ValueError("No charge found")
    return charge


@app.function(image=_pubchem_image)
async def get_compound_formula(compound_id: str) -> str:
    pubchem = await PubChem.create(compound_id)
    formula = await pubchem._get_compound_formula()
    if formula is None:
        raise ValueError("No formula found")
    return formula


@app.function(image=_pubchem_image)
async def get_number_isomers(compound_id: str) -> int:
    pubchem = await PubChem.create(compound_id)
    number_isomers = await pubchem._get_number_isomers()
    if number_isomers is None:
        raise ValueError("No isomers found")
    return number_isomers


@app.function(image=_pubchem_image)
async def get_compound_isomers(*args, **kwargs):
    pubchem = await PubChem.create(*args, **kwargs)
    data = await pubchem._get_compound_isomers()
    if data is None:
        raise ValueError("No isomers found")
    return data


@app.function(image=_pubchem_image)
async def get_number_heavy_atoms(compound_id: str) -> int:
    pubchem = await PubChem.create(compound_id)
    number_heavy_atoms = await pubchem._get_number_heavy_atoms()
    if number_heavy_atoms is None:
        raise ValueError("No heavy atoms found")
    return number_heavy_atoms


@app.function(image=_pubchem_image)
async def _get_number_chiral_atoms(compound_id: str) -> int:
    pubchem = await PubChem.create(compound_id)
    number_chiral_atoms = await pubchem._get_number_chiral_atoms()
    if number_chiral_atoms is None:
        raise ValueError("No chiral atoms found")
    return number_chiral_atoms


@app.function(image=_converters_image)
async def get_iupac_name(smiles: str) -> str:
    converter = _Smiles2Name(smiles)
    name = await converter.get_name()
    if name is None:
        return ""
    return name
