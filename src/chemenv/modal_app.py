from modal import App, Image
from chemenv.tools.cheminformatics import (
    get_tanimoto_similarity as _get_tanimoto_similarity,
    get_number_of_topologically_distinct_atoms as _get_topologically_distinct_atoms,
    get_element_info as _get_element_info,
)
from chemenv.tools.converters import (
    _converters_image,
    _safe_image,
    _Smiles2Name,
    _Name2Smiles,
    _smiles_to_selfies,
    _smiles_to_deepsmiles,
    _smiles_to_inchi,
    _smiles_to_inchikey,
    _smiles_to_safe,
    _selfies_to_smiles,
    _inchi_to_smiles,
    _deepsmiles_to_smiles,
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


@app.function(image=_converters_image)
async def get_iupac_name(smiles: str) -> str:
    converter = _Smiles2Name(smiles)
    name = await converter.get_name()
    if name is None:
        return ""
    return name


@app.function(image=_converters_image)
async def get_smiles_from_name(name: str) -> str:
    converter = _Name2Smiles(name)
    smiles = await converter.get_smiles()
    if smiles is None:
        return ""
    return smiles


@app.function(image=_converters_image)
def convert_to_selfies(smiles: str) -> str:
    """Convert SMILES to SELFIES encoding."""
    return _smiles_to_selfies(smiles)


@app.function(image=_converters_image)
def convert_to_deepsmiles(smiles: str) -> str:
    """Convert SMILES to DeepSMILES encoding."""
    return _smiles_to_deepsmiles(smiles)


@app.function(image=_converters_image)
def convert_to_inchi(smiles: str) -> str:
    """Convert SMILES to InChI."""
    return _smiles_to_inchi(smiles)


@app.function(image=_converters_image)
def convert_to_inchikey(smiles: str) -> str:
    """Convert SMILES to InChIKey."""
    return _smiles_to_inchikey(smiles)


@app.function(image=_safe_image)
def convert_to_safe(smiles: str) -> str:
    """Convert SMILES to SAFE encoding."""
    return _smiles_to_safe(smiles)


@app.function(image=_converters_image)
def selfies_to_smiles(selfies: str) -> str:
    """Convert SELFIES to SMILES."""
    return _selfies_to_smiles(selfies)


@app.function(image=_converters_image)
def inchi_to_smiles(inchi: str) -> str:
    """Convert InChI to SMILES."""
    return _inchi_to_smiles(inchi)


@app.function(image=_converters_image)
def deepsmiles_to_smiles(deepsmiles: str) -> str:
    """Convert DeepSMILES to SMILES."""
    return _deepsmiles_to_smiles(deepsmiles)
