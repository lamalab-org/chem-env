import os
from modal import App, fastapi_endpoint
from pydantic import BaseModel
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


class IupacNameInput(BaseModel):
    smiles: str
    timeout: int = 60


converters_name = os.getenv("CONVERTERS_NAME", "")
if converters_name and not converters_name.startswith("-"):
    converters_name = f"-{converters_name}"

# Create the app
converters_app = App(f"chemenv_converters{converters_name}")


@converters_app.function(image=_converters_image, timeout=60)
@fastapi_endpoint(method="POST")
async def get_iupac_name(body: IupacNameInput) -> str:
    """
    Get the IUPAC name of a molecule from its SMILES string.

    Args:
        smiles (str): The SMILES string of the molecule.
        timeout (int): The timeout in seconds for the request.

    Returns:
        str: The IUPAC name of the molecule.

    Raises:
        ValueError: If the conversion fails.
    """
    converter = _Smiles2Name(body.smiles, body.timeout)
    try:
        name = await converter.get_name()
        return name
    except Exception as e:
        raise ValueError(f"Error: {e}") from e


@converters_app.function(image=_converters_image)
async def get_smiles_from_name(name: str, timeout: int = 10) -> str:
    """
    Get the SMILES string of a molecule from its IUPAC name.

    Args:
        name (str): The IUPAC name of the molecule.
        timeout (int): The timeout in seconds for the request.

    Returns:
        str: The SMILES string of the molecule.

    Raises:
        ValueError: If the conversion fails.
    """
    converter = _Name2Smiles(name, timeout)
    try:
        smiles = await converter.get_smiles()
        return smiles
    except Exception as e:
        raise ValueError(f"Error converting name to SMILES: {e}") from e


@converters_app.function(image=_converters_image)
def convert_to_selfies(smiles: str) -> str:
    """
    Convert SMILES to SELFIES encoding.

    Args:
        smiles (str): The SMILES string to convert.

    Returns:
        str: The SELFIES encoding of the molecule.
    """
    return _smiles_to_selfies(smiles)


@converters_app.function(image=_converters_image)
def convert_to_deepsmiles(smiles: str) -> str:
    """
    Convert SMILES to DeepSMILES encoding.

    Args:
        smiles (str): The SMILES string to convert.

    Returns:
        str: The DeepSMILES encoding of the molecule.
    """
    return _smiles_to_deepsmiles(smiles)


@converters_app.function(image=_converters_image)
def convert_to_inchi(smiles: str) -> str:
    """
    Convert SMILES to InChI.

    Args:
        smiles (str): The SMILES string to convert.

    Returns:
        str: The InChI encoding of the molecule.
    """
    try:
        return _smiles_to_inchi(smiles)
    except Exception as e:
        raise ValueError(f"Error converting SMILES to InChI: {e}") from e


@converters_app.function(image=_converters_image)
def convert_to_inchikey(smiles: str) -> str:
    """
    Convert SMILES to InChIKey.

    Args:
        smiles (str): The SMILES string to convert.

    Returns:
        str: The InChIKey encoding of the molecule.
    """
    try:
        return _smiles_to_inchikey(smiles)
    except Exception as e:
        raise ValueError(f"Error converting SMILES to InChIKey: {e}") from e


@converters_app.function(image=_safe_image)
def convert_to_safe(smiles: str) -> str:
    """
    Convert SMILES to SAFE encoding.

    Args:
        smiles (str): The SMILES string to convert.

    Returns:
        str: The SAFE encoding of the molecule.
    """
    return _smiles_to_safe(smiles)


@converters_app.function(image=_converters_image)
def selfies_to_smiles(selfies: str) -> str:
    """
    Convert SELFIES to SMILES.

    Args:
        selfies (str): The SELFIES encoding to convert.

    Returns:
        str: The SMILES string of the molecule.
    """
    return _selfies_to_smiles(selfies)


@converters_app.function(image=_converters_image)
def inchi_to_smiles(inchi: str) -> str:
    """
    Convert InChI to SMILES.

    Args:
        inchi (str): The InChI encoding to convert.

    Returns:
        str: The SMILES string of the molecule.
    """
    try:
        return _inchi_to_smiles(inchi)
    except Exception as e:
        raise ValueError(f"Error converting InChI to SMILES: {e}") from e


@converters_app.function(image=_converters_image)
def deepsmiles_to_smiles(deepsmiles: str) -> str:
    """
    Convert DeepSMILES to SMILES.

    Args:
        deepsmiles (str): The DeepSMILES encoding to convert.

    Returns:
        str: The SMILES string of the molecule.
    """
    try:
        return _deepsmiles_to_smiles(deepsmiles)
    except Exception as e:
        raise ValueError(f"Error converting DeepSMILES to SMILES: {e}") from e
