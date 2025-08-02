import os
from modal import App, Volume
from typing import Any
from chemenv.tools.rxn_schema_processing import (
    decimer_image,
    _decimer_extractor,
    rxnscribe_image,
    _molscribe_extractor,
    _rxnscribe_extractor,
)
from chemenv.tools.rxn_utils import (
    rxn_mapper_image,
    _get_rxn_mapper_confidence,
    rxn_utils_image,
    _unify_rxn_smiles,
    _canonicalize_rxn_smiles,
)

MINUTES = 60  # 60 seconds

rxn_name = os.getenv("RXNENV_NAME", "")
if rxn_name and not rxn_name.startswith("-"):
    rxn_name = f"-{rxn_name}"

# Define the volume to safe checkpoints
hf_cache_vol = Volume.from_name("huggingface-cache", create_if_missing=True)

# Create the app
rxn_app = App(f"chemenv_rxnenv{rxn_name}")


@rxn_app.function(image=decimer_image)
def molecule_image_extraction_decimer(*args, **kwargs) -> str:
    """
    Predicts the SMILES from the image using DECIMER models.
    https://github.com/Kohulan/DECIMER-Image_Transformer

    Args:
        image: str: Path to the image file.

    Returns:
        str: SMILES predicted.
    """
    return _decimer_extractor(*args, **kwargs)


@rxn_app.function(
    image=rxnscribe_image,
    gpu="A10G",
    volumes={
        "/root/.cache/huggingface": hf_cache_vol,
    },
)
def molecule_image_extraction_molscribe(*args, **kwargs) -> dict[str, Any]:
    """
    Extracts the reaction scheme from the image using MolScribe model.
    https://github.com/thomas0809/MolScribe

    Args:
        image: bytes: Image in bytes.

    Returns:
        dict[str, Any]: Dictionary containing the reaction schema.
    """
    return _molscribe_extractor(*args, **kwargs)


@rxn_app.function(
    image=rxnscribe_image,
    gpu="A10G",
    volumes={
        "/root/.cache/huggingface": hf_cache_vol,
    },
)
def rxn_schema_extraction(*args, **kwargs) -> list[dict]:
    """
    Extracts the reaction scheme from the image using RxnScribe model.
    https://github.com/thomas0809/RxnScribe

    Args:
        image: bytes: Image in bytes.

    Returns:
        list[dict]: List of dictionaries containing the reaction schema.
    """
    return _rxnscribe_extractor(*args, **kwargs)


@rxn_app.function(image=rxn_mapper_image)
def get_rxn_mapping(rxns: list[str]) -> list[dict]:
    """
    Get the atom mapping of a reaction using RxnMapper.

    Args:
        rxns (list[str]): List of reaction strings

    Returns:
        list[dict]: List of atom mapping dictionaries

    Examples:
    >>> _get_rxn_mapper_confidence("CC(C)S.CN(C)C=O.Fc1cccnc1F.O=C([O-])[O-].[K+].[K+]>>CC(C)Sc1ncccc1F")
        [{'mapped_rxn': 'CN(C)C=O.F[c:5]1[n:6][cH:7][cH:8][cH:9]...',
        'confidence': 0.9565619900376546}]
    """
    return _get_rxn_mapper_confidence(rxns)


@rxn_app.function(image=rxn_utils_image)
def unify_rxn_smiles(*args, **kwargs):
    """
    Unify reaction SMILES strings.

    Args:
        rxn_smiles (str): The reaction SMILES string to unify.

    Returns:
        str: The unified reaction SMILES string.

    Examples:
    >>> _unify_rxn_smiles("CC.O.[Na+]~[Cl-]>>CCO")
        "CC.O.[Na+].[Cl-]>>CCO"
    """
    return _unify_rxn_smiles(*args, **kwargs)


@rxn_app.function(image=rxn_utils_image)
def canonicalize_rxn_smiles(*args, **kwargs):
    """
    Canonicalize reaction SMILES strings.

    Args:
        rxn_smiles (str): The reaction SMILES string to canonicalize.

    Returns:
        str: The canonicalized reaction SMILES string.

    Examples:
        >>> canonicalize_rxn_smiles("CO.O.C>>C(O) |f:1.2|")
            "CO.C.O>>CO |f:1.2|"
    """
    return _canonicalize_rxn_smiles(*args, **kwargs)
