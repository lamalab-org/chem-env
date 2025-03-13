import os
from typing import Any
from modal import App, Volume
from chemenv.tools.rxn_schema_processing import (
    decimer_image,
    _decimer_extractor,
    rxnscribe_image,
    _molscribe_extractor,
    _rxnscribe_extractor,
)

MINUTES = 60  # 60 seconds
app_name = os.getenv("RXNENV_NAME", "")
if app_name and not app_name.startswith("-"):
    app_name = f"-{app_name}"


# Define the volume to safe checkpoints
hf_cache_vol = Volume.from_name("huggingface-cache", create_if_missing=True)


# Create the app
app = App(f"rxnenv{app_name}")


@app.function(image=decimer_image)
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


@app.function(
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


@app.function(
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
