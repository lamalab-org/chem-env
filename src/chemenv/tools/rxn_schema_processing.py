from modal import Image


decimer_image = Image.debian_slim(python_version="3.10.0").pip_install("decimer")


with decimer_image.imports():
    import DECIMER as decimer


rxnscribe_image = (
    Image.debian_slim(python_version="3.11")
    .apt_install(
        "git",
        "libgl1-mesa-glx",
        "libglib2.0-0",
        "libsm6",
        "libxext6",
        "libxrender-dev",
    )
    .pip_install(
        "opencv-python",
        "rxnscribe @ git+https://github.com/lamalab-org/RxnScribe.git@main",
        "molscribe @ git+https://github.com/lamalab-org/MolScribe.git@main",
        "huggingface_hub",
        "loguru",
        "transformers",
        "torch",
    )
)


with rxnscribe_image.imports():
    from typing import Any
    import torch
    from rxnscribe import RxnScribe
    from huggingface_hub import hf_hub_download
    from molscribe import MolScribe


def _decimer_extractor(image: bytes) -> str:
    """
    Predicts the SMILES from the image using DECIMER models.
    https://github.com/Kohulan/DECIMER-Image_Transformer

    Args:
        image: bytes: Image file containing the molecule.

    Returns:
        str: SMILES predicted.

    Example:
        >>> decimer_prediction("path/to/image.png")
        'CCO'
    """
    return decimer.predict_SMILES(image)


def _molscribe_extractor(image: bytes) -> dict[str, Any]:
    """
    Predicts the molecule from the image using MolScribe model.

    Args:
        image: bytes: Image file with the molecule.

    Returns:
        dict: Dictionary containing the predicted molecule.
    """
    ckpt_path = hf_hub_download("yujieq/MolScribe", "swin_base_char_aux_1m.pth")
    model = MolScribe(ckpt_path, device=torch.device("cuda"))
    output = model.predict_image_file(
        image, return_atoms_bonds=True, return_confidence=True
    )

    if "molfile" in output:
        output.pop("molfile")

    return output


def _rxnscribe_extractor(image: bytes) -> list[dict]:
    """
    Predicts the reaction from the image using RxnScribe model.

    Args:
        image: bytes: Image file with the reaction schema.

    Returns:
        list: List of dictionaries containing the predicted reactions.
    """
    ckpt_path = hf_hub_download("yujieq/RxnScribe", "pix2seq_reaction_full.ckpt")
    model = RxnScribe(ckpt_path, device=torch.device("cuda"))
    results = model.predict_image_file(image, molscribe=True, ocr=True)

    for result in results:
        for key, value in result.items():
            for v in value:
                if "molfile" in v:
                    v.pop("molfile")

    return results
