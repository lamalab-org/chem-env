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
        "rxnscribe @ git+https://github.com/thomas0809/RxnScribe.git@ad6b1c75d40e563e68deca0491918885948d69c7",
        "molscribe @ git+https://github.com/thomas0809/MolScribe.git@main",
        "huggingface_hub",
        "loguru",
        "transformers",
        "torch",
    )
)


with rxnscribe_image.imports():
    import os
    from typing import Any
    import torch
    from rxnscribe import RxnScribe
    from huggingface_hub import hf_hub_download
    from molscribe import MolScribe


def _decimer_extractor(image_name: str) -> str:
    """
    Predicts the SMILES from the image using DECIMER models.
    https://github.com/Kohulan/DECIMER-Image_Transformer

    Args:
        image_path: str: Path to the image file.

    Returns:
        str: SMILES predicted.

    Example:
        >>> decimer_prediction("path/to/image.png")
        'CCO'
    """
    return decimer.predict_SMILES(f"/data/images/{image_name}")


def _molscribe_extractor(image_name: str) -> dict[str, Any]:
    """
    Predicts the molecule from the image using MolScribe model.

    Args:
        image_name: str: Name of the image file.

    Returns:
        dict: Dictionary containing the predicted molecule.
    """
    image_path = f"/data/images/{image_name}"
    ckpt_path = hf_hub_download("yujieq/MolScribe", "swin_base_char_aux_1m.pth")
    model = MolScribe(ckpt_path, device=torch.device("cpu"))
    output = model.predict_image_file(
        "assets/example.png", return_atoms_bonds=True, return_confidence=True
    )

    if "molfile" in output:
        output.pop("molfile")

    if os.path.exists(image_path):
        os.remove(image_path)

    return output


def _rxnscribe_extractor(image_name: str) -> list[dict]:
    """
    Predicts the reaction from the image using RxnScribe model.

    Args:
        image_name: str: Name of the image file.

    Returns:
        list: List of dictionaries containing the predicted reactions.
    """
    image_path = f"/data/images/{image_name}"
    ckpt_path = hf_hub_download("yujieq/RxnScribe", "pix2seq_reaction_full.ckpt")
    model = RxnScribe(ckpt_path, device=torch.device("cuda"))
    results = model.predict_image_file(image_path, molscribe=True, ocr=True)

    for result in results:
        for key, value in result.items():
            for v in value:
                if "molfile" in v:
                    v.pop("molfile")

    if os.path.exists(image_path):
        os.remove(image_path)

    return results
