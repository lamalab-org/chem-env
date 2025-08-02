import os
from modal import App
from chemenv.tools.util_tool import (
    exomol_image,
    _get_functional_groups,
    is_patented_image,
    _is_patented,
    _is_buyable,
)

utils_name = os.getenv("UTILS_NAME", "")
if utils_name and not utils_name.startswith("-"):
    utils_name = f"-{utils_name}"

# Create the app
utils_app = App(f"chemenv_utils{utils_name}")


@utils_app.function(image=exomol_image)
def get_functional_groups(*args, **kwargs):
    """
    Get the functional groups of a molecule using Exomol

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        list[str]: List of functional groups of the molecule

    Raises:
        ValueError: If an error occurs while getting the functional groups

    Example:
        >>> _get_functional_groups("CCO")
            ['Alcohol']
    """
    return _get_functional_groups(*args, **kwargs)


@utils_app.function(image=is_patented_image)
def is_patented(*args, **kwargs):
    """
    Check if a molecule is patented using Molbloom

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        str: "Patented" if the molecule is patented, "Novel" otherwise

    Raises:
        ValueError: If an error occurs while checking if the molecule is patented

    Examples:
        >>> _is_patented("CCO")
            False
    """
    return _is_patented(*args, **kwargs)


@utils_app.function(image=is_patented_image)
def is_buyable(*args, **kwargs):
    """
    Check if a molecule is buyable using Molbloom

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        str: "Buyable" if the molecule is buyable, "Not buyable" otherwise

    Raises:
        ValueError: If an error occurs while checking if the molecule is buyable

    Examples:
        >>> _is_buyable("CCO")
            True
    """
    return _is_buyable(*args, **kwargs)
