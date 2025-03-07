from modal import Image

is_patented_image = Image.debian_slim(python_version=3.12).pip_install(
    ["molbloom", "loguru"]
)

with is_patented_image.imports():
    import molbloom
    from loguru import logger


def _is_patented(smiles: str) -> bool:
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
    logger.debug(f"Checking if {smiles} is patented")
    try:
        r = molbloom.buy(smiles, canonicalize=True, catalog="surechembl")
    except Exception as e:
        raise ValueError(f"Error while checking if {smiles} is patented: {e}")
    if r:
        return True
    else:
        return False


def _is_buyable(smiles: str) -> bool:
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
    logger.debug(f"Checking if {smiles} is buyable")
    try:
        r = molbloom.buy(smiles, canonicalize=True)
    except Exception as e:
        raise ValueError(f"Error while checking if {smiles} is buyable: {e}")
    if r:
        return True
    else:
        return False
