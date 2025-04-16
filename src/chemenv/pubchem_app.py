import os
from modal import App
from chemenv.tools.pubchem import (
    PubChem,
    pubchem_image as _pubchem_image,
)

pubchem_name = os.getenv("PUBCHEM_NAME", "")
if pubchem_name and not pubchem_name.startswith("-"):
    pubchem_name = f"-{pubchem_name}"

# Create the app
pubchem_app = App(f"chemenv_pubchem{pubchem_name}")


@pubchem_app.function(image=_pubchem_image)
async def get_pubchem_full_record(*args, **kwargs) -> dict:
    """
    Get the full PubChem record for a compound

    Args:
        compound (str): Any type of compound identifier (CID, SMILES, InChI, etc.)

    Returns:
        dict: Full PubChem record for the compound

    Raises:
        ValueError: If no record is found
    """
    pubchem = await PubChem.create(*args, **kwargs)
    return await pubchem._get_pubchem_full_record() if pubchem.cid else None


@pubchem_app.function(image=_pubchem_image)
async def get_number_atoms_pubchem(*args, **kwargs) -> int:
    """
    Get the number of atoms in a compound using RDKit.

    Args:
        compound (str): Any type of compound identifier (CID, SMILES, InChI, etc.)

    Returns:
        int: Number of atoms in the compound.

    Raises:
        ValueError: If no atoms are found

    Example:
        >>> await self._get_number_atoms()
            21
    """
    pubchem = await PubChem.create(*args, **kwargs)
    number_atoms = await pubchem._get_number_atoms()
    if number_atoms is None:
        raise ValueError("No number of atoms found")
    return number_atoms


@pubchem_app.function(image=_pubchem_image)
async def get_isomeric_smiles_pubchem(*args, **kwargs) -> str:
    """
    Get the isomeric SMILES for a compound from PubChem.

    Args:
        compound (str): Any type of compound identifier (CID, SMILES, InChI, etc.)

    Returns:
        str: Isomeric SMILES of the compound.

    Raises:
        ValueError: If the isomeric SMILES could not be retrieved

    Example:
        >>> await self._get_isomeric_smiles()
            'CCO'
    """
    pubchem = await PubChem.create(*args, **kwargs)
    isomeric_smiles = await pubchem._get_isomeric_smiles()
    if isomeric_smiles is None:
        raise ValueError("No isomeric SMILES found")
    return isomeric_smiles


@pubchem_app.function(image=_pubchem_image)
async def get_canonical_smiles_pubchem(*args, **kwargs) -> str:
    """
    Get the canonical SMILES for a compound from PubChem.

    Args:
        compound (str): Any type of compound identifier (CID, SMILES, InChI, etc.)

    Returns:
        str: Canonical SMILES of the compound.

    Raises:
        ValueError: If the canonical SMILES could not be retrieved

    Example:
        >>> await self._get_canonical_smiles()
            'CCO'
    """
    pubchem = await PubChem.create(*args, **kwargs)
    smiles = await pubchem._get_canonical_smiles()
    if smiles is None:
        raise ValueError("No canonical SMILES found")
    return smiles


@pubchem_app.function(image=_pubchem_image)
async def get_compound_mass_pubchem(*args, **kwargs) -> float:
    """
    Get the molecular weight of a compound from PubChem.

    Args:
        compound (str): Any type of compound identifier (CID, SMILES, InChI, etc.)

    Returns:
        float: Molecular weight of the compound.

    Raises:
        ValueError: If the molecular weight could not be retrieved

    Example:
        >>> await self._get_compound_mass()
            46.07
    """
    pubchem = await PubChem.create(*args, **kwargs)
    molecular_weight = await pubchem._get_compound_mass()
    if molecular_weight is None:
        raise ValueError("No molecular weight found")
    return molecular_weight


@pubchem_app.function(image=_pubchem_image)
async def get_compound_charge_pubchem(compound: str) -> int:
    """
    Get the charge of a compound from PubChem.

    Args:
        compound (str): Any type of compound identifier (CID, SMILES, InChI, etc.)

    Returns:
        int: Charge of the compound.

    Raises:
        ValueError: If the charge could not be retrieved

    Example:
        >>> await self._get_compound_charge()
            0
    """
    pubchem = await PubChem.create(compound)
    charge = await pubchem._get_compound_charge()
    if charge is None:
        raise ValueError("No charge found")
    return charge


@pubchem_app.function(image=_pubchem_image)
async def get_compound_formula_pubchem(compound: str) -> str:
    """
    Get the formula of a compound from PubChem.

    Args:
        compound (str): Any type of compound identifier (CID, SMILES, InChI, etc.)

    Returns:
        str: Formula of the compound.

    Raises:
        ValueError: If the formula could not be retrieved

    Example:
        >>> await self._get_compound_formula()
            'C2H6O'
    """
    pubchem = await PubChem.create(compound)
    formula = await pubchem._get_compound_formula()
    if formula is None:
        raise ValueError("No formula found")
    return formula


@pubchem_app.function(image=_pubchem_image)
async def get_number_isomers_pubchem(compound: str) -> int:
    """
    Get the number of compound isomers for a compound from PubChem.

    Args:
        compound (str): Any type of compound identifier (CID, SMILES, InChI, etc.)

    Returns:
        int: Number of compound isomers.

    Raises:
        ValueError: If the number of isomers could not be retrieved

    Example:
        >>> await self._get_number_isomers()
            2
    """
    pubchem = await PubChem.create(compound)
    number_isomers = await pubchem._get_number_isomers()
    if number_isomers is None:
        raise ValueError("No isomers found")
    return number_isomers


@pubchem_app.function(image=_pubchem_image, timeout=86399)
async def get_compound_isomers_pubchem(compound: str, limit: int = 10) -> list:
    """
    Get the compound isomers for a compound from PubChem.
    This function can take some time depending on the number of isomers.
    Returns a maximum of 10 isomers by default.

    Args:
        compound (str): Any type of compound identifier (CID, SMILES, InChI, etc.)
        limit (int, optional): Maximum number of isomers to return. Defaults to 10.

    Returns:
        list: List of compound isomers (limited to the specified number).

    Raises:
        ValueError: If the isomers could not be retrieved

    Example:
        >>> await self._get_compound_isomers()
            ['CCO', 'COC']
        >>> await self._get_compound_isomers(limit=5)
            ['CCO', 'COC', 'CC[O-]', 'C[O-]C', 'C[CH-]O']
    """
    pubchem = await PubChem.create(compound)
    data = await pubchem._get_compound_isomers(limit)
    if data is None:
        raise ValueError("No isomers found")
    return data


@pubchem_app.function(image=_pubchem_image, timeout=86399)
async def get_compound_isomers_pubchem_by_formula(
    formula: str, limit: int = 10
) -> list:
    """
    Get the compound isomers for a compound from PubChem.
    This function can take some time depending on the number of isomers.
    Returns a maximum of 10 isomers by default.

    Args:
        formulat(str): The empirical formula of the compound.
        limit (int, optional): Maximum number of isomers to return. Defaults to 10.

    Returns:
        list: List of compound isomers (limited to the specified number).

    Raises:
        ValueError: If the isomers could not be retrieved

    Example:
        >>> await self._get_compound_isomers()
            ['CCO', 'COC']
        >>> await self._get_compound_isomers(limit=5)
            ['CCO', 'COC', 'CC[O-]', 'C[O-]C', 'C[CH-]O']
    """
    data = await PubChem._get_compound_isomers_from_formula(formula, limit)
    if data is None:
        raise ValueError("No isomers found")
    return data


@pubchem_app.function(image=_pubchem_image)
async def get_number_heavy_atoms_pubchem(*args, **kwargs) -> int:
    """
    Get the number of heavy atoms in a compound from PubChem.

    Args:
        compound (str): Any type of compound identifier (CID, SMILES, InChI, etc.)

    Returns:
        int: Number of heavy atoms in the compound.

    Raises:
        ValueError: If the number of heavy atoms could not be retrieved

    Example:
        >>> await self._get_number_heavy_atoms()
            3
    """
    pubchem = await PubChem.create(*args, **kwargs)
    number_heavy_atoms = await pubchem._get_number_heavy_atoms()
    if number_heavy_atoms is None:
        raise ValueError("No heavy atoms found")
    return number_heavy_atoms


@pubchem_app.function(image=_pubchem_image)
async def get_number_chiral_atoms_pubchem(*args, **kwargs) -> int:
    """
    Get the number of chiral atoms in a compound from PubChem.

    Args:
        compound (str): Any type of compound identifier (CID, SMILES, InChI, etc.)

    Returns:
        int: Number of chiral atoms in the compound.

    Raises:
        ValueError: If the number of chiral atoms could not be retrieved

    Example:
        >>> await self._get_number_chiral_atoms()
            1
    """
    pubchem = await PubChem.create(*args, **kwargs)
    number_chiral_atoms = await pubchem._get_number_chiral_atoms()
    if number_chiral_atoms is None:
        raise ValueError("No chiral atoms found")
    return number_chiral_atoms


@pubchem_app.function(image=_pubchem_image)
async def get_c_nmr_spectra_pubchem(*args, **kwargs) -> dict:
    """
    Get the C-NMR spectra for a compound from PubChem.

    Args:
        compound (str): Any type of compound identifier (CID, SMILES, InChI, etc.)

    Returns:
        dict: C-NMR spectra data

    Raises:
        ValueError: If the data could not be retrieved

    Example:
        >>> await self._get_c_nmr_spectra()
            {'1': {'instrument': 'Bruker', 'frequency': '400 MHz',...}}
    """
    pubchem = await PubChem.create(*args, **kwargs)
    c_nmr_spectra = await pubchem._get_c_nmr_spectra()
    if c_nmr_spectra is None:
        raise ValueError("No C-NMR spectra found")
    return c_nmr_spectra


@pubchem_app.function(image=_pubchem_image)
async def get_h_nmr_spectra_pubchem(*args, **kwargs) -> dict:
    """
    Get the 1H-NMR spectra for a compound from PubChem.

    Args:
        compound (str): Any type of compound identifier (CID, SMILES, InChI, etc.)

    Returns:
        dict: 1H NMR spectra data

    Raises:
        ValueError: If the data could not be retrieved

    Example:
        >>> await self._get_h_nmr_spectra()
            {'1': {'instrument': 'Bruker', 'frequency': '400 MHz',...}}
    """
    pubchem = await PubChem.create(*args, **kwargs)
    h_nmr_spectra = await pubchem._get_h_nmr_spectra()
    if h_nmr_spectra is None:
        raise ValueError("No H-NMR spectra found")
    return h_nmr_spectra


@pubchem_app.function(image=_pubchem_image)
async def get_uv_spectra_pubchem(*args, **kwargs) -> dict:
    """
    Get the UV spectra for a compound from PubChem.

    Args:
        compound (str): Any type of compound identifier (CID, SMILES, InChI, etc.)

    Returns:
        str: UV spectra data

    Raises:
        ValueError: If the data could not be retrieved

    Example:
        >>> await self._get_uv_spectra()
            'Reference 1:\nMAX ABSORPTION: 210 nm\nReference 2:\nUV MAX: 210 nm'
    """
    pubchem = await PubChem.create(*args, **kwargs)
    uv_spectra = await pubchem._get_uv_spectra()
    if uv_spectra is None:
        raise ValueError("No UV spectra found")
    return uv_spectra


@pubchem_app.function(image=_pubchem_image)
async def get_ms_spectra_pubchem(*args, **kwargs) -> dict:
    """
    Get the MS spectra for a compound from PubChem.

    Args:
        compound (str): Any type of compound identifier (CID, SMILES, InChI, etc.)

    Returns:
        dict: MS spectra data

    Raises:
        ValueError: If the data could not be retrieved

    Example:
        >>> await self._get_ms_spectra()
            {'Mass bank ID': 'MoNA ID', ...}
    """
    pubchem = await PubChem.create(*args, **kwargs)
    ms_spectra = await pubchem._get_ms_spectra()
    if ms_spectra is None:
        raise ValueError("No MS spectra found")
    return ms_spectra


@pubchem_app.function(image=_pubchem_image)
async def get_ghs_classification_pubchem(*args, **kwargs) -> str:
    """
    Get the GHS classification for a compound from PubChem.

    Args:
        compound (str): Any type of compound identifier (CID, SMILES, InChI, etc.)

    Returns:
        dict: GHS classification data

    Raises:
        ValueError: If the data could not be retrieved

    Example:
        >>> await self._get_ghs_classification()
            {'H225': ['Highly flammable liquid and vapour'], ...}
    """
    pubchem = await PubChem.create(*args, **kwargs)
    ghs_classification = await pubchem._get_ghs_classification()
    if ghs_classification is None:
        raise ValueError("No GHS classification found")
    return ghs_classification


@pubchem_app.function(image=_pubchem_image)
async def get_patent_count(*args, **kwargs):
    """
    Get the number of patents for a compound from PubChem.

    Args:
        compound (str): Any type of compound identifier (CID, SMILES, InChI, etc.)

    Returns:
        int: Number of patents

    Raises:
        ValueError: If the data could not be retrieved

    Example:
        >>> await self._get_patent_count()
            2
    """
    pubchem = await PubChem.create(*args, **kwargs)
    patent_count = await pubchem._get_patent_count()
    if patent_count is None:
        raise ValueError("No patent count found")
    return patent_count


@pubchem_app.function(image=_pubchem_image)
async def return_physical_property(*args, **kwargs):
    """
    Get some physical properties of a compound from PubChem.

    Args:
        compound (str): Any type of compound identifier (CID, SMILES, InChI, etc.)

    Returns:
        dict: Physical property data

    Raises:
        ValueError: If the data could not be retrieved

    Example:
        >>> await self.return_physical_property()
            {'Experimental Properties': ['Appearance: clear colorless liquid', 'Boiling Point: 78.37...']}
    """
    pubchem = await PubChem.create(*args, **kwargs)
    physical_property = await pubchem._return_physical_property()
    if physical_property is None:
        raise ValueError("No physical property found")
    return physical_property
