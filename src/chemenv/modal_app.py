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
from chemenv.tools.cheminformatics import (
    rdkit_image,
    mendeleev_image,
    get_tanimoto_similarity as _get_tanimoto_similarity,
    get_number_of_topologically_distinct_atoms as _get_topologically_distinct_atoms,
    get_element_info as _get_element_info,
    _get_number_atoms,
    _get_number_heavy_atoms,
    _get_canonical_smiles,
    _get_compound_charge,
    _get_number_rings,
    _get_number_aromatic_rings,
    _get_aromatic_rings,
    _get_ring_sizes,
    _get_chiral_centers,
    _get_number_chiral_centers,
    _get_number_cis_bonds,
    _get_number_trans_bonds,
    _get_molecular_properties,
    _has_substructure,
    _get_substructure_count,
    _pka_from_smiles,
)
from chemenv.tools.util_tool import (
    exomol_image,
    _get_functional_groups,
    is_patented_image,
    _is_patented,
    _is_buyable,
)
from chemenv.tools.pubchem import (
    PubChem,
    pubchem_image as _pubchem_image,
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
from chemenv.tools.rxn_utils import (
    rxn_mapper_image,
    _get_rxn_mapper_confidence,
    rxn_utils_image,
    _unify_rxn_smiles,
    _canonicalize_rxn_smiles,
)

MINUTES = 60  # 60 seconds

chemenv_name = os.getenv("CHEMENV_NAME", "")
if chemenv_name and not chemenv_name.startswith("-"):
    chemenv_name = f"-{chemenv_name}"

app_name = os.getenv("RXNENV_NAME", "")
if app_name and not app_name.startswith("-"):
    app_name = f"-{app_name}"

# Define the volume to safe checkpoints
hf_cache_vol = Volume.from_name("huggingface-cache", create_if_missing=True)

# Create the app
rxn_app = App(f"rxnenv{app_name}")

# Create the app
app = App(f"chemenv{chemenv_name}")


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


@app.function(image=rdkit_image)
def get_tanimoto_similarity(*args, **kwargs):
    """
    Calculate the Tanimoto similarity of two SMILES strings.

    Args:
        s1 (str): The first SMILES string.
        s2 (str): The second SMILES string.

    Returns:
        float: The Tanimoto similarity between the two SMILES strings.

    Raises:
        ValueError: If either of the input SMILES strings is invalid.

    This function calculates the Tanimoto similarity between two SMILES strings.
    It uses the RDKit library to calculate the Morgan fingerprints of the molecules.
    The Tanimoto similarity is then calculated using the DataStructs module of the RDKit library.

    Example:
        >>> get_tanimoto_similarity("CCO", "CC")
            0.143
    """
    return _get_tanimoto_similarity(*args, **kwargs)


@app.function(image=rdkit_image)
def get_number_of_topologically_distinct_atoms(*args, **kwargs):
    """Return the number of unique `element` environments based on environmental topology.
    This corresponds to the number of peaks one could maximally observe in an NMR spectrum.
    Args:
        smiles (str): SMILES string
        atomic_number (int, optional): Atomic number. Defaults to 1.

    Returns:
        int: Number of unique environments.

    Raises:
        ValueError: If not a valid SMILES string

    Example:
        >>> get_number_of_topologically_distinct_atoms("CCO", 1)
            3

        >>> get_number_of_topologically_distinct_atoms("CCO", 6)
            2
    """
    return _get_topologically_distinct_atoms(*args, **kwargs)


@app.function(image=rdkit_image)
def get_number_atoms(*args, **kwargs):
    """
    Get the number of atoms in a molecule given its SMILES string.

    Args:
        smiles (str): The SMILES string of the molecule.

    Returns:
        int: The number of atoms in the molecule.

    Raises:
        ValueError: If the SMILES string is invalid.

    Example:
        >>> _get_number_atoms("CCO")
            3
    """
    return _get_number_atoms(*args, **kwargs)


@app.function(image=rdkit_image)
def get_number_heavy_atoms(*args, **kwargs):
    """
    Get the number of heavy atoms in a molecule given its SMILES string.

    Args:
        smiles (str): The SMILES string of the molecule.

    Returns:
        int: The number of heavy atoms in the molecule.

    Raises:
        ValueError: If the SMILES string is invalid.

    Example:
        >>> _get_number_heavy_atoms("CCO")
            3
    """
    return _get_number_heavy_atoms(*args, **kwargs)


@app.function(image=rdkit_image)
def get_canonical_smiles(*args, **kwargs):
    """
    Get the canonical SMILES string of a molecule given its SMILES string.

    Args:
        smiles (str): The SMILES string of the molecule.

    Returns:
        str: The canonical SMILES string of the molecule.

    Raises:
        ValueError: If the SMILES string is invalid.

    Example:
        >>> _canonical_smiles("CCO")
            'CCO'
    """
    return _get_canonical_smiles(*args, **kwargs)


@app.function(image=rdkit_image)
def get_compound_charge(*args, **kwargs):
    """
    Get the charge of a molecule given its SMILES string.

    Args:
        smiles (str): The SMILES string of the molecule.

    Returns:
        int: The charge of the molecule.

    Raises:
        ValueError: If the SMILES string is invalid.

    Example:
        >>> _get_compound_charge("CCO")
            0
    """
    return _get_compound_charge(*args, **kwargs)


@app.function(image=rdkit_image)
def get_number_rings(*args, **kwargs):
    """
    Get the number of aromatic rings in a molecule given its SMILES string.

    Args:
        smiles (str): The SMILES string of the molecule.

    Returns:
        int: The number of aromatic rings in the molecule.

    Raises:
        ValueError: If the SMILES string is invalid.

    Example:
        >>> _get_number_aromatic_rings("c1ccccc1")
            1
    """
    return _get_number_rings(*args, **kwargs)


@app.function(image=rdkit_image)
def get_ring_sizes(*args, **kwargs):
    """
    Get the sizes of rings in a molecule given its SMILES string.

    Args:
        smiles (str): The SMILES string of the molecule.

    Returns:
        list: A list of ring sizes in the molecule.

    Raises:
        ValueError: If the SMILES string is invalid.

    Example:
        >>> _get_ring_sizes("C1CCCCC1")
            [6]
    """
    return _get_ring_sizes(*args, **kwargs)


@app.function(image=rdkit_image)
def get_number_aromatic_rings(*args, **kwargs):
    """
    Get the number of aromatic rings in a molecule given its SMILES string.

    Args:
        smiles (str): The SMILES string of the molecule.

    Returns:
        int: The number of aromatic rings in the molecule.

    Raises:
        ValueError: If the SMILES string is invalid.

    Example:
        >>> _get_number_aromatic_rings("c1ccccc1")
            1
    """
    return _get_number_aromatic_rings(*args, **kwargs)


@app.function(image=rdkit_image)
def get_aromatic_rings(*args, **kwargs):
    """
    Get the aromatic rings in a molecule given its SMILES string.

    Args:
        smiles (str): The SMILES string of the molecule.

    Returns:
        list: A list of aromatic rings in the molecule.

    Raises:
        ValueError: If the SMILES string is invalid.

    Example:
        >>> _get_aromatic_rings("c1ccccc1")
            [[0, 1, 2, 3, 4, 5]]
    """
    return _get_aromatic_rings(*args, **kwargs)


@app.function(image=rdkit_image)
def get_chiral_centers(*args, **kwargs):
    """
    Get the number of chiral centers in a molecule given its SMILES string.

    Args:
        smiles (str): The SMILES string of the molecule.

    Returns:
        int: The number of chiral centers in the molecule.

    Raises:
        ValueError: If the SMILES string is invalid.

    Example:
        >>> _get_chiral_centers("CCO")
            0
    """
    return _get_chiral_centers(*args, **kwargs)


@app.function(image=rdkit_image)
def get_number_chiral_centers(*args, **kwargs):
    """
    Get the number of chiral centers in a molecule given its SMILES string.

    Args:
        smiles (str): The SMILES string of the molecule.

    Returns:
        int: The number of chiral centers in the molecule.

    Raises:
        ValueError: If the SMILES string is invalid.

    Example:
        >>> _get_number_chiral_centers("CCO")
            0
    """
    return _get_number_chiral_centers(*args, **kwargs)


@app.function(image=rdkit_image)
def get_number_cis_bonds(*args, **kwargs):
    """
    Get the number of cis bonds in a molecule given its SMILES string.

    Args:
        smiles (str): The SMILES string of the molecule.

    Returns:
        int: The number of cis bonds in the molecule.

    Raises:
        ValueError: If the SMILES string is invalid.

    Example:
        >>> _get_number_cis_bonds("C/C=C/C")
            1
    """
    return _get_number_cis_bonds(*args, **kwargs)


@app.function(image=rdkit_image)
def get_number_trans_bonds(*args, **kwargs):
    """
    Get the number of trans bonds in a molecule given its SMILES string.

    Args:
        smiles (str): The SMILES string of the molecule.

    Returns:
        int: The number of trans bonds in the molecule.

    Raises:
        ValueError: If the SMILES string is invalid.

    Example:
        >>> _get_number_trans_bonds("C/C=C/C")
            1
    """
    return _get_number_trans_bonds(*args, **kwargs)


@app.function(image=rdkit_image)
def get_molecular_properties(*args, **kwargs):
    """Get basic molecular properties.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        dict: Dictionary containing properties like LogP, TPSA, etc.

    Raises:
        ValueError: If the SMILES string is invalid

    Example:
        >>> _get_molecular_properties("CCO")
            {'logp': 0.22399999999999998, 'tpsa': 20.23, 'molecular_weight': 46.069, 'rotatable_bonds': 1, 'hbd': 1, 'hba': 1}
    """
    return _get_molecular_properties(*args, **kwargs)


@app.function(image=rdkit_image)
def has_substructure(*args, **kwargs):
    """Check if a molecule contains a specific substructure.

    Args:
        smiles (str): SMILES string of the molecule
        substructure_smarts (str): SMARTS pattern of the substructure

    Returns:
        bool: True if substructure is present

    Raises:
        ValueError: If the SMILES or SMARTS pattern is invalid

    Example:
        >>> _has_substructure("CCO", "CO")
            True
    """
    return _has_substructure(*args, **kwargs)


@app.function(image=rdkit_image)
def get_substructure_count(*args, **kwargs):
    """Get the count of a specific substructure in a molecule.

    Args:
        smiles (str): SMILES string of the molecule
        substructure_smarts (str): SMARTS pattern of the substructure

    Returns:
        int: Number of occurrences of the substructure

    Raises:
        ValueError: If the SMILES string or SMARTS pattern is invalid

    Example:
        >>> _get_substructure_count("CCO", "C")
            2
    """
    return _get_substructure_count(*args, **kwargs)


@app.function(image=rdkit_image)
def pka_from_smiles(*args, **kwargs):
    """
    Calculate the pKa of a molecule given its SMILES string.

    Args:
        smiles (str): The SMILES string of the molecule.

    Returns:
        float: The pKa of the molecule.

    Raises:
        ValueError: If the SMILES string is invalid.

    Example:
        >>> _pka_from_smiles("CCO")
            16.0
    """
    return _pka_from_smiles(*args, **kwargs)


@app.function(image=mendeleev_image)
def get_element_info(*args, **kwargs):
    """
    A function to retrieve basic information about a chemical element based on its identifier.

    Args:
        identifier (str): The identifier of the chemical element.

    Returns:
        dict: A dictionary containing basic information about the chemical element including its name,
            symbol, atomic number, mass, electron configuration, electronegativity, group, period, and block.

    Raises:
        ValueError: If the identifier is not a valid element identifier.

    Example:
        >>> get_element_info("H")["name"]
            'Hydrogen'
    """
    return _get_element_info(*args, **kwargs)


@app.function(image=exomol_image)
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


@app.function(image=is_patented_image)
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


@app.function(image=is_patented_image)
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


@app.function(image=_pubchem_image)
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


@app.function(image=_pubchem_image)
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


@app.function(image=_pubchem_image)
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


@app.function(image=_pubchem_image)
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


@app.function(image=_pubchem_image)
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


@app.function(image=_pubchem_image)
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


@app.function(image=_pubchem_image)
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


@app.function(image=_pubchem_image)
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


@app.function(image=_pubchem_image, timeout=86399)
async def get_compound_isomers_pubchem(*args, **kwargs):
    """
    Get the compound isomers for a compound from PubChem.
    This function can take some time depending on the number of isomers.

    Args:
        compound (str): Any type of compound identifier (CID, SMILES, InChI, etc.)

    Returns:
        list: List of compound isomers.

    Raises:
        ValueError: If the isomers could not be retrieved

    Example:
        >>> await self._get_compound_isomers()
            ['CCO', 'COC']
    """
    pubchem = await PubChem.create(*args, **kwargs)
    data = await pubchem._get_compound_isomers()
    if data is None:
        raise ValueError("No isomers found")
    return data


@app.function(image=_pubchem_image)
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


@app.function(image=_pubchem_image)
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


@app.function(image=_pubchem_image)
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


@app.function(image=_pubchem_image)
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


@app.function(image=_pubchem_image)
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


@app.function(image=_pubchem_image)
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


@app.function(image=_pubchem_image)
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


@app.function(image=_pubchem_image)
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


@app.function(image=_pubchem_image)
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


@app.function(image=_converters_image)
async def get_iupac_name(smiles: str, timeout: int = 10) -> str:
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
    converter = _Smiles2Name(smiles, timeout)
    try:
        name = await converter.get_name()
        return name
    except Exception as e:
        raise ValueError(f"Error: {e}") from e


@app.function(image=_converters_image)
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


@app.function(image=_converters_image)
def convert_to_selfies(smiles: str) -> str:
    """
    Convert SMILES to SELFIES encoding.

    Args:
        smiles (str): The SMILES string to convert.

    Returns:
        str: The SELFIES encoding of the molecule.
    """
    return _smiles_to_selfies(smiles)


@app.function(image=_converters_image)
def convert_to_deepsmiles(smiles: str) -> str:
    """
    Convert SMILES to DeepSMILES encoding.

    Args:
        smiles (str): The SMILES string to convert.

    Returns:
        str: The DeepSMILES encoding of the molecule.
    """
    return _smiles_to_deepsmiles(smiles)


@app.function(image=_converters_image)
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


@app.function(image=_converters_image)
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


@app.function(image=_safe_image)
def convert_to_safe(smiles: str) -> str:
    """
    Convert SMILES to SAFE encoding.

    Args:
        smiles (str): The SMILES string to convert.

    Returns:
        str: The SAFE encoding of the molecule.
    """
    return _smiles_to_safe(smiles)


@app.function(image=_converters_image)
def selfies_to_smiles(selfies: str) -> str:
    """
    Convert SELFIES to SMILES.

    Args:
        selfies (str): The SELFIES encoding to convert.

    Returns:
        str: The SMILES string of the molecule.
    """
    return _selfies_to_smiles(selfies)


@app.function(image=_converters_image)
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


@app.function(image=_converters_image)
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


@app.function(image=rxn_mapper_image)
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


@app.function(image=rxn_utils_image)
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


@app.function(image=rxn_utils_image)
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
