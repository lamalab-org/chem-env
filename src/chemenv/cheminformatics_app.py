import os
from modal import App
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
    _get_molecular_formula,
)

cheminf_name = os.getenv("CHEMINF_NAME", "")
if cheminf_name and not cheminf_name.startswith("-"):
    cheminf_name = f"-{cheminf_name}"

# Create the app
cheminf_app = App(f"chemenv_cheminf{cheminf_name}")


@cheminf_app.function(image=rdkit_image)
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


@cheminf_app.function(image=rdkit_image)
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


@cheminf_app.function(image=rdkit_image)
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


@cheminf_app.function(image=rdkit_image)
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


@cheminf_app.function(image=rdkit_image)
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


@cheminf_app.function(image=rdkit_image)
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


@cheminf_app.function(image=rdkit_image)
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


@cheminf_app.function(image=rdkit_image)
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


@cheminf_app.function(image=rdkit_image)
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


@cheminf_app.function(image=rdkit_image)
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


@cheminf_app.function(image=rdkit_image)
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


@cheminf_app.function(image=rdkit_image)
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


@cheminf_app.function(image=rdkit_image)
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


@cheminf_app.function(image=rdkit_image)
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


@cheminf_app.function(image=rdkit_image)
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


@cheminf_app.function(image=rdkit_image)
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


@cheminf_app.function(image=rdkit_image)
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


@cheminf_app.function(image=rdkit_image)
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


@cheminf_app.function(image=rdkit_image)
def get_molecular_formula(*args, **kwargs):
    """
    Get the molecular formula of a molecule given its SMILES string.

    Args:
        smiles (str): The SMILES string of the molecule.

    Returns:
        str: The molecular formula of the molecule.

    Raises:
        ValueError: If the SMILES string is invalid.

    Example:
        >>> _get_molecular_formula("CCO")
            'C2H6O'
    """
    return _get_molecular_formula(*args, **kwargs)


@cheminf_app.function(image=mendeleev_image)
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
