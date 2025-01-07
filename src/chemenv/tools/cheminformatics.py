from modal import Image

rdkit_image = (
    Image.debian_slim(python_version="3.12").pip_install("rdkit").pip_install("numpy")
)
mendeleev_image = Image.debian_slim().pip_install("mendeleev")

with rdkit_image.imports():
    from rdkit import Chem, DataStructs
    from rdkit.Chem import AllChem
    import numpy as np

with mendeleev_image.imports():
    from mendeleev import element


def get_tanimoto_similarity(s1: str, s2: str) -> float:
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
    try:
        mol1 = Chem.MolFromSmiles(s1)
        mol2 = Chem.MolFromSmiles(s2)
        fp1 = AllChem.GetMorganFingerprintAsBitVect(mol1, 2, nBits=2048)
        fp2 = AllChem.GetMorganFingerprintAsBitVect(mol2, 2, nBits=2048)
        return DataStructs.TanimotoSimilarity(fp1, fp2)
    except (TypeError, ValueError, AttributeError) as e:
        raise ValueError("Invalid SMILES strings") from e


def get_number_of_topologically_distinct_atoms(smiles: str, atomic_number: int = 1):
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
    try:
        molecule = Chem.MolFromSmiles(smiles)

        if atomic_number == 1:
            # add hydrogen
            mol = Chem.AddHs(molecule)
        else:
            mol = molecule

        # Get unique canonical atom rankings
        atom_ranks = list(Chem.rdmolfiles.CanonicalRankAtoms(mol, breakTies=False))

        # Select the unique element environments
        atom_ranks = np.array(atom_ranks)

        # Atom indices
        atom_indices = np.array(
            [
                atom.GetIdx()
                for atom in mol.GetAtoms()
                if atom.GetAtomicNum() == atomic_number
            ]
        )
        # Count them
        return len(set(atom_ranks[atom_indices]))
    except (TypeError, ValueError, AttributeError):
        return "Error: Not a valid SMILES string"


def get_element_info(identifier: str) -> dict:
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
    try:
        # Try to get the element
        if isinstance(identifier, int) or identifier.isdigit():
            el = element(int(identifier))
        else:
            el = element(identifier)

        # Collect basic information
        info = {
            "name": el.name,
            "symbol": el.symbol,
            "atomic_number": el.atomic_number,
            "mass": el.mass,
            "electron_configuration": str(el.ec.conf),
            "electronegativity": el.electronegativity,
            "group": el.group.name if el.group else None,
            "period": el.period,
            "block": el.block,
        }

        return info

    except ValueError:
        raise ValueError(f"Error: '{identifier}' is not a valid element identifier.")


def _get_number_atoms(smiles: str) -> int:
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
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError("Invalid SMILES string")
        return mol.GetNumAtoms()
    except Exception as e:
        raise ValueError(f"Invalid SMILES string: {e}")


def _get_number_heavy_atoms(smiles: str) -> int:
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
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError("Invalid SMILES string")
        return mol.GetNumHeavyAtoms()
    except Exception as e:
        raise ValueError(f"Invalid SMILES string: {e}")


def _get_canonical_smiles(smiles: str) -> str:
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
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError("Invalid SMILES string")
        return Chem.MolToSmiles(mol)
    except Exception as e:
        raise ValueError(f"Invalid SMILES string: {e}")


def _get_compound_charge(smiles: str) -> int:
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
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError("Invalid SMILES string")
        return Chem.GetFormalCharge(mol)
    except Exception as e:
        raise ValueError(f"Invalid SMILES string: {e}")


def _get_number_rings(smiles: str) -> int:
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
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError("Invalid SMILES string")
        rings = Chem.GetRingInfo(mol)
        return rings.numRings()
    except Exception as e:
        raise ValueError(f"Invalid SMILES string: {e}")


def _get_ring_sizes(smiles: str) -> list:
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
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError("Invalid SMILES string")
        rings = Chem.GetRingInfo(mol)
        return [len(ring) for ring in rings.AtomRings()]
    except Exception as e:
        raise ValueError(f"Invalid SMILES string: {e}")


def _get_number_aromatic_rings(smiles: str) -> int:
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
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError("Invalid SMILES string")
        rings = Chem.GetRingInfo(mol)
        return len(
            [
                ring
                for ring in rings.AtomRings()
                if all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring)
            ]
        )
    except Exception as e:
        raise ValueError(f"Invalid SMILES string: {e}")


def _get_aromatic_rings(smiles: str) -> list:
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
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError("Invalid SMILES string")
        rings = Chem.GetRingInfo(mol)
        return [
            ring
            for ring in rings.AtomRings()
            if all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring)
        ]
    except Exception as e:
        raise ValueError(f"Invalid SMILES string: {e}")


def _get_chiral_centers(smiles: str) -> int:
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
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError("Invalid SMILES string")
        Chem.AssignStereochemistry(mol)
        return Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    except Exception as e:
        raise ValueError(f"Invalid SMILES string: {e}")


def _get_number_chiral_centers(smiles: str) -> int:
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
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError("Invalid SMILES string")
        Chem.AssignStereochemistry(mol)
        return len(Chem.FindMolChiralCenters(mol, includeUnassigned=True))
    except Exception as e:
        raise ValueError(f"Invalid SMILES string: {e}")


def _get_number_cis_bonds(smiles: str) -> int:
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
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError("Invalid SMILES string")
        Chem.AssignStereochemistry(mol)
        return len(
            [
                bond
                for bond in mol.GetBonds()
                if bond.GetStereo() == Chem.BondStereo.STEREOCIS
            ]
        )
    except Exception as e:
        raise ValueError(f"Invalid SMILES string: {e}")


def _get_number_trans_bonds(smiles: str) -> int:
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
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError("Invalid SMILES string")
        Chem.AssignStereochemistry(mol)
        return len(
            [
                bond
                for bond in mol.GetBonds()
                if bond.GetStereo() == Chem.BondStereo.STEREOTRANS
            ]
        )
    except Exception as e:
        raise ValueError(f"Invalid SMILES string: {e}")


def _get_molecular_properties(smiles: str) -> dict:
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
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError("Invalid SMILES string")
        return {
            "logp": Chem.Descriptors.MolLogP(mol),
            "tpsa": Chem.Descriptors.TPSA(mol),
            "molecular_weight": Chem.Descriptors.ExactMolWt(mol),
            "rotatable_bonds": Chem.Descriptors.NumRotatableBonds(mol),
            "hbd": Chem.Descriptors.NumHDonors(mol),
            "hba": Chem.Descriptors.NumHAcceptors(mol),
        }
    except Exception as e:
        raise ValueError(f"Invalid SMILES string: {e}")


def _has_substructure(smiles: str, substructure_smarts: str) -> bool:
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
    try:
        mol = Chem.MolFromSmiles(smiles)
        pattern = Chem.MolFromSmarts(substructure_smarts)
        if mol is None or pattern is None:
            raise ValueError("Invalid SMILES or SMARTS pattern")
        return mol.HasSubstructMatch(pattern)
    except Exception as e:
        raise ValueError(f"Error in substructure matching: {e}")


def _get_substructure_count(smiles: str, substructure_smarts: str) -> int:
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
    try:
        mol = Chem.MolFromSmiles(smiles)
        pattern = Chem.MolFromSmarts(substructure_smarts)
        if mol is None or pattern is None:
            raise ValueError("Invalid SMILES or SMARTS pattern")
        return len(mol.GetSubstructMatches(pattern))
    except Exception as e:
        raise ValueError(f"Error in substructure matching: {e}")
