from modal import Image


rxn_mapper_image = Image.debian_slim(python_version="3.6").pip_install(
    ["rdkit", "rxn_mapper"]
)

with rxn_mapper_image.imports():
    from rxnmapper import RXNMapper

rxn_utils_image = Image.debian_slim(python_version="3.12").pip_install(
    ["rxn-chem-utils"]
)

with rxn_utils_image.imports():
    from rxn.chemutils.reaction_smiles import (
        ReactionFormat,
        to_reaction_smiles,
        parse_any_reaction_smiles,
    )
    from rxn.chemutils.miscellaneous import canonicalize_any


def _get_rxn_mapper_confidence(rxns: list[str]) -> list[dict]:
    """
    Get the confidence of a reaction using RxnMapper

    Args:
        rxn (str): Reaction string

    Returns:
        float: Confidence of the reaction

    Examples:
        >>> _get_rxn_mapper_confidence("CCO>>CC")
            0.95
    """
    rxn_mapper = RXNMapper()
    return rxn_mapper.get_attention_guided_atom_maps(rxns)


def _unify_rxn_smiles(rxn_smiles: str) -> str:
    """
    Unify the SMILES of a reaction using RxnMapper

    Args:
        rxn_smiles (str): Reaction SMILES

    Returns:
        str: Unified SMILES of the reaction

    Examples:
        >>> _unify_rxn_smiles("CCO>>CC")
            "CCO>>CC"
    """
    return to_reaction_smiles(
        parse_any_reaction_smiles(rxn_smiles), ReactionFormat.STANDARD
    )


def _canonicalize_rxn_smiles(rxn_smiles: str) -> str:
    """
    Canonicalize the SMILES of a reaction using RxnMapper

    Args:
        rxn_smiles (str): Reaction SMILES

    Returns:
        str: Canonicalized SMILES of the reaction

    Examples:
        >>> canonicalize_rxn_smiles("CCO>>CC")
            "CCO>>CC"
    """
    return canonicalize_any(rxn_smiles)
