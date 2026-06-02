from modal import Image


rxn_mapper_image = Image.debian_slim(python_version="3.9").pip_install(
    ["rdkit", "rxnmapper"]
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
        >>> _get_rxn_mapper_confidence("CC(C)S.CN(C)C=O.Fc1cccnc1F.O=C([O-])[O-].[K+].[K+]>>CC(C)Sc1ncccc1F")
            [{'mapped_rxn': 'CN(C)C=O.F[c:5]1[n:6][cH:7][cH:8][cH:9]...',
            'confidence': 0.9565619900376546}]
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
        >>> _unify_rxn_smiles("CC.O.[Na+]~[Cl-]>>CCO")
            "CC.O.[Na+].[Cl-]>>CCO"
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
        >>> canonicalize_rxn_smiles("CO.O.C>>C(O) |f:1.2|")
            "CO.C.O>>CO |f:1.2|"
    """
    return canonicalize_any(rxn_smiles)
