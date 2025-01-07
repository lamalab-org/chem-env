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
)
from chemenv.tools.pubchem import (
    PubChem,
    pubchem_image as _pubchem_image,
)
import os

chemenv_name = os.getenv("CHEMENV_NAME", "")
if chemenv_name and not chemenv_name.startswith("-"):
    chemenv_name = f"-{chemenv_name}"

# Create the app
app = App(f"chemenv{chemenv_name}")


@app.function(image=rdkit_image)
def get_tanimoto_similarity(*args, **kwargs):
    return _get_tanimoto_similarity(*args, **kwargs)


@app.function(image=rdkit_image)
def get_number_of_topologically_distinct_atoms(*args, **kwargs):
    return _get_topologically_distinct_atoms(*args, **kwargs)


@app.function(image=rdkit_image)
def get_number_atoms(*args, **kwargs):
    return _get_number_atoms(*args, **kwargs)


@app.function(image=rdkit_image)
def get_number_heavy_atoms(*args, **kwargs):
    return _get_number_heavy_atoms(*args, **kwargs)


@app.function(image=rdkit_image)
def get_canonical_smiles(*args, **kwargs):
    return _get_canonical_smiles(*args, **kwargs)


@app.function(image=rdkit_image)
def get_compound_charge(*args, **kwargs):
    return _get_compound_charge(*args, **kwargs)


@app.function(image=rdkit_image)
def get_number_rings(*args, **kwargs):
    return _get_number_rings(*args, **kwargs)


@app.function(image=rdkit_image)
def get_ring_sizes(*args, **kwargs):
    return _get_ring_sizes(*args, **kwargs)


@app.function(image=rdkit_image)
def get_number_aromatic_rings(*args, **kwargs):
    return _get_number_aromatic_rings(*args, **kwargs)


@app.function(image=rdkit_image)
def get_aromatic_rings(*args, **kwargs):
    return _get_aromatic_rings(*args, **kwargs)


@app.function(image=rdkit_image)
def get_chiral_centers(*args, **kwargs):
    return _get_chiral_centers(*args, **kwargs)


@app.function(image=rdkit_image)
def get_number_chiral_centers(*args, **kwargs):
    return _get_number_chiral_centers(*args, **kwargs)


@app.function(image=rdkit_image)
def get_number_cis_bonds(*args, **kwargs):
    return _get_number_cis_bonds(*args, **kwargs)


@app.function(image=rdkit_image)
def get_number_trans_bonds(*args, **kwargs):
    return _get_number_trans_bonds(*args, **kwargs)


@app.function(image=rdkit_image)
def get_molecular_properties(*args, **kwargs):
    return _get_molecular_properties(*args, **kwargs)


@app.function(image=rdkit_image)
def has_substructure(*args, **kwargs):
    return _has_substructure(*args, **kwargs)


@app.function(image=rdkit_image)
def get_substructure_count(*args, **kwargs):
    return _get_substructure_count(*args, **kwargs)


@app.function(image=mendeleev_image)
def get_element_info(*args, **kwargs):
    return _get_element_info(*args, **kwargs)


@app.function(image=_pubchem_image)
async def get_number_atoms_pubchem(compound_id: str) -> int:
    pubchem = await PubChem.create(compound_id)
    number_atoms = await pubchem._get_number_atoms()
    if number_atoms is None:
        raise ValueError("No number of atoms found")
    return number_atoms


@app.function(image=_pubchem_image)
async def get_isomeric_smiles_pubchem(compound_id: str) -> str:
    pubchem = await PubChem.create(compound_id)
    isomeric_smiles = await pubchem._get_isomeric_smiles()
    if isomeric_smiles is None:
        raise ValueError("No isomeric SMILES found")
    return isomeric_smiles


@app.function(image=_pubchem_image)
async def get_canonical_smiles_pubchem(compound_id: str) -> str:
    pubchem = await PubChem.create(compound_id)
    smiles = await pubchem._get_canonical_smiles()
    if smiles is None:
        raise ValueError("No canonical SMILES found")
    return smiles


@app.function(image=_pubchem_image)
async def get_compound_mass_pubchem(compound_id: str) -> float:
    pubchem = await PubChem.create(compound_id)
    molecular_weight = await pubchem._get_compound_mass()
    if molecular_weight is None:
        raise ValueError("No molecular weight found")
    return molecular_weight


@app.function(image=_pubchem_image)
async def get_compound_charge_pubchem(compound_id: str) -> int:
    pubchem = await PubChem.create(compound_id)
    charge = await pubchem._get_compound_charge()
    if charge is None:
        raise ValueError("No charge found")
    return charge


@app.function(image=_pubchem_image)
async def get_compound_formula_pubchem(compound_id: str) -> str:
    pubchem = await PubChem.create(compound_id)
    formula = await pubchem._get_compound_formula()
    if formula is None:
        raise ValueError("No formula found")
    return formula


@app.function(image=_pubchem_image)
async def get_number_isomers_pubchem(compound_id: str) -> int:
    pubchem = await PubChem.create(compound_id)
    number_isomers = await pubchem._get_number_isomers()
    if number_isomers is None:
        raise ValueError("No isomers found")
    return number_isomers


@app.function(image=_pubchem_image, timeout=86399)
async def get_compound_isomers_pubchem(*args, **kwargs):
    pubchem = await PubChem.create(*args, **kwargs)
    data = await pubchem._get_compound_isomers()
    if data is None:
        raise ValueError("No isomers found")
    return data


@app.function(image=_pubchem_image)
async def get_number_heavy_atoms_pubchem(compound_id: str) -> int:
    pubchem = await PubChem.create(compound_id)
    number_heavy_atoms = await pubchem._get_number_heavy_atoms()
    if number_heavy_atoms is None:
        raise ValueError("No heavy atoms found")
    return number_heavy_atoms


@app.function(image=_pubchem_image)
async def get_number_chiral_atoms_pubchem(compound_id: str) -> int:
    pubchem = await PubChem.create(compound_id)
    number_chiral_atoms = await pubchem._get_number_chiral_atoms()
    if number_chiral_atoms is None:
        raise ValueError("No chiral atoms found")
    return number_chiral_atoms


@app.function(image=_pubchem_image)
async def get_c_nmr_spectra_pubchem(compound_id: str) -> dict:
    pubchem = await PubChem.create(compound_id)
    c_nmr_spectra = await pubchem._get_c_nmr_spectra()
    if c_nmr_spectra is None:
        raise ValueError("No C-NMR spectra found")
    return c_nmr_spectra


@app.function(image=_pubchem_image)
async def get_h_nmr_spectra_pubchem(compound_id: str) -> dict:
    pubchem = await PubChem.create(compound_id)
    h_nmr_spectra = await pubchem._get_h_nmr_spectra()
    if h_nmr_spectra is None:
        raise ValueError("No H-NMR spectra found")
    return h_nmr_spectra


@app.function(image=_pubchem_image)
async def get_uv_spectra_pubchem(compound_id: str) -> dict:
    pubchem = await PubChem.create(compound_id)
    uv_spectra = await pubchem._get_uv_spectra()
    if uv_spectra is None:
        raise ValueError("No UV spectra found")
    return uv_spectra


@app.function(image=_pubchem_image)
async def get_ms_spectra_pubchem(compound_id: str) -> dict:
    pubchem = await PubChem.create(compound_id)
    ms_spectra = await pubchem._get_ms_spectra()
    if ms_spectra is None:
        raise ValueError("No MS spectra found")
    return ms_spectra


@app.function(image=_pubchem_image)
async def get_ghs_classification_pubchem(compound_id: str) -> str:
    pubchem = await PubChem.create(compound_id)
    ghs_classification = await pubchem._get_ghs_classification()
    if ghs_classification is None:
        raise ValueError("No GHS classification found")
    return ghs_classification
