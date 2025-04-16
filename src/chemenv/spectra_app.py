import os
from modal import App
from chemenv.tools.spectra_simulation import (
    spectra_simulation_image,
    SpectraAPI,
)

spectra_name = os.getenv("CHEMENV_NAME", "")
if spectra_name and not spectra_name.startswith("-"):
    spectra_name = f"-{spectra_name}"

# Create the app
spectra_app = App(f"chemenv_spectra{spectra_name}")


@spectra_app.function(image=spectra_simulation_image)
def simulate_spectra(*args, **kwargs):
    """
    Simulate 1H NMR, 13C NMR, and IR spectra for a given molecule using its SMILES string.
    If some of the spectra are not available, the function will return None for those spectra.

    Args:
        smiles (str): The SMILES string of the molecule.

    Returns:
        dict[str,str]: Simulated spectra data.

    Raises:
        ValueError: If the simulation fails.
    """
    return SpectraAPI.get_all_predictions(*args, **kwargs)
