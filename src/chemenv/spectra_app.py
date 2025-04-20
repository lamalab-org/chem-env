import os
import asyncio
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


@spectra_app.function(image=spectra_simulation_image, is_generator=True)
async def simulate_spectra_async(smiles: str, **kwargs):
    """
    Simulate 1H NMR, 13C NMR, and IR spectra for a given molecule using its SMILES string.
    If some of the spectra are not available, the function will return None for those spectra.

    Args:
        smiles (str): The SMILES string of the molecule.
        **kwargs: Additional arguments to pass to the prediction API.

    Returns:
        dict[str,str]: Simulated spectra data.

    Raises:
        ValueError: If the simulation fails.
    """
    return await SpectraAPI.get_all_predictions(smiles, **kwargs)


# Keep the synchronous version as well for backward compatibility
@spectra_app.function(image=spectra_simulation_image)
def simulate_spectra(*args, **kwargs):
    """
    Simulate 1H NMR, 13C NMR, and IR spectra for a given molecule using its SMILES string.
    This is a synchronous wrapper around the async version.

    Args:
        smiles (str): The SMILES string of the molecule.

    Returns:
        dict[str,str]: Simulated spectra data.

    Raises:
        ValueError: If the simulation fails.
    """
    return asyncio.run(SpectraAPI.get_all_predictions(*args, **kwargs))
