import os
from modal import App

from chemenv.cheminformatics_app import cheminf_app
from chemenv.clinical_trials_app import clinical_app
from chemenv.converters_app import converters_app
from chemenv.rxn_app import rxn_app
from chemenv.pubchem_app import pubchem_app
from chemenv.spectra_app import spectra_app
from chemenv.utils_app import utils_app


chemenv_name = os.getenv("CHEMENV_NAME", "")
if chemenv_name and not chemenv_name.startswith("-"):
    chemenv_name = f"-{chemenv_name}"

# Create the app
app = App(f"chemenv{chemenv_name}")
app.include(cheminf_app)
app.include(clinical_app)
app.include(converters_app)
app.include(rxn_app)
app.include(pubchem_app)
app.include(spectra_app)
app.include(utils_app)
