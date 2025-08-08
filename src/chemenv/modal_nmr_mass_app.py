# app.py
from modal import App, Image, fastapi_endpoint
import subprocess, json, textwrap
from pydantic import BaseModel
from rdkit import Chem
from typing import Optional

app = App("nmr-prediction-api")
p = App("nmr-predictor")

# Build an image with Node 20 + build tools, then install your JS deps.
image = (
    Image.debian_slim()
    .run_commands(
        # Install Node from NodeSource + build tools
        "apt-get update && apt-get install -y curl gnupg build-essential python3 make g++",
        "curl -fsSL https://deb.nodesource.com/setup_20.x | bash -",
        "apt-get install -y nodejs",
        # Prepare the tiny node project and install deps
        "mkdir -p /srv/js && cd /srv/js && "
        "npm init -y && npm pkg set type=module && "
        "npm i --omit=dev nmr-processing openchemlib isotopic-distribution",
    )
    .pip_install(["rdkit", "pydantic"])
)


def get_1d_spectrum_from_output(out, core: str = "1H"):
    return list(
        filter(lambda s: s.get("info", {}).get("nucleus") == core, out["spectra"])
    )


def get_2d_spectrum_from_output(out, pulse_sequence: str = "cosy"):
    return list(
        filter(
            lambda s: s.get("info", {}).get("pulseSequence") == pulse_sequence,
            out["spectra"],
        )
    )


class NMRPredictInput(BaseModel):
    smiles: str


class IsotopicDistributionInput(BaseModel):
    smiles: str
    ionization: Optional[str] = None


@app.function(image=image, timeout=60)
@fastapi_endpoint(method="POST")
def predict_isotopic_distribution(body: IsotopicDistributionInput):
    from rdkit.Chem import MolFromSmiles
    from rdkit.Chem.rdMolDescriptors import CalcMolFormula

    mol = MolFromSmiles(body.smiles)
    formula = CalcMolFormula(mol)
    print(f"Calculating isotopic distribution for formula: {formula}")
    js = textwrap.dedent(f"""
    import {{ IsotopicDistribution }} from "isotopic-distribution";


    const isotopicDistribution = new IsotopicDistribution("{formula}");
    console.log(JSON.stringify(isotopicDistribution.getPeaks()));
    """)
    out = subprocess.check_output(
        ["node", "--input-type=module", "-e", js],
        cwd="/srv/js",
        text=True,
        stderr=subprocess.STDOUT,
    )
    return json.loads(out)


@app.function(image=image, timeout=60)
@fastapi_endpoint(method="POST")
def predict_nmr(body: NMRPredictInput):
    smiles = body.smiles
    js = textwrap.dedent(f"""
        import {{ predictSpectra }} from "nmr-processing";
        import {{ Molecule }} from "openchemlib";

        const mol = Molecule.fromSmiles("{smiles}");
        const result = await predictSpectra(mol);
        console.log(JSON.stringify(result));
    """)
    out = subprocess.check_output(
        ["node", "--input-type=module", "-e", js],
        cwd="/srv/js",
        text=True,
        stderr=subprocess.STDOUT,
    )
    return json.loads(out)
