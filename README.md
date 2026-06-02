# chem-env

`chem-env` is a small Python package for chemistry utilities that can be run locally or deployed as a Modal app.

It includes tools for cheminformatics, molecular conversion, reaction utilities, PubChem lookups, clinical trial queries, and spectra simulation.

## Install

```bash
pip install -e .
```

## Usage

Run the CLI:

```bash
chemenv --help
```

Deploy the Modal app:

```bash
chemenv deploy
```

Call a deployed function:

```python
import modal

app_name = "chemenv"
fxn = modal.Function.lookup(app_name, "get_tanimoto_similarity")
result = await fxn.remote.aio("CCO", "CC")
```

## Development

```bash
pip install -e ".[dev]"
pytest
```

## Citation

```bibtex
@article{ríos-garcía2026ai,
  title   = {AI scientists produce results without reasoning scientifically},
  author  = {Martiño Ríos-García and Nawaf Alampara and Chandan Gupta and Indrajeet Mandal and Sajid Mannan and Ali Asghar Aghajani and N. M. Anoop Krishnan and Kevin Maik Jablonka},
  year    = {2026},
  journal = {arXiv preprint arXiv: 2604.18805}
}
```
