import os
from modal import App, Volume
from chemenv.tools.clinical_trials import _clinical_image, ClinicalTrialsAPI


MINUTES = 60  # 60 seconds


# Define the volume to safe checkpoints
hf_cache_vol = Volume.from_name("huggingface-cache", create_if_missing=True)


clinical_name = os.getenv("CLINICAL_NAME", "")
if clinical_name and not clinical_name.startswith("-"):
    clinical_name = f"-{clinical_name}"


clinical_app = App(f"chemenv_clinical{clinical_name}")


@clinical_app.function(
    image=_clinical_image,
    gpu="A10G",
    volumes={
        "/root/.cache/huggingface": hf_cache_vol,
    },
    retries=3,
)
async def clinical_trials_search(*args, **kwargs) -> list[dict]:
    """
    Fetches clinical trial data and performs semantic search.

    This function can search in two modes:
    1. Search by query only (uses query term for both search and relevance ranking)
    2. Search by drug name with a query (finds all trials for a drug, then ranks by query relevance)

    Args:
        query (str): Text query to find relevant clinical trials
        drug_name (str, optional): Name of specific drug to search for. If provided,
                                  search will be drug-focused. Defaults to None.
        top_k (int, optional): Number of top results to return. Defaults to 10.

    Returns:
        list[dict]: List of dictionaries with relevant clinical trial data
    """
    clinical = await ClinicalTrialsAPI.create(*args, **kwargs)
    return await clinical.search_clinical_trials(*args, **kwargs)
