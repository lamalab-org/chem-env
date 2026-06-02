from modal import Image
from typing import Any


_clinical_image = Image.debian_slim(python_version="3.12").pip_install(
    [
        "asyncio",
        "aiohttp",
        "numpy",
        "scikit-learn",
        "sentence-transformers",
        "tenacity",
        "loguru",
    ]
)


with _clinical_image.imports():
    from sklearn.metrics.pairwise import cosine_similarity
    from sentence_transformers import SentenceTransformer
    import aiohttp
    import asyncio
    import numpy as np
    import time
    import gc
    from typing import Any
    from loguru import logger


class ClinicalTrialsAPI:
    """A class for interacting with the ClinicalTrials.gov API and performing semantic search on trial data."""

    def __init__(self, model_name: str = "intfloat/multilingual-e5-large-instruct"):
        """Initialize the ClinicalTrialsAPI.

        Args:
            model_name (str): Name of the sentence transformer model to use for embeddings.
        """
        self.model_name = model_name
        self.base_url = "https://clinicaltrials.gov/api/v2/studies"

    @classmethod
    async def create(cls, model_name: str = "intfloat/multilingual-e5-large-instruct"):
        """Factory method to create a ClinicalTrialsAPI instance.

        Args:
            model_name (str, optional): Name of the sentence transformer model to use for embeddings. Defaults to 'intfloat/multilingual-e5-large-instruct'.

        Returns:
            ClinicalTrialsAPI: A configured instance of the API client
        """
        instance = cls(model_name)
        return instance

    def get_model(self):
        """Load the sentence transformer model on demand."""
        return SentenceTransformer(self.model_name)

    async def make_api_request(
        self,
        url: str,
        method: str,
        params: dict[str, Any] | None = None,
        verbose: bool = False,
    ) -> dict[str, Any]:
        """Make an asynchronous API request with retry logic for robustness.

        Args:
            url (str): The API endpoint URL
            method (str): The HTTP method (only "GET" is supported)
            params (dict[str, Any], optional): Query parameters for the request. Defaults to None.
            verbose (bool, optional): Whether to log detailed information. Defaults to False.

        Returns:
            Dict[str, Any]: The parsed JSON response

        Raises:
            ValueError: If method is not "GET"
            aiohttp.ClientError: If all retry attempts fail
        """
        if method.upper() != "GET":
            raise ValueError("Only GET method is supported in this implementation")

        if verbose:
            logger.info(f"Making {method} request to {url} with params: {params}")

        timeout = aiohttp.ClientTimeout(total=30, connect=10)

        headers = {
            "User-Agent": "ChemEnv-Research-Tool/1.0",
            "Accept": "application/json",
        }

        async with aiohttp.ClientSession(timeout=timeout) as session:
            try:
                start_time = time.time()

                async with session.get(url, params=params, headers=headers) as response:
                    elapsed_time = time.time() - start_time

                    if verbose:
                        logger.info(
                            f"Request completed in {elapsed_time:.2f}s with status code: {response.status}"
                        )

                    response.raise_for_status()

                    if response.status == 429:
                        retry_after = int(response.headers.get("Retry-After", 60))
                        logger.warning(
                            f"Rate limited. Retrying after {retry_after} seconds."
                        )
                        await asyncio.sleep(retry_after)
                        raise aiohttp.ClientResponseError(
                            request_info=response.request_info,
                            history=response.history,
                            status=429,
                            message="Rate limited",
                        )

                    data = await response.json()
                    return data

            except aiohttp.ContentTypeError as e:
                response_text = await response.text()
                logger.error(f"Failed to decode JSON response: {str(e)}")
                logger.error(f"Response text: {response_text[:200]}...")
                raise aiohttp.ClientError(f"JSON decode error: {str(e)}")
            except Exception as e:
                logger.error(f"Request failed: {str(e)}")
                raise

    async def fetch_all_studies(self, drug_name: str) -> list[dict[str, Any]]:
        """Fetch all studies related to a specific drug from ClinicalTrials.gov

        Args:
            drug_name (str): The name of the drug to search for

        Returns:
            list[dict[str, Any]]: A list of dictionaries containing study data
        """
        params = {"query.term": drug_name, "pageSize": 100, "format": "json"}
        all_studies = []
        next_page_token = None

        try:
            while True:
                if next_page_token:
                    params["pageToken"] = next_page_token

                data = await self.make_api_request(
                    url=self.base_url, method="GET", params=params, verbose=True
                )

                studies = data.get("studies", [])
                all_studies.extend(studies)

                next_page_token = data.get("nextPageToken")
                if not next_page_token:
                    break

            return all_studies
        except Exception as e:
            logger.warning(f"Warning: Exception during fetching studies: {e!s}")
            return []

    def parse_study_data(self, studies: list[dict[str, Any]]) -> list[dict[str, Any]]:
        """Parse the study data to extract relevant information

        Args:
            studies (list[dict[str, Any]]): A list of dictionaries containing study data

        Returns:
            list[dict[str, Any]]: A list of dictionaries with parsed study data
        """
        parsed_data = []
        for study in studies:
            nct_id = (
                study.get("protocolSection", {})
                .get("identificationModule", {})
                .get("nctId", "N/A")
            )
            official_title = (
                study.get("protocolSection", {})
                .get("identificationModule", {})
                .get("officialTitle", "N/A")
            )
            status = (
                study.get("protocolSection", {})
                .get("statusModule", {})
                .get("overallStatus", "N/A")
            )
            conditions = (
                study.get("protocolSection", {})
                .get("conditionsModule", {})
                .get("conditions", [])
            )
            interventions = (
                study.get("protocolSection", {})
                .get("armsInterventionsModule", {})
                .get("interventions", [])
            )
            brief_summary = (
                study.get("protocolSection", {})
                .get("descriptionModule", {})
                .get("briefSummary", "N/A")
            )
            detailed_description = (
                study.get("protocolSection", {})
                .get("descriptionModule", {})
                .get("detailedDescription", "N/A")
            )
            primary_outcomes = (
                study.get("protocolSection", {})
                .get("outcomesModule", {})
                .get("primaryOutcomes", [])
            )
            secondary_outcomes = (
                study.get("protocolSection", {})
                .get("outcomesModule", {})
                .get("secondaryOutcomes", [])
            )

            parsed_data.append(
                {
                    "NCT ID": nct_id,
                    "Official Title": official_title,
                    "Status": status,
                    "Conditions": conditions,
                    "Interventions": interventions,
                    "Brief Summary": brief_summary,
                    "Detailed Description": detailed_description,
                    "Primary Outcomes": primary_outcomes,
                    "Secondary Outcomes": secondary_outcomes,
                }
            )

        return parsed_data

    def _tokenize_and_split_chunks(
        self, chunks: list[str], chunk_size: int, batch_size: int = 5
    ) -> list[str]:
        """
        Tokenize each chunk and split any chunks that exceed the specified token limit.
        Process chunks in batches to reduce memory usage.
        The chunking is done with a 20% chunk size overlap to ensure no information is lost.
        In addition the chunking is done naively, i.e. it does not look for sentence boundaries.

        Args:
            chunks (list[str]): List of text chunks to process
            chunk_size (int): Maximum number of tokens per chunk
            batch_size (int, optional): Number of chunks to process in each batch. Default is 5.

        Returns:
            list[str]: List of processed chunks, each within the specified token limit

        Raises:
            ValueError: If chunk_size is less than 1
        """
        logger.info(
            f"Tokenizing and splitting {len(chunks)} chunks with max size {chunk_size} tokens"
        )
        chunks = [str(chunk) for chunk in chunks]

        target_size = int(0.8 * chunk_size)
        processed_chunks = []

        try:
            for batch_idx in range(0, len(chunks), batch_size):
                batch = chunks[batch_idx : batch_idx + batch_size]
                logger.debug(
                    f"Processing batch {batch_idx//batch_size + 1}/{(len(chunks) + batch_size - 1)//batch_size}"
                )

                batch_results = []
                for chunk_idx, chunk in enumerate(batch):
                    tokens = self.model.tokenizer.encode(chunk, add_special_tokens=True)

                    if len(tokens) <= chunk_size:
                        batch_results.append(chunk)
                        logger.debug(
                            f"Chunk {batch_idx + chunk_idx + 1} is within token limit ({len(tokens)}/{chunk_size})"
                        )
                        del tokens
                        continue

                    logger.debug(
                        f"Chunk {batch_idx + chunk_idx + 1} exceeds token limit ({len(tokens)}/{chunk_size}), splitting..."
                    )

                    start_idx = 0
                    sub_chunk_count = 0

                    while start_idx < len(tokens):
                        end_idx = min(start_idx + target_size, len(tokens))

                        logger.debug(
                            f"Sub-chunk {sub_chunk_count+1}: Processing from token {start_idx} to {end_idx} ({end_idx-start_idx} tokens)"
                        )

                        sub_chunk_tokens = tokens[start_idx:end_idx]
                        sub_chunk = self.model.tokenizer.decode(
                            sub_chunk_tokens, skip_special_tokens=True
                        )
                        batch_results.append(sub_chunk)
                        sub_chunk_count += 1

                        if end_idx >= len(tokens):
                            break

                        old_start_idx = start_idx
                        overlap = min(100, end_idx - start_idx)
                        start_idx = end_idx - overlap

                        logger.debug(
                            f"Sub-chunk {sub_chunk_count}: Added {len(sub_chunk)} chars, moved start_idx from {old_start_idx} to {start_idx} (overlap: {overlap} tokens)"
                        )

                        if start_idx <= old_start_idx:
                            logger.error(
                                f"Loop not progressing! start_idx={start_idx}, old_start_idx={old_start_idx}, end_idx={end_idx}"
                            )
                            logger.error(
                                "Breaking infinite loop, please check the algorithm logic"
                            )
                            break

                    logger.debug(
                        f"Split chunk {batch_idx + chunk_idx + 1} into {sub_chunk_count} smaller chunks"
                    )
                    del tokens
                processed_chunks.extend(batch_results)
                del batch_results
                del batch
                gc.collect()

        except Exception as e:
            logger.error(f"Error during tokenization and splitting: {e}")
            raise

        logger.info(
            f"Completed tokenization and splitting: {len(chunks)} input chunks → {len(processed_chunks)} output chunks"
        )
        return processed_chunks

    async def search_clinical_trials(
        self, query: str, drug_name: str = None, top_k: int = 10
    ) -> list[dict]:
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
            List[Dict]: List of dictionaries with relevant clinical trial data
        """
        search_term = drug_name if drug_name else query

        studies = await self.fetch_all_studies(search_term)
        parsed_studies = self.parse_study_data(studies)

        try:
            chunks = []
            for study in parsed_studies:
                study_text = (
                    f"Title: {study['Official Title']}\n"
                    f"Summary: {study['Brief Summary']}\n"
                    f"Status: {study['Status']}\n"
                    f"Conditions: {', '.join(study['Conditions'])}\n"
                    f"Description: {study['Detailed Description']}"
                )
                chunks.append(study_text)

            if not chunks:
                return []

            # The E5 model has a maximum context length of 512 tokens
            chunks = self._tokenize_and_split_chunks(
                chunks, chunk_size=512, batch_size=2
            )

            original_study_indices = [
                i // (len(chunks) // len(parsed_studies)) for i in range(len(chunks))
            ]
            original_study_indices = [
                min(idx, len(parsed_studies) - 1) for idx in original_study_indices
            ]

            task = "Given a clinical trials database, retrieve relevant trials that answer the query"
            formatted_query = f"Instruct: {task}\nQuery: {query}"

            chunk_embeddings = self.model.encode(chunks, normalize_embeddings=True)
            query_embedding = self.model.encode(
                [formatted_query], normalize_embeddings=True
            )[0]

            similarities = cosine_similarity([query_embedding], chunk_embeddings)[0]

            top_indices = np.argsort(similarities)[-top_k:][::-1]

            seen_study_indices = set()
            results = []

            for idx in top_indices:
                study_idx = (
                    original_study_indices[idx]
                    if idx < len(original_study_indices)
                    else 0
                )

                if study_idx in seen_study_indices:
                    continue

                seen_study_indices.add(study_idx)

                results.append(
                    {
                        "score": float(similarities[idx]),
                        "study": parsed_studies[study_idx],
                    }
                )

                if len(results) >= top_k:
                    break

            return results

        except Exception as e:
            return [{"error": f"Error searching clinical trials: {e!s}"}]
