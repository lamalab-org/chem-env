from modal import Image
import os


_converters_image = (
    Image.debian_slim(python_version="3.12")
    .pip_install(
        [
            "rdkit",
            "selfies",
            "deepsmiles",
            "aiohttp",
            "backoff",
            "loguru",
        ]
    )
    .env({"PRIVATE_API_URL": os.environ.get("PRIVATE_API_URL", "")})
)

with _converters_image.imports():
    import os
    from loguru import logger
    import aiohttp
    import backoff
    import asyncio
    from typing import Optional
    from urllib.parse import quote, unquote
    from rdkit import Chem
    import deepsmiles
    import selfies


_safe_image = Image.debian_slim().pip_install("safe-mol")

with _safe_image.imports():
    pass


class Name2Smiles:
    """Convert chemical names to SMILES notation using multiple chemical APIs.

    This class attempts to convert chemical compound names to SMILES notation
    by querying multiple chemical databases APIs in parallel until a valid
    result is found.

    Args:
        name: The chemical compound name to convert to SMILES notation.

    Raises:
        ValueError: If the name cannot be URL-encoded or contains invalid characters.
    """

    def __init__(self, name: str):
        """Initialize converter with chemical name.

        Args:
            name: Chemical compound name to convert

        Raises:
            ValueError: If name cannot be URL-encoded

        Example:
        # Basic usage with IUPAC name
        >>> converter = Name2Smiles("2-propanone")
        >>> await converter.get_smiles()
        'CC(=O)C'
        """
        try:
            self.name = quote(name)
        except Exception as e:
            logger.error(f"Error encoding name: {e}")
            raise ValueError(f"Invalid chemical name: {name}")

        self.timeout = 10  # seconds

    @backoff.on_exception(
        backoff.expo,
        (aiohttp.ClientError, asyncio.TimeoutError),
        max_time=10,
        logger=logger,
    )
    async def opsin_api(self) -> Optional[str]:
        """
        Queries the OPSIN (Open Parser for Systematic IUPAC Nomenclature) API
        to convert chemical name to SMILES.

        Returns:
            str: SMILES notation if successful, None otherwise

        Raises:
            aiohttp.ClientError: On API connection errors
            asyncio.TimeoutError: If request times out
        """
        url = f"https://opsin.ch.cam.ac.uk/opsin/{self.name}"

        async with aiohttp.ClientSession() as session:
            try:
                async with session.get(url, timeout=self.timeout) as response:
                    if response.status == 200:
                        data = await response.json()
                        return data["smiles"]
                    raise ValueError(f"OPSIN API failed with status {response.status}")
            except Exception as e:
                raise e

    @backoff.on_exception(
        backoff.expo,
        (aiohttp.ClientError, asyncio.TimeoutError),
        max_time=10,
        logger=logger,
    )
    async def cactus(self) -> Optional[str]:
        """
        Queries the NCI CACTUS Chemical Identifier Resolver API
        to convert chemical name to SMILES.

        Returns:
            str: SMILES notation if successful, None otherwise

        Raises:
            aiohttp.ClientError: On API connection errors
            asyncio.TimeoutError: If request times out
        """
        url = f"https://cactus.nci.nih.gov/chemical/structure/{self.name}/smiles"

        async with aiohttp.ClientSession() as session:
            try:
                async with session.get(url, timeout=self.timeout) as response:
                    if response.status == 200:
                        return await response.text()
                    raise ValueError(f"CACTUS API failed with status {response.status}")
            except Exception as e:
                raise e

    @backoff.on_exception(
        backoff.expo,
        (aiohttp.ClientError, asyncio.TimeoutError),
        max_time=10,
        logger=logger,
    )
    async def pubchem(self) -> Optional[str]:
        """
        Queries the PubChem REST API to convert chemical name to
        isomeric SMILES notation.

        Returns:
            str: SMILES notation if successful, None otherwise

        Raises:
            aiohttp.ClientError: On API connection errors
            asyncio.TimeoutError: If request times out
        """
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{self.name}/property/IsomericSMILES/JSON"

        async with aiohttp.ClientSession() as session:
            try:
                async with session.get(url, timeout=self.timeout) as response:
                    if response.status == 200:
                        data = await response.json()
                        return data["PropertyTable"]["Properties"][0]["IsomericSMILES"]
                    raise ValueError(
                        f"PubChem API failed with status {response.status}"
                    )
            except Exception as e:
                raise e

    @backoff.on_exception(
        backoff.expo,
        (aiohttp.ClientError, asyncio.TimeoutError),
        max_time=10,
        logger=logger,
    )
    async def unknown(self) -> Optional[str]:
        """
        Queries the Unknown API to convert chemical name to SMILES.

        Returns:
            str: SMILES notation if successful, None otherwise

        Raises:
            aiohttp.ClientError: On API connection errors
            asyncio.TimeoutError: If request times out
        """
        base_url = os.environ.get("PRIVATE_API_URL")
        if not base_url:
            logger.error("Unknown API URL not set")
            return None

        url = f"{base_url}name=%7B{self.name}%7D&what=smiles"

        async with aiohttp.ClientSession() as session:
            try:
                async with session.get(url, timeout=self.timeout) as response:
                    if response.status == 200:
                        message_content = await response.text()
                        message_text = message_content.split("Message:")[1].strip()
                        return message_text
                    raise ValueError(
                        f"Unknown API failed with status {response.status}"
                    )
            except Exception as e:
                raise e

    async def get_smiles(self) -> Optional[str]:
        """Query all APIs in parallel until a valid SMILES is found.

        Attempts to convert name to SMILES using multiple APIs concurrently,
        returning the first successful result.

        Returns:
            str: First valid SMILES notation found

        Raises:
            ValueError: If no API returns a valid SMILES notation
        """
        tasks = [
            self.opsin_api(),
            self.cactus(),
            self.pubchem(),
            self.unknown(),
        ]

        for result in asyncio.as_completed(tasks):
            try:
                smiles = await result
                if smiles:
                    return smiles.strip()
            except Exception:
                continue

        logger.error(f"Could not find SMILES for {unquote(self.name)}")
        raise ValueError(f"Could not find SMILES for {unquote(self.name)}")


class Smiles2Name:
    """
    Convert SMILES chemical notation to IUPAC chemical names.

    This class provides methods to query different chemical databases (PubChem, CACTUS)
    to obtain IUPAC names for chemical compounds using their SMILES representation.

    Args:
        smiles (str): The SMILES string representing the chemical compound.

    Raises:
        ValueError: If the SMILES string is invalid or cannot be encoded.
    """

    def __init__(self, smiles):
        """Initialize Name2Smiles converter with a chemical compound name.

        Takes a chemical compound name and prepares it for API queries by URL-encoding.
        Sets default timeout for API requests to 10 seconds.

        Args:
            name (str): Chemical compound name to convert to SMILES notation.
                Should be a valid IUPAC or common chemical name.

        Raises:
            ValueError: If the name cannot be URL-encoded or contains invalid characters.

        Example:
            >>> converter = Name2Smiles("ethanol")
            >>> await converter.get_smiles()
            'CCO'
        """
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError(f"Invalid SMILES: {smiles}")

        self.smiles = smiles
        self.timeout = 10  # seconds

    @backoff.on_exception(
        backoff.expo,
        (aiohttp.ClientError, asyncio.TimeoutError),
        max_time=10,
        logger=logger,
    )
    async def pubchem(self) -> Optional[str]:
        """
        Query PubChem API to get IUPAC name from SMILES.

        Returns:
            Optional[str]: IUPAC name if found, None if the query failed.

        Raises:
            aiohttp.ClientError: If the API request fails.
            asyncio.TimeoutError: If the request times out.
        """
        smiles = quote(self.smiles)
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/{smiles}/property/IUPACName/TXT"
        async with aiohttp.ClientSession() as session:
            try:
                async with session.get(url, timeout=self.timeout) as response:
                    if response.status == 200:
                        return await response.text()
                    raise ValueError(
                        f"PubChem API failed with status {response.status}"
                    )
            except Exception as e:
                raise e

    @backoff.on_exception(
        backoff.expo,
        (aiohttp.ClientError, asyncio.TimeoutError),
        max_time=10,
        logger=logger,
    )
    async def cactus(self) -> Optional[str]:
        """
        Query CACTUS API to get IUPAC name from SMILES.

        Returns:
            Optional[str]: IUPAC name if found, None if the query failed.

        Raises:
            aiohttp.ClientError: If the API request fails.
            asyncio.TimeoutError: If the request times out.
        """
        inchi = Chem.MolToInchi(Chem.MolFromSmiles(self.smiles))
        url = f"https://cactus.nci.nih.gov/chemical/structure/{inchi}/iupac_name"

        async with aiohttp.ClientSession() as session:
            try:
                async with session.get(url, timeout=self.timeout) as response:
                    if response.status == 200:
                        return await response.text()
                    raise ValueError(f"CACTUS API failed with status {response.status}")
            except Exception as e:
                raise e

    async def get_name(self) -> Optional[str]:
        """
        Query multiple chemical APIs in parallel to get IUPAC name.

        Attempts to retrieve the IUPAC name by querying multiple chemical databases
        concurrently (CACTUS and PubChem). Returns the first successful result.

        Returns:
            str: The IUPAC name of the chemical compound.

        Raises:
            ValueError: If no name could be found in any of the chemical databases.
        """
        tasks = [
            self.cactus(),
            self.pubchem(),
        ]

        for result in asyncio.as_completed(tasks):
            try:
                name = await result
                if name:
                    return name.strip()
            except Exception:
                continue

        logger.error(f"Could not find name for {self.smiles}")
        raise ValueError(f"Could not find name for {self.smiles}")


def smiles_to_selfies(smiles: str) -> str:
    """
    Takes a SMILES and return the SELFIES encoding.

    Args:
        smiles: SMILES string

    Returns:
        str: SELFIES of the input SMILES
    """

    return selfies.encoder(smiles)


def smiles_to_deepsmiles(smiles: str) -> str:
    """
    Takes a SMILES and return the DeepSMILES encoding.

    Args:
        smiles: SMILES string

    Returns:
        str: DeepSMILES of the input SMILES
    """
    converter = deepsmiles.Converter(rings=True, branches=True)
    return converter.encode(smiles)


def smiles_to_canoncial(smiles: str) -> str:
    """
    Takes a SMILES and return the canonical SMILES.

    Args:
        smiles: SMILES string

    Returns:
        str: Canonical SMILES of the input SMILES
    """
    mol = Chem.MolFromSmiles(smiles)
    return Chem.MolToSmiles(mol)


def smiles_to_inchi(smiles: str) -> str:
    """
    Takes a SMILES and return the InChI.

    Args:
        smiles: SMILES string

    Returns:
        str: InChI of the input SMILES
    """
    mol = Chem.MolFromSmiles(smiles)
    return Chem.MolToInchi(mol)


def smiles_to_safe(smiles: str) -> str:
    """
    Takes a SMILES and return the SAFE (https://github.com/datamol-io/safe).

    Args:
        smiles: SMILES string

    Returns:
        str: SAFE of the input SMILES
    """
    import safe

    return safe.encode(smiles, seed=42, canonical=True, randomize=False)
