import os
from modal import Image

converters_image = (
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

with converters_image.imports():
    import backoff
    from rdkit import Chem
    from urllib.parse import quote
    from typing import Optional, List
    from loguru import logger
    import aiohttp
    import asyncio


pubchem_image = Image.debian_slim().pip_install(
    "backoff",
    "asyncio",
    "aiohttp",
    "pubchempy",
    "loguru",
)
with pubchem_image.imports():
    import backoff
    import aiohttp
    import asyncio
    from typing import Dict, Any, Optional
    import pubchempy as pcp
    from urllib.parse import quote
    from loguru import logger
"""
Full
{base_url}/compound/cid/{self.cid}/record/JSON

Precise
{base_url}/compound/cid/{self.cid}/record/JSON
"""


class PubChem:
    def __init__(self):
        """
        Initialize compound handler with any type of compound identifier.
        Automatically converts the input to a CID and retrieves the data for that compound from PubChem.

        Example:
            >>> pubchem = await PubChem.create("2244")
            >>> isomers = await pubchem._get_number_atoms()
            21
        """
        self.cid = None
        self.base_url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
        logger.info(f"Initialized PubChem handler for CID {self.cid}")

    @classmethod
    async def create(cls, compound: str):
        """
        Factory method to create a PubChem handler with a compound identifier.
        Automatically converts the input to a CID and retrieves the data for that compound from PubChem.

        Args:
            compound (str): Any type of compound identifier (CID, SMILES, InChI, etc.)

        Returns:
            PubChem: Instance of PubChem handler with the compound data
        """
        self = cls()
        self.cid = await self._normalize_to_cid(compound)
        logger.info(f"Initialized PubChem handler for CID {self.cid}")
        return self

    @backoff.on_exception(
        backoff.expo, (aiohttp.ClientError, ValueError), max_tries=10, max_time=30
    )
    async def _normalize_to_cid(self, compound: str) -> Optional[str]:
        """
        Convert any compound identifier to PubChem CID asynchronously.
        Tries different approaches to identify the compound format and get its CID.

        Returns:
            int: PubChem CID if found, None if not found
        """
        try:
            if compound.isdigit():
                results = pcp.get_compounds(compound, "cid")
                if results:
                    return compound
        except Exception as e:
            logger.error(f"Invalid compound CID: {compound}, {e}")
            raise ValueError("Invalid compound CID")

        for namespace in ["name", "smiles", "inchi"]:
            try:
                results = pcp.get_compounds(compound, namespace)
                if results:
                    return results[0].cid
            except Exception:
                continue

        logger.error(f"Invalid compound identifier: {compound}")
        raise ValueError(
            "Invalid compound identifier. Only name, smiles or InChI are supported."
        )

    @backoff.on_exception(
        backoff.expo,
        (aiohttp.ClientError, asyncio.TimeoutError),
        max_tries=10,
        max_time=30,
    )
    async def get_data_from_url(
        self, record_url: Optional[str] = None
    ) -> Dict[Any, Any]:
        if not record_url:
            raise ValueError("No record URL provided.")
        timeout = aiohttp.ClientTimeout(total=10)
        async with aiohttp.ClientSession(timeout=timeout) as session:
            try:
                async with session.get(record_url) as record_response:
                    record_response.raise_for_status()
                    record_data = await record_response.json()
                return record_data
            except Exception as e:
                raise e

    @backoff.on_exception(
        backoff.expo,
        (aiohttp.ClientError, asyncio.TimeoutError),
        max_tries=10,
        max_time=30,
    )
    async def get_compound_data(
        self, record_url: Optional[str] = None
    ) -> Dict[Any, Any]:
        """
        Fetch compound data from PubChem REST API combining basic record and detailed data

        Returns:
            dict: Combined compound data in JSON format containing both record and detailed information
        """
        if not record_url:
            raise ValueError("No record URL provided.")
        record_url = (
            self.base_url + record_url[0] + str(self.cid) + record_url[1]
            if record_url
            else None
        )
        timeout = aiohttp.ClientTimeout(total=10)

        async with aiohttp.ClientSession(timeout=timeout) as session:
            try:
                async with session.get(record_url) as record_response:
                    record_response.raise_for_status()
                    record_data = await record_response.json()
                return record_data

            except Exception as e:
                raise e

    async def _get_number_atoms(self) -> Optional[int]:
        url = [
            "/compound/cid/",
            "/record/JSON",
        ]
        logger.info(f"Getting number of atoms for CID {self.cid}")
        try:
            data = (await self.get_compound_data(url))["PC_Compounds"][0]["atoms"][
                "aid"
            ]
        except Exception as e:
            logger.error(f"No atoms found. {e}")
            raise ValueError(f"No atoms found. {e}")
        logger.info(f"Number of atoms for CID {self.cid}: {len(data)}")
        return len(data)

    async def _get_isomeric_smiles(self) -> Optional[str]:
        url = [
            "/compound/cid/",
            "/property/IsomericSMILES/JSON",
        ]
        logger.info(f"Getting isomeric SMILES for CID {self.cid}")
        try:
            data = (await self.get_compound_data(url))["PropertyTable"]["Properties"][
                0
            ]["IsomericSMILES"]
        except Exception as e:
            logger.error(f"Failed to extract Isomeric SMILES: {str(e)}")
            raise ValueError(f"Failed to extract Isomeric SMILES: {str(e)}")
        logger.info(f"Isomeric SMILES for CID {self.cid}: {data}")
        return data

    async def _get_canonical_smiles(self) -> Optional[str]:
        url = [
            "/compound/cid/",
            "/property/CanonicalSMILES/JSON",
        ]
        logger.info(f"Getting canonical SMILES for CID {self.cid}")
        try:
            data = (await self.get_compound_data(url))["PropertyTable"]["Properties"][
                0
            ]["CanonicalSMILES"]
        except Exception as e:
            logger.error(f"Failed to extract Canonical SMILES: {str(e)}")
            raise ValueError(f"Failed to extract Canonical SMILES: {str(e)}")
        logger.info(f"Canonical SMILES for CID {self.cid}: {data}")
        return data

    async def _get_compound_mass(self) -> Optional[float]:
        url = [
            "/compound/cid/",
            "/property/MolecularWeight/JSON",
        ]
        logger.info(f"Getting compound mass for CID {self.cid}")
        try:
            data = (await self.get_compound_data(url))["PropertyTable"]["Properties"][
                0
            ]["MolecularWeight"]
        except Exception as e:
            logger.error(f"No compound mass found. {e}")
            raise ValueError(f"No compound mass found. {e}")
        logger.info(f"Compound mass for CID {self.cid}: {data}")
        return data

    async def _get_compound_charge(self) -> Optional[int]:
        url = [
            "/compound/cid/",
            "/property/Charge/JSON",
        ]
        logger.info(f"Getting compound charge for CID {self.cid}")
        try:
            data = (await self.get_compound_data(url))["PropertyTable"]["Properties"][
                0
            ]["Charge"]
        except Exception as e:
            logger.error(f"No compound charge found. {e}")
            raise ValueError(f"No compound charge found. {e}")
        logger.info(f"Compound charge for CID {self.cid}: {data}")
        return data

    async def _get_compound_formula(self) -> Optional[str]:
        url = [
            "/compound/cid/",
            "/property/MolecularFormula/JSON",
        ]
        logger.info(f"Getting compound formula for CID {self.cid}")
        try:
            data = (await self.get_compound_data(url))["PropertyTable"]["Properties"][
                0
            ]["MolecularFormula"]
        except Exception as e:
            logger.error(f"No compound formula found. {e}")
            raise ValueError(f"No compound formula found. {e}")
        logger.info(f"Compound formula for CID {self.cid}: {data}")
        return data

    async def _get_number_isomers(self) -> Optional[int]:
        url = [
            "/compound/fastformula/",
            "/cids/JSON",
        ]
        logger.info(f"Getting number of compound isomers for CID {self.cid}")
        formula = await self._get_compound_formula()
        url = self.base_url + url[0] + quote(str(formula)) + url[1]
        try:
            data = len((await self.get_data_from_url(url))["IdentifierList"]["CID"])
        except Exception as e:
            logger.error(f"No compound isomers found. {e}")
            raise ValueError(f"No compound isomers found. {e}")
        logger.info(f"Number of compound isomers for CID {self.cid}: {data}")
        return data

    async def _get_compound_isomers(self) -> List[Optional[str]]:
        url = [
            "/compound/fastformula/",
            "/cids/JSON",
        ]
        logger.info(f"Getting compound isomers for CID {self.cid}")
        formula = await self._get_compound_formula()
        url = self.base_url + url[0] + quote(str(formula)) + url[1]
        try:
            isomers_cids = (await self.get_data_from_url(url))["IdentifierList"]["CID"][
                :4
            ]
            data = []
            for i in isomers_cids:
                self.cid = i
                data.append(await self._get_isomeric_smiles())
        except Exception as e:
            logger.error(f"No compound isomers found. {e}")
            raise ValueError(f"No compound isomers found. {e}")
        logger.info(f"Compound isomers for CID {self.cid}: {data[5:]}")
        return data

    async def _get_number_heavy_atoms(self) -> Optional[int]:
        url = [
            "/compound/cid/",
            "/property/HeavyAtomCount/JSON",
        ]
        logger.info(f"Getting number of heavy atoms for CID {self.cid}")
        try:
            data = (await self.get_compound_data(url))["PropertyTable"]["Properties"][
                0
            ]["HeavyAtomCount"]
        except Exception as e:
            logger.error(f"No heavy atoms found. {e}")
            raise ValueError(f"No heavy atoms found. {e}")
        logger.info(f"Number of heavy atoms for CID {self.cid}: {data}")
        return data

    async def _get_number_chiral_atoms(self) -> Optional[int]:
        url = [
            "/compound/cid/",
            "/record/JSON",
        ]
        logger.info(f"Getting number of chiral atoms for CID {self.cid}")
        try:
            data = (await self.get_compound_data(url))["PC_Compounds"][0]["count"][
                "atom_chiral"
            ]
        except Exception as e:
            logger.error(f"No chiral atoms found. {e}")
            raise ValueError(f"No chiral atoms found. {e}")
        logger.info(f"Number of chiral atoms for CID {self.cid}: {data}")
        return data


class Smiles2Name:
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
                    return name
            except Exception:
                continue

        logger.error(f"Could not find name for {self.smiles}")
        raise ValueError(f"Could not find name for {self.smiles}")


# if __name__ == "__main__":
# async def main():
# try:
# pubchem = await PubChem.create("2244")
# isomers = await pubchem._get_number_atoms()
# print(isomers)
# if isomers is None:
# return ""
# return isomers
# except Exception as e:
# logger.error(f"Error processing compound: {e}")
# return None

# asyncio.run(main())
