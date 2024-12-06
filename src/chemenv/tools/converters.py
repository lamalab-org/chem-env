
import backoff
import safe
import deepsmiles
import pubchempy as pcp
import requests
import selfies
from rdkit import Chem
from STOUT import (
    translate_forward,
    translate_reverse
)
from urllib.parse import (
    quote,
    unquote,
)
from typing import Optional
from loguru import logger
import aiohttp
import asyncio

class Name2Smiles:
    """Convert chemical names to SMILES notation using various APIs."""

    def __init__(self, name: str):
        """Initialize with chemical name.
        
        Args:
            name: Chemical compound name
        """
        try:
            self.name = quote(name)
        except Exception as e:
            logger.error(f"Error encoding name: {e}")
            raise ValueError(f"Invalid chemical name: {name}")
        
        self.timeout = 10  # seconds

    async def open_molecules(self) -> Optional[str]:
        """Query Open Molecules API for SMILES."""
        # Implementation pending API details
        return None
    
    @backoff.on_exception(backoff.expo, (aiohttp.ClientError, asyncio.TimeoutError), max_time=10, logger=logger)
    async def opsin(self) -> Optional[str]:
        """Query OPSIN API for SMILES."""
        url = f"https://opsin.ch.cam.ac.uk/opsin/{self.name}"
        
        async with aiohttp.ClientSession() as session:
            try:
                async with session.get(url, timeout=self.timeout) as response:
                    if response.status == 200:
                        data = await response.json()
                        return data["smiles"]
                    logger.warning(f"OPSIN API failed with status {response.status}")
                    return None
            except Exception as e:
                logger.error(f"OPSIN API error: {e}")
                return None

    @backoff.on_exception(backoff.expo, (aiohttp.ClientError, asyncio.TimeoutError), max_time=10, logger=logger)
    async def cactus(self) -> Optional[str]:
        """Query CACTUS API for SMILES."""
        url = f"https://cactus.nci.nih.gov/chemical/structure/{self.name}/smiles"
        
        async with aiohttp.ClientSession() as session:
            try:
                async with session.get(url, timeout=self.timeout) as response:
                    if response.status == 200:
                        return await response.text()
                    logger.warning(f"CACTUS API failed with status {response.status}")
                    return None
            except Exception as e:
                logger.error(f"CACTUS API error: {e}")
                return None

    @backoff.on_exception(backoff.expo, (aiohttp.ClientError, asyncio.TimeoutError), max_time=10, logger=logger)
    async def pubchem(self) -> Optional[str]:
        """Query PubChem API for SMILES."""
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{self.name}/property/IsomericSMILES/JSON"
        
        async with aiohttp.ClientSession() as session:
            try:
                async with session.get(url, timeout=self.timeout) as response:
                    if response.status == 200:
                        data = await response.json()
                        return data["PropertyTable"]["Properties"][0]["IsomericSMILES"]
                    logger.warning(f"PubChem API failed with status {response.status}")
                    return None
            except Exception as e:
                logger.error(f"PubChem API error: {e}")
                return None

    @backoff.on_exception(backoff.expo, (aiohttp.ClientError, asyncio.TimeoutError), max_time=10, logger=logger)
    async def unknown(self) -> Optional[str]:
        """Query unknown API for SMILES."""
        url = f"http://46.4.119.202:8082/?name=%7B{self.name}%7D&what=smiles"
        
        async with aiohttp.ClientSession() as session:
            try:
                async with session.get(url, timeout=self.timeout) as response:
                    if response.status == 200:
                        message_content = await response.text()
                        message_text = message_content.split('Message:')[1].strip()
                        return message_text
                    logger.warning(f"Unknown API failed with status {response.status}")
                    return None
            except Exception as e:
                logger.error(f"Unknown API error: {e}")
                return None
    

    async def get_smiles(self) -> Optional[str]:
        """Try all APIs in parallel until we get a result."""
        tasks = [
            self.opsin(),
            self.cactus(),
            self.pubchem(),
            self.unknown()
        ]
        
        for result in asyncio.as_completed(tasks):
            try:
                smiles = await result
                if smiles:
                    return smiles
            except Exception as e:
                self.logger.error(f"Error in get_smiles: {e}")
                continue
        
        return translate_reverse(unquote(self.name))
    
CACTUS = "https://cactus.nci.nih.gov/chemical/structure/{0}/{1}"


@backoff.on_exception(backoff.expo, requests.exceptions.RequestException, max_time=10)
def cactus_request_w_backoff(smiles, rep="iupac_name"):
    url = CACTUS.format(smiles, rep)
    response = requests.get(url, allow_redirects=True, timeout=10)
    response.raise_for_status()
    name = response.text
    if "html" in name:
        return None
    return name


def smiles_to_iupac_name(smiles: str) -> str:
    """Use the chemical name resolver https://cactus.nci.nih.gov/chemical/structure.
    If this does not work, use pubchem.
    """
    try:
        name = cactus_request_w_backoff(smiles, rep="iupac_name")
        if name is None:
            raise Exception
        return name
    except Exception:
        try:
            compound = pcp.get_compounds(smiles, "smiles")
            return compound[0].iupac_name
        except Exception:
            return None

def smiles2iupac(smiles: str) -> str:
    """
    Convert a SMILES string to an IUPAC name.

    Args:
        smiles: The SMILES string to convert.

    Returns:
        The IUPAC name of the molecule.
    """
    iupac_name = smiles_to_iupac_name(smiles)
    if iupac_name is None:
        return translate_forward(smiles)
    else:
        return iupac_name


def smiles_to_selfies(smiles: str) -> str:
    """
    Takes a SMILES and return the selfies encoding.
    """

    return selfies.encoder(smiles)


def smiles_to_deepsmiles(smiles: str) -> str:
    """
    Takes a SMILES and return the DeepSMILES encoding.
    """
    converter = deepsmiles.Converter(rings=True, branches=True)
    return converter.encode(smiles)


def smiles_to_canoncial(smiles: str) -> str:
    """
    Takes a SMILES and return the canoncial SMILES.
    """
    mol = Chem.MolFromSmiles(smiles)
    return Chem.MolToSmiles(mol)


def smiles_to_inchi(smiles: str) -> str:
    """
    Takes a SMILES and return the InChI.
    """
    mol = Chem.MolFromSmiles(smiles)
    return Chem.MolToInchi(mol)


def smiles_to_safe(smiles: str) -> str:
    """
    Takes a SMILES and return the SAFE.
    """
    return safe.encode(smiles, seed=42, canonical=True, randomize=False)
