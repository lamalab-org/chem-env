from modal import Image
from typing import Optional, List, Dict, Any
import backoff
import aiohttp
import asyncio
import pubchempy as pcp
from rdkit import Chem
from urllib.parse import quote
from time import sleep
from loguru import logger


pubchem_image = Image.debian_slim(python_version="3.12").pip_install(
    "backoff",
    "asyncio",
    "aiohttp",
    "pubchempy",
    "loguru",
    "rdkit",
)
with pubchem_image.imports():
    import backoff
    import aiohttp
    import asyncio
    from typing import Dict, Any, Optional
    import pubchempy as pcp
    from rdkit import Chem
    from urllib.parse import quote
    from time import sleep
    from loguru import logger


class PubChem:
    """
    PubChem handler to retrieve compound data from PubChem using PubChemPy and RDKit.
    Automatically converts the input to a CID and retrieves the data for that compound from PubChem.

    Example:
        >>> pubchem = await PubChem.create("2244")
        >>> isomers = await pubchem._get_number_atoms()
            21
    """

    def __init__(self):
        self.cid = None
        self.base_url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
        self.long_url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/{cid}/JSON/?response_type=display&heading="
        self.complete_url = (
            "https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/{cid}/JSON/"
        )
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

        Args:
            compound (str): Any type of compound identifier (CID, SMILES, InChI, etc.)

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
        """
        Fetch compound data from PubChem REST API using a specific URL

        Args:
            record_url (str): URL to fetch compound data from

        Returns:
            dict: Compound data in JSON format
        """
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

        Args:
            record_url (str): URL to fetch compound data from

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

    async def _get_compound_cid(self) -> Optional[str]:
        """
        Get the PubChem CID for a compound.

        Returns:
            int: PubChem CID of the compound.
        """
        return self.cid

    async def _get_pubchem_full_record(self) -> dict:
        """
        Get the full PubChem record for a compound.

        Returns:
            dict: Full PubChem record for the compound.
        """
        return await self.get_data_from_url(self.complete_url.format(cid=self.cid))

    async def _get_number_atoms(self) -> Optional[int]:
        """
        Get the number of atoms in a compound using RDKit.

        Returns:
            int: Number of atoms in the compound.

        Raises:
            ValueError: If no atoms are found

        Example:
            >>> await self._get_number_atoms()
                21
        """
        try:
            mol = Chem.MolFromSmiles(await self._get_canonical_smiles())
            return mol.GetNumAtoms()
        except Exception as e:
            logger.error(f"No atoms found. {e}")
            raise ValueError(f"No atoms found. {e}")

    async def _get_isomeric_smiles(self) -> Optional[str]:
        """
        Get the isomeric SMILES for a compound from PubChem.

        Returns:
            str: Isomeric SMILES of the compound.

        Raises:
            ValueError: If the isomeric SMILES could not be retrieved

        Example:
            >>> await self._get_isomeric_smiles()
                'CCO'
        """
        url = [
            "/compound/cid/",
            "/property/IsomericSMILES/JSON",
        ]
        logger.info(f"Getting isomeric SMILES for CID {self.cid}")
        try:
            data = (await self.get_compound_data(url))["PropertyTable"]["Properties"][
                0
            ]["SMILES"]
        except Exception as e:
            logger.error(f"Failed to extract Isomeric SMILES: {str(e)}")
            raise ValueError(f"Failed to extract Isomeric SMILES: {str(e)}")
        logger.info(f"Isomeric SMILES for CID {self.cid}: {data}")
        return data

    async def _get_canonical_smiles(self) -> Optional[str]:
        """
        Get the canonical SMILES for a compound from PubChem.

        Returns:
            str: Canonical SMILES of the compound.

        Raises:
            ValueError: If the canonical SMILES could not be retrieved

        Example:
            >>> await self._get_canonical_smiles()
                'CCO'
        """
        url = [
            "/compound/cid/",
            "/property/CanonicalSMILES/JSON",
        ]
        logger.info(f"Getting canonical SMILES for CID {self.cid}")
        try:
            data = (await self.get_compound_data(url))["PropertyTable"]["Properties"][
                0
            ]["ConnectivitySMILES"]
        except Exception as e:
            logger.error(f"Failed to extract Canonical SMILES: {str(e)}")
            raise ValueError(f"Failed to extract Canonical SMILES: {str(e)}")
        logger.info(f"Canonical SMILES for CID {self.cid}: {data}")
        return data

    async def _get_compound_mass(self) -> Optional[float]:
        """
        Get the molecular weight of a compound from PubChem.

        Returns:
            float: Molecular weight of the compound.

        Raises:
            ValueError: If the molecular weight could not be retrieved

        Example:
            >>> await self._get_compound_mass()
                46.07
        """
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
        """
        Get the charge of a compound from PubChem.

        Returns:
            int: Charge of the compound.

        Raises:
            ValueError: If the charge could not be retrieved

        Example:
            >>> await self._get_compound_charge()
                0
        """
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
        """
        Get the formula of a compound from PubChem.

        Returns:
            str: Formula of the compound.

        Raises:
            ValueError: If the formula could not be retrieved

        Example:
            >>> await self._get_compound_formula()
                'C2H6O'
        """
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
        """
        Get the number of compound isomers for a compound from PubChem.

        Returns:
            int: Number of compound isomers.

        Raises:
            ValueError: If the number of isomers could not be retrieved

        Example:
            >>> await self._get_number_isomers()
                2
        """
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

    async def _get_compound_isomers(self, limit: int = 10) -> List[str]:
        """
        Get the compound isomers for a compound from PubChem.
        This function can take some time depending on the number of isomers.
        Returns a maximum of 10 isomers by default.
        Shuffles the results each time to provide different isomers.

        Args:
            limit (int, optional): Maximum number of isomers to return. Defaults to 10.

        Returns:
            list: List of compound isomers (limited to the specified number).

        Raises:
            ValueError: If the isomers could not be retrieved

        Example:
            >>> await self._get_compound_isomers()
                ['CCO', 'COC']
            >>> await self._get_compound_isomers(limit=5)
                ['COC', 'CCO', 'CC[O-]', 'C[CH-]O', 'C[O-]C']
        """
        import random

        url = [
            "/compound/fastformula/",
            "/cids/JSON",
        ]
        logger.info(f"Getting compound isomers for CID {self.cid} (limit: {limit})")
        formula = await self._get_compound_formula()
        url = self.base_url + url[0] + quote(str(formula)) + url[1]
        try:
            isomers_cids = (await self.get_data_from_url(url))["IdentifierList"]["CID"]
            # Shuffle the isomers to get different ones each time
            random.shuffle(isomers_cids)
            # Limit the number of isomers to process
            isomers_cids = isomers_cids[:limit]
            data = []
            for i in isomers_cids:
                sleep(0.5)
                self.cid = i
                try:
                    # Isomeric SMILES to capture enantiomers
                    smiles = await self._get_isomeric_smiles()
                    if smiles:
                        data.append(smiles)
                except ValueError as ve:
                    logger.warning(f"Could not get SMILES for CID {i}: {ve}")
                    continue
                except Exception as e:
                    logger.error(f"Unexpected error getting SMILES for CID {i}: {e}")
                    continue
        except Exception as e:
            logger.error(f"No compound isomers found. {e}")
            raise ValueError(f"No compound isomers found. {e}")
        logger.info(f"Found {len(data)} compound isomers for the given formula")
        return data

    @classmethod
    async def _get_compound_isomers_from_formula(
        cls, formula: str, limit: int = 10
    ) -> List[str]:
        """
        Get the compound isomers for a compound from PubChem from a given formula.
        This function can take some time depending on the number of isomers.
        Returns a maximum of 10 isomers by default.
        Shuffles the results each time to provide different isomers.

        Args:
            formula (str): The empirical formula of the compound.
            limit (int, optional): Maximum number of isomers to return. Defaults to 10.

        Returns:
            list: List of compound isomers (limited to the specified number).

        Raises:
            ValueError: If the isomers could not be retrieved
        """
        instance = cls()
        url = [
            "/compound/fastformula/",
            "/cids/JSON",
        ]
        import random

        url = instance.base_url + url[0] + quote(str(formula)) + url[1]
        try:
            isomers_cids = (await instance.get_data_from_url(url))["IdentifierList"][
                "CID"
            ]
            # Shuffle the isomers to get different ones each time
            random.shuffle(isomers_cids)
            # Limit the number of isomers to process
            isomers_cids = isomers_cids[:limit]
            data = []
            for i in isomers_cids:
                sleep(0.5)
                instance.cid = i
                try:
                    # Isomeric SMILES to capture enantiomers
                    smiles = await instance._get_isomeric_smiles()
                    if smiles:
                        data.append(smiles)
                except ValueError as ve:
                    logger.warning(f"Could not get SMILES for CID {i}: {ve}")
                    continue
                except Exception as e:
                    logger.error(f"Unexpected error getting SMILES for CID {i}: {e}")
                    continue
        except Exception as e:
            logger.error(f"No compound isomers found. {e}")
            raise ValueError(f"No compound isomers found. {e}")
        logger.info(f"Found {len(data)} compound isomers for the given formula")
        return data

    async def _get_number_heavy_atoms(self) -> Optional[int]:
        """
        Get the number of heavy atoms in a compound from PubChem.

        Returns:
            int: Number of heavy atoms in the compound.

        Raises:
            ValueError: If the number of heavy atoms could not be retrieved

        Example:
            >>> await self._get_number_heavy_atoms()
                3
        """
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
        """
        Get the number of chiral atoms in a compound from PubChem.

        Returns:
            int: Number of chiral atoms in the compound.

        Raises:
            ValueError: If the number of chiral atoms could not be retrieved

        Example:
            >>> await self._get_number_chiral_atoms()
                1
        """
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

    async def _format_long_url(self, heading: str) -> Dict[str, Any]:
        """
        Format the long URL to get specific information from PubChem and return the data.

        Args:
            heading (str): Heading of the information to retrieve from PubChem

        Returns:
            dict: Data retrieved from PubChem

        Raises:
            ValueError: If the data could not be retrieved

        Example:
            >>> await self._format_long_url("Mass Spectrometry")
                {'Information': [{'Name': 'Mass bank ID', 'ReferenceNumber': 1, 'Value': {'StringWithMarkup': [{'String...
        """
        url = self.long_url.format(cid=self.cid) + quote(heading)
        logger.info(f"Getting spectral information for CID {self.cid}")
        try:
            return await self.get_data_from_url(url)
        except Exception as e:
            logger.error(f"Failed to get spectral information: {str(e)}")
            raise ValueError(f"Failed to get spectral information: {str(e)}")

    async def _format_ms_spectra(self, data):
        """
        Format the MS spectra data retrieved from PubChem.

        Args:
            data (dict): Data retrieved from PubChem

        Returns:
            dict: Formatted MS spectra data

        Raises:
            ValueError: If the data could not be retrieved

        Example:
            >>> await self._format_ms_spectra(data)
                {'Mass bank ID': 'MoNA ID', ...}
        """
        try:
            information = data["Information"]
        except KeyError:
            raise ValueError("No MS spectra found")

        field_mapping = {
            "Mass bank ID": "MoNA ID",
            "Spectra type": "MS Category",
            "MS Type": "MS Type",
            "MS Level": "MS Level",
            "Instrument": "Instrument",
            "Instrument Type": "Instrument Type",
            "Ionization Mode": "Ionization Mode",
            "Top Peaks": "StringWithMarkup",
        }

        results = {}

        counter = 0
        for info in information:
            if counter > 5:
                return results
            ref_num = info.get("ReferenceNumber")
            if not ref_num:
                continue

            value = info.get("Value", {})

            for orig_key, mapped_key in field_mapping.items():
                if mapped_key == "StringWithMarkup":
                    string_with_markup = value.get("StringWithMarkup", [])
                    peaks = [item["String"] for item in string_with_markup]
                    results[orig_key] = peaks
                else:
                    string_with_markup = value.get("StringWithMarkup", [])
                    if string_with_markup:
                        results[orig_key] = string_with_markup[0]["String"]
                    else:
                        results[orig_key] = None

            counter += 1

        return results

    async def _format_nmr_spectra(self, data):
        """
        Format the NMR spectra data retrieved from PubChem.

        Args:
            data (dict): Data retrieved from PubChem

        Returns:
            dict: Formatted NMR spectra data

        Raises:
            ValueError: If the data could nt be retrieved

        Example:
            >>> await self._format_nmr_spectra(data)
                {'1': {'instrument': 'Bruker', 'frequency': '400 MHz',...}}
        """
        try:
            information = data["Information"]
        except KeyError:
            raise ValueError("No 1H NMR spectra found")

        # Define field mapping
        field_mapping = {
            "Instrument Type": "instrument",
            "Frequency": "frequency",
            "Solvent": "solvent",
            "pH": "ph",
            "Shifts [ppm]:Intensity": "shifts",
        }

        results = {}

        for info in information:
            ref_num = info.get("ReferenceNumber")
            if not ref_num:
                continue

            name = info.get("Name")
            string_value = (
                info.get("Value", {}).get("StringWithMarkup", [{}])[0].get("String")
            )

            if not (name and string_value):
                continue

            if ref_num not in results:
                results[ref_num] = {}

            if name in field_mapping:
                results[ref_num][field_mapping[name]] = string_value

        return results

    async def _get_c_nmr_spectra(self) -> Dict[str, Any]:
        """
        Get the C-NMR spectra for a compound from PubChem.

        Returns:
            dict: C-NMR spectra data

        Raises:
            ValueError: If the data could not be retrieved

        Example:
            >>> await self._get_c_nmr_spectra()
                {'1': {'instrument': 'Bruker', 'frequency': '400 MHz',...}}
        """
        try:
            data = await self._format_long_url("13C NMR Spectra")
            return self._format_nmr_spectra(data)
        except Exception:
            raise ValueError("No C-NMR spectra found. {e}")

    async def _get_h_nmr_spectra(self) -> Dict[str, Any]:
        """
        Get the 1H NMR spectra for a compound from PubChem.

        Returns:
            dict: 1H NMR spectra data

        Raises:
            ValueError: If the data could not be retrieved

        Example:
            >>> await self._get_h_nmr_spectra()
                {'1': {'instrument': 'Bruker', 'frequency': '400 MHz',...}}
        """
        try:
            data = await self._format_long_url("1H NMR Spectra")
            return self._format_nmr_spectra(data)
        except Exception:
            raise ValueError("No 1H NMR spectra found. {e}")

    async def _get_uv_spectra(self) -> str:
        """
        Get the UV spectra for a compound from PubChem.

        Returns:
            str: UV spectra data

        Raises:
            ValueError: If the data could not be retrieved

        Example:
            >>> await self._get_uv_spectra()
                'Reference 1:\nMAX ABSORPTION: 210 nm\nReference 2:\nUV MAX: 210 nm'
        """
        data = await self._format_long_url("UV Spectra")
        results = {}
        try:
            for info in data["Information"]:
                ref_num = info["ReferenceNumber"]
                string_value = info["Value"]["StringWithMarkup"][0]["String"]

                if ref_num not in results:
                    results[ref_num] = ""

                if (
                    "MAX ABSORPTION" in string_value.upper()
                    or "UV MAX" in string_value.upper()
                ):
                    if results[ref_num]:
                        results[ref_num] += "\n"
                    results[ref_num] += string_value

            if not results:
                raise ValueError("No UV spectra found")

            output = []
            for ref_num, value in sorted(results.items()):
                output.append(f"Reference {ref_num}:\n{value}")

            return "\n\n".join(output)

        except Exception:
            raise ValueError("No UV spectra found")

    async def _get_ms_spectra(self) -> Dict[str, Any]:
        """
        Get the MS spectra for a compound from PubChem.

        Returns:
            dict: MS spectra data

        Raises:
            ValueError: If the data could not be retrieved

        Example:
            >>> await self._get_ms_spectra()
                {'Mass bank ID': 'MoNA ID', ...}
        """
        try:
            data = await self._format_long_url("Mass Spectrometry")
            results = {}
            for section in data["Record"]["Section"]:
                title = section["TOCHeading"]
                spectra = self._format_ms_spectra(section)
                if spectra:
                    results[title] = spectra

            if not results:
                raise ValueError("No MS spectra found")

            return results

        except Exception:
            raise ValueError("No MS spectra found")

    async def _get_ghs_classification(self) -> Dict[str, List[str]]:
        """
        Get the GHS classification for a compound from PubChem.

        Returns:
            dict: GHS classification data

        Raises:
            ValueError: If the data could not be retrieved

        Example:
            >>> await self._get_ghs_classification()
                {'H225': ['Highly flammable liquid and vapour'], ...}
        """
        data = await self._format_long_url("GHS%20Classification")
        logger.info(f"Getting GHS classification for CID {self.cid}")
        try:
            information_list = data["Record"]["Section"][0]["Section"][0]["Section"][0][
                "Information"
            ]

            hazard_statements = {}

            for info in information_list:
                if info.get("Name") == "GHS Hazard Statements":
                    ref_number = info.get("ReferenceNumber")
                    string_values = [
                        markup["String"] for markup in info["Value"]["StringWithMarkup"]
                    ]
                    hazard_statements[ref_number] = string_values

            return hazard_statements
        except Exception:
            logger.error("Failed to get GHS classification")
            raise ValueError("Failed to get GHS classification")

    async def _get_patent_count(self) -> int:
        """
        Get the number of patents for a compound from PubChem.

        Returns:
            int: Number of patents for the compound.

        Raises:
            ValueError: If the data could not be retrieved

        Example:
            >>> await self._get_patent_count()
                10
        """
        url = [
            "/compound/cid/",
            "/property/PatentCount/JSON",
        ]
        logger.info(f"Getting number of patents for CID {self.cid}")
        try:
            data = (await self.get_compound_data(url))["PropertyTable"]["Properties"][
                0
            ]["PatentCount"]
        except Exception as e:
            logger.error(f"No patents found. {e}")
            raise ValueError(f"No patents found. {e}")
        logger.info(f"Number of patents for CID {self.cid}: {data}")
        return data

    async def _return_physical_property(self):
        """
        Get the physical properties for a compound from PubChem.

        Returns:
            dict: Physical properties data

        Raises:
            ValueError: If the data could not be retrieved

        Example:
            >>> await self.return_physical_property()
                {'Experimental Properties': ['Appearance: clear colorless liquid', 'Boiling Point: 78.37...']}
        """
        data = await self._format_long_url("Experimental%20Properties")
        results = {}
        for section in data["Record"]["Section"]:
            heading = section["TOCHeading"]

            results[heading] = []
            for info in data["Information"]:
                for string_markup in info["Value"]["StringWithMarkup"]:
                    results[heading].append(string_markup["String"])

        return results
