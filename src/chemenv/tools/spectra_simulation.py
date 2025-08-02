from modal import Image
from typing import Dict

spectra_simulation_image = Image.debian_slim(python_version="3.12").pip_install(
    "asyncio",
    "aiohttp",
    "backoff",
    "loguru",
)

with spectra_simulation_image.imports():
    import asyncio
    import aiohttp
    import backoff
    from loguru import logger
    from typing import Dict


class SpectraAPI:
    """Class for handling spectral predictions from the NMR and IR APIs

    Attributes:
        IR_BASE_URL: str - Base URL for the IR prediction API
        NMR_BASE_URL: str - Base URL for the NMR prediction API
    """

    IR_BASE_URL: str = "https://ir.cheminfo.org/v1/ir"
    NMR_BASE_URL: str = "https://nmr-prediction.service.zakodium.com/v1/predict"
    VALID_SPECTRUM_TYPES = {"carbon", "proton"}

    _request_lock = asyncio.Lock()
    _last_request_time = 0
    _min_request_interval = 1
    _timeout = 10

    @staticmethod
    def format_c13_nmr(json_response: Dict) -> str:
        """Format C13 NMR response to literature format

        Args:
            json_response: C13 NMR prediction response

        Returns:
            str: Formatted C13 NMR prediction
        """
        shifts = sorted(
            [signal["delta"] for signal in json_response["data"]["signals"]],
            reverse=True,
        )
        return f"Deltas: {', '.join(f'{shift:.2f}' for shift in shifts)}"

    @staticmethod
    def format_h_nmr(json_response: Dict) -> str:
        """Format 1H NMR response to literature format

        Args:
            json_response: 1H NMR prediction response

        Returns:
            str: Formatted 1H NMR prediction
        """
        formatted_signals = []

        for range_data in json_response["data"]["ranges"]:
            signal = range_data["signals"][0]
            delta = signal["delta"]
            integration = range_data["integration"]
            multiplicity = signal.get("multiplicity", "m")

            j_values = [
                round(j["coupling"], 1) for j in signal.get("js", []) if "coupling" in j
            ]

            signal_str = f"{delta:.2f}"
            if multiplicity != "m" and j_values:
                j_str = ", ".join(f"{j:.1f}" for j in j_values)
                signal_str += f" ({multiplicity}, J = {j_str} Hz, {integration}H)"
            else:
                signal_str += f" ({multiplicity}, {integration}H)"

            formatted_signals.append(signal_str)

        return f"Deltas {', '.join(formatted_signals)}."

    @staticmethod
    def format_ir(json_response: Dict) -> str:
        """Format IR response to literature format

        Args:
            json_response: IR prediction response

        Returns:
            str: Formatted IR prediction

        Raises:
            ValueError: If invalid method provided
        """
        wavenumbers = sorted(
            [
                round(mode["wavenumber"])
                for mode in json_response["modes"]
                if not mode["imaginary"]
            ],
            reverse=True,
        )
        wavenumbers = sorted(
            [
                round(mode["wavenumber"])
                for mode in json_response["modes"]
                if not mode["imaginary"] and mode["wavenumber"] >= 1500
            ],
            reverse=True,
        )
        return f"Wavenumbers (cm-1): {', '.join(map(str, wavenumbers))}"

    @staticmethod
    @backoff.on_exception(
        backoff.expo,
        (aiohttp.ClientError, asyncio.TimeoutError),
        max_tries=1,
        max_time=10,
        giveup=lambda e: isinstance(e, aiohttp.ClientResponseError)
        and e.status in {400, 401, 403, 404},
        jitter=backoff.full_jitter,
        base=2,
        logger=logger,
    )
    async def get_prediction_async(
        session: aiohttp.ClientSession,
        smiles: str,
        prediction_type: str,
        spectrum_type: str = "carbon",
        method: str = "GFN2xTB",
    ) -> Dict:
        """
        Get spectral prediction for a given SMILES string

        Args:
            session: aiohttp.ClientSession - Aiohttp client session
            smiles: str - SMILES string of the molecule
            prediction_type: str - Type of prediction ("nmr" or "ir")
            spectrum_type: str - Type of NMR spectrum (carbon or proton), only for NMR
            method: str - IR prediction method, only for IR

        Returns:
            Dict: Prediction response

        Raises:
            ValueError: If invalid prediction_type or spectrum_type provided
            aiohttp.ClientError: If the request fails after all retries
        """
        async with SpectraAPI._request_lock:
            current_time = asyncio.get_event_loop().time()
            time_since_last_request = current_time - SpectraAPI._last_request_time
            if time_since_last_request < SpectraAPI._min_request_interval:
                await asyncio.sleep(
                    SpectraAPI._min_request_interval - time_since_last_request
                )
            SpectraAPI._last_request_time = asyncio.get_event_loop().time()

        if prediction_type == "nmr":
            if spectrum_type not in SpectraAPI.VALID_SPECTRUM_TYPES:
                raise ValueError(
                    f"spectrum_type must be one of {SpectraAPI.VALID_SPECTRUM_TYPES}"
                )
            url = f"{SpectraAPI.NMR_BASE_URL}/{spectrum_type}"
            payload = {"smiles": smiles}
        elif prediction_type == "ir":
            url = SpectraAPI.IR_BASE_URL
            payload = {"smiles": smiles, "method": method}
        else:
            raise ValueError("prediction_type must be 'nmr' or 'ir'")

        async with session.post(
            url,
            headers={"Content-Type": "application/json"},
            json=payload,
            timeout=SpectraAPI._timeout,
        ) as response:
            try:
                response.raise_for_status()
                return await response.json()
            except Exception as e:
                logger.error(
                    f"{prediction_type.upper()} prediction failed - Status: {response.status}"
                )
                logger.error(f"Headers: {response.headers}")
                logger.error(f"Response body: {await response.text()}")
                raise e

    @classmethod
    async def get_all_predictions(cls, smiles: str) -> Dict[str, str]:
        """
        Get all spectral predictions for a molecule.
        If some predictions fail, still returns the successful ones.

        Args:
            smiles: SMILES string of the molecule

        Returns:
            Dict containing formatted spectral predictions
        """
        async with aiohttp.ClientSession() as session:
            logger.debug(f"Getting predictions for {smiles}")
            tasks = [
                cls.get_prediction_async(session, smiles, "nmr", "carbon"),
                cls.get_prediction_async(session, smiles, "nmr", "proton"),
                cls.get_prediction_async(session, smiles, "ir"),
            ]
            results = await asyncio.gather(*tasks, return_exceptions=True)

            predictions = {}

            # Process carbon NMR result
            if not isinstance(results[0], Exception):
                predictions["c13_nmr"] = SpectraAPI.format_c13_nmr(results[0])
            else:
                logger.error(f"Failed to retrieve C13 NMR: {str(results[0])}")
                predictions["c13_nmr"] = "C13 NMR prediction failed"

            # Process proton NMR result
            if not isinstance(results[1], Exception):
                predictions["h_nmr"] = SpectraAPI.format_h_nmr(results[1])
            else:
                logger.error(f"Failed to retrieve H NMR: {str(results[1])}")
                predictions["h_nmr"] = "H NMR prediction failed"

            # Process IR result
            if not isinstance(results[2], Exception):
                predictions["ir"] = SpectraAPI.format_ir(results[2])
            else:
                logger.error(f"Failed to retrieve IR: {str(results[2])}")
                predictions["ir"] = "IR prediction failed"

            return predictions

    @classmethod
    async def get_c13_nmr_prediction(cls, smiles: str) -> str:
        """
        Get only the C13 NMR spectral prediction for a molecule.

        Args:
            smiles: SMILES string of the molecule

        Returns:
            Formatted C13 NMR prediction or error message
        """
        async with aiohttp.ClientSession() as session:
            try:
                result = await cls.get_prediction_async(
                    session, smiles, "nmr", "carbon"
                )
                return SpectraAPI.format_c13_nmr(result)
            except Exception as e:
                logger.error(f"Failed to retrieve C13 NMR: {str(e)}")
                return "C13 NMR prediction failed"

    @classmethod
    async def get_h_nmr_prediction(cls, smiles: str) -> str:
        """
        Get only the H NMR spectral prediction for a molecule.

        Args:
            smiles: SMILES string of the molecule

        Returns:
            Formatted H NMR prediction or error message
        """
        async with aiohttp.ClientSession() as session:
            try:
                result = await cls.get_prediction_async(
                    session, smiles, "nmr", "proton"
                )
                return SpectraAPI.format_h_nmr(result)
            except Exception as e:
                logger.error(f"Failed to retrieve H NMR: {str(e)}")
                return "H NMR prediction failed"

    @classmethod
    async def get_ir_prediction(cls, smiles: str) -> str:
        """
        Get only the IR spectral prediction for a molecule.

        Args:
            smiles: SMILES string of the molecule

        Returns:
            Formatted IR prediction or error message
        """
        async with aiohttp.ClientSession() as session:
            try:
                result = await cls.get_prediction_async(session, smiles, "ir")
                return SpectraAPI.format_ir(result)
            except Exception as e:
                logger.error(f"Failed to retrieve IR: {str(e)}")
                return "IR prediction failed"


if __name__ == "__main__":
    import asyncio

    async def main():
        smiles = "CCO"  # Example SMILES for ethanol
        predictions = await SpectraAPI.get_ir_prediction(smiles)
        predictions = await SpectraAPI.get_c13_nmr_prediction(smiles)
        print("H NMR Prediction:")
        print(predictions)

    asyncio.run(main())
