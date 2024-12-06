
import pytest
import asyncio
import aiohttp
from unittest.mock import patch, AsyncMock

from chemenv.tools.converters import (
    Name2Smiles,
    Smiles2Name,
    smiles_to_selfies,
    smiles_to_deepsmiles,
    smiles_to_canoncial,
    smiles_to_inchi,
    smiles_to_safe,
)

def test_smiles_to_selfies():
    assert smiles_to_selfies("CCO") == "[C][C][O]"
    assert smiles_to_selfies("CC") == "[C][C]"

def test_smiles_to_deepsmiles():
    assert smiles_to_deepsmiles("CCO") == "CCO"
    assert smiles_to_deepsmiles("CC") == "CC"

def test_smiles_to_canoncial():
    assert smiles_to_canoncial("CCO") == "CCO"
    assert smiles_to_canoncial("CC") == "CC"

def test_smiles_to_inchi():
    assert smiles_to_inchi("CCO") == "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3"
    assert smiles_to_inchi("CC") == "InChI=1S/C2H6/c1-2/h1-2H3"

@pytest.mark.asyncio
async def test_name2smiles_failures():
    # Test failed API calls return None
    converter = Name2Smiles("ethanol")
    
    # Test non-200 response
    mock_response = AsyncMock()
    mock_response.status = 404
    
    with patch("aiohttp.ClientSession.get", return_value=AsyncMock(__aenter__=AsyncMock(return_value=mock_response))):
        result = await converter.opsin()
        assert result is None
        
        result = await converter.cactus()
        assert result is None
        
        result = await converter.pubchem() 
        assert result is None
        
        result = await converter.unknown()
        assert result is None

    # Test invalid JSON response
    mock_response.status = 200
    mock_response.json.side_effect = ValueError
    
    with patch("aiohttp.ClientSession.get", return_value=AsyncMock(__aenter__=AsyncMock(return_value=mock_response))):
        result = await converter.opsin()
        assert result is None
        
        result = await converter.pubchem()
        assert result is None

@pytest.mark.asyncio
async def test_name2smiles_retries():
    converter = Name2Smiles("ethanol")
    
    # Mock that fails twice then succeeds
    mock_response = AsyncMock()
    mock_response.status = 200
    mock_response.json.return_value = {"smiles": "CCO"}
    
    with patch("aiohttp.ClientSession.get", 
              side_effect=[aiohttp.ClientError, aiohttp.ClientError, 
                          AsyncMock(__aenter__=AsyncMock(return_value=mock_response))]):
        result = await converter.opsin()
        assert result == "CCO"
        
async def test_cactus_api(name2smiles):
    mock_response = AsyncMock()
    mock_response.status = 200
    mock_response.text.return_value = "CCO"
    
    with patch("aiohttp.ClientSession.get", return_value=AsyncMock(__aenter__=AsyncMock(return_value=mock_response))):
        result = await name2smiles.cactus()
        assert result == "CCO"

@pytest.mark.asyncio
async def test_pubchem_api(name2smiles):
    mock_response = AsyncMock()
    mock_response.status = 200
    mock_response.json.return_value = {
        "PropertyTable": {
            "Properties": [{"IsomericSMILES": "CCO"}]
        }
    }
    
    with patch("aiohttp.ClientSession.get", return_value=AsyncMock(__aenter__=AsyncMock(return_value=mock_response))):
        result = await name2smiles.pubchem()
        assert result == "CCO"

@pytest.mark.asyncio
async def test_unknown_api(name2smiles):
    mock_response = AsyncMock()
    mock_response.status = 200
    mock_response.text.return_value = "Message:CCO"
    
    with patch("aiohttp.ClientSession.get", return_value=AsyncMock(__aenter__=AsyncMock(return_value=mock_response))):
        result = await name2smiles.unknown()
        assert result == "CCO"

@pytest.mark.asyncio
async def test_api_errors(name2smiles):
    # Test timeout error
    with patch("aiohttp.ClientSession.get", side_effect=asyncio.TimeoutError):
        result = await name2smiles.opsin()
        assert result is None

    # Test client error
    with patch("aiohttp.ClientSession.get", side_effect=aiohttp.ClientError):
        result = await name2smiles.cactus()
        assert result is None

@pytest.mark.asyncio
async def test_get_smiles_parallel(name2smiles):
    mock_opsin = AsyncMock(return_value=None)
    mock_cactus = AsyncMock(return_value="CCO")
    mock_pubchem = AsyncMock(return_value=None)
    mock_unknown = AsyncMock(return_value=None)
    
    with patch.multiple(name2smiles, 
                        opsin=mock_opsin,
                        cactus=mock_cactus,
                        pubchem=mock_pubchem,
                        unknown=mock_unknown):
        result = await name2smiles.get_smiles()
        assert result == "CCO"
        assert isinstance(result, str)

@pytest.fixture
def smiles2name():
    return Smiles2Name("CCO")

@pytest.mark.asyncio
async def test_smiles2name_init():
    # Valid SMILES
    converter = Smiles2Name("CCO")
    assert converter.smiles == "CCO"
    
    # Invalid SMILES
    with pytest.raises(ValueError, match="Invalid SMILES"):
        Smiles2Name("invalid_smiles")

@pytest.mark.asyncio
async def test_smiles2name_pubchem(smiles2name):
    # Test successful response
    mock_response = AsyncMock()
    mock_response.status = 200
    mock_response.text.return_value = "ethanol"
    
    with patch("aiohttp.ClientSession.get", return_value=AsyncMock(__aenter__=AsyncMock(return_value=mock_response))):
        result = await smiles2name.pubchem()
        assert result == "ethanol"

    # Test failed response
    mock_response.status = 404
    with patch("aiohttp.ClientSession.get", return_value=AsyncMock(__aenter__=AsyncMock(return_value=mock_response))):
        result = await smiles2name.pubchem()
        assert result is None

@pytest.mark.asyncio
async def test_smiles2name_cactus(smiles2name):
    # Test successful response
    mock_response = AsyncMock()
    mock_response.status = 200
    mock_response.text.return_value = "ethanol"
    
    with patch("aiohttp.ClientSession.get", return_value=AsyncMock(__aenter__=AsyncMock(return_value=mock_response))):
        result = await smiles2name.cactus()
        assert result == "ethanol"

    # Test failed response
    mock_response.status = 404
    with patch("aiohttp.ClientSession.get", return_value=AsyncMock(__aenter__=AsyncMock(return_value=mock_response))):
        result = await smiles2name.cactus()
        assert result is None

@pytest.mark.asyncio
async def test_smiles2name_api_errors(smiles2name):
    # Test timeout error
    with patch("aiohttp.ClientSession.get", side_effect=asyncio.TimeoutError):
        result = await smiles2name.pubchem()
        assert result is None

    # Test client error
    with patch("aiohttp.ClientSession.get", side_effect=aiohttp.ClientError):
        result = await smiles2name.cactus()
        assert result is None

@pytest.mark.asyncio
async def test_smiles2name_get_name_parallel(smiles2name):
    # Test when first API succeeds
    mock_cactus = AsyncMock(return_value="ethanol")
    mock_pubchem = AsyncMock(return_value=None)
    
    with patch.multiple(smiles2name,
                        cactus=mock_cactus,
                        pubchem=mock_pubchem):
        result = await smiles2name.get_name()
        assert result == "ethanol"

    # Test when all APIs fail
    mock_cactus = AsyncMock(return_value=None)
    mock_pubchem = AsyncMock(return_value=None)
    
    with patch.multiple(smiles2name,
                        cactus=mock_cactus,
                        pubchem=mock_pubchem):
        with pytest.raises(ValueError, match="Could not find name"):
            await smiles2name.get_name()