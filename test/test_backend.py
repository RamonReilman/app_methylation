"""
test_backend.py
Author: Ramon Reilman
Version: 1.0
Year: BFV2

Usage:
This script will test certain functions from the backend

run:
python3 -m pytest
"""

import pytest
import configparser
import asyncio
import pytest_asyncio
from unittest.mock import patch
from src import backend as be

config = be.parse_config()
main_data = be.read_data(config=config)

@pytest.fixture
def fake_config():
    """ Will create a fake config file
    This function creates a fake config file that backend functions need to work.
    
    Parameters
    ----------
    None

    Returns
    -------
    config.configparser
        Fake config file

    """
    config = configparser.ConfigParser()
    config.add_section("PATHS")
    config.set("PATHS", "group_data", "mock_group_data.csv")
    config.set("PATHS", "data_folder", "mock_data_folder")
    return config


@pytest.mark.parametrize(
    "mock_data, expected_columns",
    [
        (be.pl.DataFrame({" description": ["controle", "behandeling"], "barcode": ["1","2"]}),
         [" description", "barcode", "group_and_n"]),
        (be.pl.DataFrame({" description": ["Other", "Other_test"], "barcode": ["1","2"]}),
         [" description", "barcode", "group_and_n"])
    ]
)
def test_process_group(fake_config, mock_data, expected_columns):
    """Test the process group function
    This test will test the process group function from the backend
    
    Parameters
    ----------
    fake_config
        Fake config file that the function requires
    mock_data
        fake data to test the columns
    expected_columns
        Columns that should be returned

    """
    with patch("polars.read_csv", return_value = mock_data):
        df = be.process_groups(fake_config)
        for col in expected_columns:
            assert col in df.columns


def test_process_group_bad_column(fake_config):
    """Test the process group function
    This test will test the process group function from the backend
    This test tests what happens when a bad header is given
    should raise a pl.exceptions.ColumnNotFoundError
    
    Parameters
    ----------
    fake_config
        Fake config file that the function requires
    """
    mock_data = be.pl.DataFrame({"beshrijving": ["controle", "behandeling"], "barcode": ["1","2"]})
    with patch("polars.read_csv", return_value = mock_data):
        with pytest.raises(be.pl.exceptions.ColumnNotFoundError):
            df = be.process_groups(fake_config)


def test_process_group_false_path(fake_config):
    """Test the process group function
    This test will test the process group function from the backend
    This test tests what happens when an incorrect path is given
    Should raise a filenotfound error
    
    Parameters
    ----------
    fake_config
        Fake config file that the function requires
    """
    with pytest.raises(FileNotFoundError):
        df = be.process_groups(fake_config)

@pytest.mark.parametrize(
    "mock_data, expected_columns",
    [
        (be.pl.DataFrame({"chr": ["chr1", "chr2"], "start": [1,2],
                          "end": [3,5], "frac":[1,1],
                          "valid": [1,1], "group_name": ["group_1", "group_2"]}), 
                         ["n methylations", "group_name"]),
        (be.pl.DataFrame({"chr": ["safs", "fasa"], "start": [124,132],
                          "end": [52135,4234], "frac":[13,32321321],
                          "valid": [15332,152532], "group_name": ["group_1", "group_2"]}), 
                         ["n methylations", "group_name"]),
        (be.pl.DataFrame({"chr": ["safs", "fasa"], "start": ["1", "3"],
                          "end": ["1", "2"], "frac":[13,32321321],
                          "valid": [15332,152532], "group_name": ["group_1", "group_2"]}), 
                         ["n methylations", "group_name"])
    ]
)
def test_count_data(mock_data,expected_columns):
    """Test the count methylation data function
    This test will test the process group function from the backend
    this tests if the function will work with different values
    Like if start is defined as a string
    
    Parameters
    ----------
    mock_data
        fake data to test the columns
    expected_columns
        Columns that should be returned
    """
    df = be.count_methylation_data(mock_data)
    for column in expected_columns:
        assert(column in df.columns)


@pytest.mark.asyncio
async def test_plotting_empty_df():
    """Tests the plotting functionality
    Tests what happens when an empty df is given
    """
    empty_df = be.pl.DataFrame(
        {"chr":[],
         "start":[],
         "end":[],
         "frac":[],
         "valid":[],
         "group_name":[]})
    plots_task = be.plot_plots(empty_df, True)
    plots = await plots_task

    assert isinstance(plots, be.pn.widgets.indicators.LoadingSpinner)


@pytest.mark.asyncio
async def test_plotting_df():
    """Tests the plotting functionality
    Tests what happens when an filled df is given
    """
    df = be.pl.DataFrame(
        {"chr":["chr1"],
         "start":[2],
         "end":[5],
         "frac":[1],
         "valid":[1],
         "group_name":["Test"]})

    plots_task = be.plot_plots(df, True)
    plots = await plots_task

    assert isinstance(plots, list)



@pytest.mark.parametrize(
    "chr_ls, df",
    [
        pytest.param(
            ["chr1", "chr2"], main_data
        ),
        pytest.param(
            ["dhhds", "djdsjd"], main_data
        ),
        pytest.param(
            [], main_data
        ),
    ]
)
def test_filter_chr(chr_ls, df):
    """Test the filter chr function
    what happens when weird chromosomes are given, or an empty list
    
    Parameters
    ----------
    chr_ls: list
        list with chromosomes wanted
    df: pl.dataframe
        Data to filter on
    """
    assert isinstance(be.filter_chr(chr_ls, df), be.pl.DataFrame)
