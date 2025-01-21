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
    config = configparser.ConfigParser()
    config.add_section("PATHS")
    config.set("PATHS", "group_data", "mock_group_data.csv")
    config.set("PATHS", "data_folder", "mock_data_folder")
    return config


@pytest.mark.parametrize(
    "mock_data, expected_columns",
    [
        (be.pl.DataFrame({" description": ["controle", "behandeling"], "barcode": ["1","2"]}), [" description", "barcode", "group_and_n"]),
        (be.pl.DataFrame({" description": ["Other", "Other_test"], "barcode": ["1","2"]}), [" description", "barcode", "group_and_n"])
    ]
)
def test_process_group(fake_config, mock_data, expected_columns):
    with patch("polars.read_csv", return_value = mock_data):
        df = be.process_groups(fake_config)
        for col in expected_columns:
            assert col in df.columns


def test_process_group_bad_column(fake_config):
    mock_data = be.pl.DataFrame({"beshrijving": ["controle", "behandeling"], "barcode": ["1","2"]})
    with patch("polars.read_csv", return_value = mock_data):
        with pytest.raises(be.pl.exceptions.ColumnNotFoundError):
            df = be.process_groups(fake_config)


def test_process_group_false_path(fake_config):
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
                         ["n methylations", "group_name"])
    ]
)
def test_count_data(mock_data,expected_columns):
    df = be.count_methylation_data(mock_data)
    for column in expected_columns:
        assert(column in df.columns)


@pytest.mark.asyncio
async def test_plotting_empty_df():
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
    "chr_ls, df",  # Fix: Parameters should be a single string with comma-separated names
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
    assert isinstance(be.filter_chr(chr_ls, df), be.pl.DataFrame)