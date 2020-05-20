import logging
import os
import pathlib
import typing

import covid19geo.census

logger = logging.getLogger(__name__)


def main(
        download_dir: typing.Optional[typing.AnyStr] = None,
        output_dir: typing.Optional[typing.AnyStr] = None,
        checksum: bool = True,
        plot_da_stats: bool = False,
):
    """Downloads and extracts useful data from 2016 census records/boundaries.

    By default, the data will be downloaded & exported in the project's ``data``
    directory, two levels up from this file's path.
    """
    logging.basicConfig()
    logging.getLogger().setLevel(logging.NOTSET)
    if download_dir is None:
        # by default, we will go up two directories to find the C++ project root;
        # this might not be the correct default place to download stuff if you
        # installed the Python package directly. In that case, make sure you provide
        # a default argument!
        download_dir = pathlib.Path(__file__).parent.absolute().parents[2] / "data"
    logger.debug(f"path where data might be downloaded/unpacked: {download_dir}")
    os.makedirs(download_dir, exist_ok=True)
    records_dir_path, boundaries_dir_path = \
        covid19geo.census.download_and_extract_census_data(download_dir, checksum)
    extractor = covid19geo.census.CensusDataExtractor(
        census_records_dir_path=records_dir_path,
        census_boundaries_dir_path=boundaries_dir_path,
    )
    if plot_da_stats:
        extractor.plot_da_stats()
    if output_dir is None:
        output_dir = download_dir
    output_full_hdf5_path = os.path.join(output_dir, "full.hdf5")
    extractor.export_full_hdf5(output_full_hdf5_path)
    covid19geo.census.validate_exported_hdf5(output_full_hdf5_path)
    output_cd_hdf5_path = os.path.join(output_dir, "cd.hdf5")
    extractor.export_cd_hdf5(output_cd_hdf5_path)
    output_da_hdf5_dir_path = os.path.join(output_dir, "divisions")
    extractor.export_da_hdf5s(output_da_hdf5_dir_path)


if __name__ == "__main__":
    main()
