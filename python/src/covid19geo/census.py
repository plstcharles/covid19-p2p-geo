import logging
import os
import typing
import zipfile

import h5py
import numpy as np
import pandas as pd
import tqdm

import covid19geo.utils

logger = logging.getLogger(__name__)

census_geo_code = "GEO_CODE (POR)"
"""Name of the census column used to identify the geographic level's code."""

census_da_member_id = "Member ID: Profile of Dissemination Areas (2247)"
"""Name of the census column used to identify the geographic level's attribute category."""

census_da_attrib_value = "Dim: Sex (3): Member ID: [1]: Total - Sex"
"""Name of the census column used to identify the geographic level's attribute value."""

census_cols_of_interest = {
    census_geo_code: str,
    "GEO_LEVEL": np.int32,
    census_da_member_id: np.int32,
    census_da_attrib_value: np.float32,
}
"""Mapping of census names to data types used for parsing and exportation."""

census_member_ids_of_interest = {  # these should hopefully rarely have NaNs
    1: "pop",  # Population, 2016
    4: "dwellings",  # Total private dwellings
    7: "area",  # Land area in square kilometres
}
"""Mapping of census attribute row indices to name that will be extracted from the CSV files."""

census_2016_csv_da_full_url = "https://www12.statcan.gc.ca/census-recensement/2016/dp-pd/prof/details/download-telecharger/comp/GetFile.cfm?Lang=E&FILETYPE=CSV&GEONO=044"
"""URL to the country-wide (1.7GB) 2016 census CSV data, split into dissemination areas."""

census_2016_csv_da_file_name = "98-401-X2016044_eng_CSV.zip"
"""File name for the country-wide (1.7GB) 2016 census CSV data."""

census_2016_csv_da_file_md5 = "1d4d14b0207858f74b4922580ba70687"
"""MD5 checksum value for the country-wide (1.7GB) 2016 census CSV data."""

census_2016_boundaries_url = "http://www12.statcan.gc.ca/census-recensement/2011/geo/bound-limit/files-fichiers/2016/lda_000a16a_e.zip"
"""URL to the dissemination area boundary shapefile (from 2011, but compatible with 2016 census)."""

census_2016_boundaries_file_name = "lada000b16a_e.zip"
"""File name for the 2016 dissemination area boundary shapefile."""

census_2016_boundaries_md5 = "f4c1d262ec498c469cd7603f20eea44f"
"""MD5 checksum value for the 2016 dissemination area boundary shapefile."""


class CensusDataExtractor:
    """Extracts useful census data and exports it in a variety of HDF5 archives.

    This interface offers three alternatives for exportation. It can:
      - export all data to a single HDF5 file (``export_full_hdf5``);
      - export high-level data (census divisions envelopes) to a single HDF5 file (``export_cd_hdf5``);
      - export low-level data (dissemination areas) to multiple HDF5s based on census divisions (``export_da_hdf5s``).
    The former is useful for debugging, visualization, and general manipulation. The latter
    two provide a hierarchical way to separate and query data dynamically on mobile devices.

    See the constructor's docstring for information on the input files from which the data
    is extracted.
    """

    def __init__(
            self,
            census_records_dir_path: typing.AnyStr,
            census_boundaries_dir_path: typing.AnyStr,
    ):
        """Parses and prepares census records & boundaries data for HDF5 exportation.

        Args:
            census_records_dir_path: the path to the directory containing CSV records of the
                2016 census split into dissemination areas. This data can be downloaded from:
                https://www12.statcan.gc.ca/census-recensement/2016/dp-pd/prof/details/download-telecharger/comp/page_dl-tc.cfm?Lang=E
                (pick the ``Canada, provinces, territories, census divisions (CDs), census
                subdivisions (CSDs) and dissemination areas (DAs)`` category)
            census_boundaries_dir_path: the path to the directory containing an Esri shapefile
                of the 2016 census dissemination areas. It can be downloaded from:
                https://www12.statcan.gc.ca/census-recensement/2011/geo/bound-limit/bound-limit-2016-eng.cfm
                (pick ``English`` for the language, ``ArcGIS (.shp)`` for the format,
                ``Dissemination areas``"`` for the geographic area, ``Digital boundary file``
                for the type, and click continue at the bottom for the download link)
        """
        assert os.path.isdir(census_boundaries_dir_path), \
            f"invalid census boundaries directory path: {census_boundaries_dir_path}"
        self.dauid_map, self.cduid_map, self.ptuid_map = \
            self._parse_geometries(census_boundaries_dir_path)
        assert os.path.isdir(census_records_dir_path), \
            f"invalid census record directory path: {census_records_dir_path}"
        dauid_stats_map = self._parse_records(
            census_records_dir_path,
            expected_da_count=len(self.dauid_map),
        )
        assert all([k in self.dauid_map for k in dauid_stats_map]), \
            "unexpected mis-overlap between da stats & geometries"
        for dauid, stats in dauid_stats_map.items():
            for k, v in stats.items():
                self.dauid_map[dauid][k] = v

    def export_full_hdf5(
            self,
            output_hdf5_path: typing.AnyStr,
            compress_wkb_lz4: bool = False,
    ):
        """Exports all provided data (and geoboundaries) to a single HDF5 file.

        Args:
            output_hdf5_path: filepath of the HDF5 archive that should be created.
            compress_wkb_lz4: toggles LZ4 compression on/off inside the HDF5 archive.
        """
        self._export_hdf5(self.dauid_map, self.cduid_map, self.ptuid_map,
                          output_hdf5_path, compress_wkb_lz4)

    def export_cd_hdf5(
            self,
            output_hdf5_path: typing.AnyStr,
    ):
        """Exports high-level regions ("census divisions envelopes") to an HDF5 file.

        The generated file should be very lightweight compared to the 'full' or
        'low-level' exports. Its purpose is to help determine which low-level regions
        to query in a hierarchical parsing approach.

        Args:
            output_hdf5_path: filepath of the HDF5 archive that should be created.
        """
        logger.debug(f"exporting cd data to: {output_hdf5_path}")
        with h5py.File(output_hdf5_path, "w") as fd:
            fd.attrs["cd_count"] = len(self.cduid_map)
            fd.attrs["stat_cols"] = list(census_member_ids_of_interest.values())
            fd.create_dataset(
                "cduid",
                data=np.asarray([np.uint32(id) for id in self.cduid_map.keys()]),
            )
            stats = np.asarray([
                [np.float32(sum([da[mid] for da in das]))
                 for mid in census_member_ids_of_interest.values()]
                for das in self.cduid_map.values()])
            fd.create_dataset("cd_stats", data=stats)
            cd_envelopes = []
            import ogr  # we import gdal components just before use in case someone can't install
            for cd_children in self.cduid_map.values():
                cd_geom = ogr.Geometry(ogr.wkbMultiPolygon)
                for cd_child in cd_children:
                    cd_geom.AddGeometry(ogr.CreateGeometryFromWkb(cd_child["wkb"]))
                cd_envelope = cd_geom.GetEnvelope()
                cd_envelopes.append([np.float32(v) for v in cd_envelope])
            fd.create_dataset("cd_envelope", data=np.asarray(cd_envelopes))

    def export_da_hdf5s(
            self,
            output_hdf5_directory: typing.AnyStr,
            compress_wkb_lz4: bool = False,
    ):
        """Exports low-level regions ("dissemination areas") to various HDF5 files.

        The generated files should individially be fairly lightweight compared to the
        'full' export. One file will be created for each census division. All files
        will be created in the same directory and using the census division geo ID (cduid)
        for naming.

        Args:
            output_hdf5_directory: directory path where the HDF5 files should be created.
                If the directory does not exist, it will be created.
            compress_wkb_lz4: toggles LZ4 compression on/off inside the HDF5 archive.
        """
        os.makedirs(output_hdf5_directory, exist_ok=True)
        for cduid, cd_children in self.cduid_map.items():
            output_hdf5_path = os.path.join(output_hdf5_directory, f"{cduid}.hdf5")
            dauid_map = {cd_child["dauid"]: cd_child for cd_child in cd_children}
            cduid_map = {cduid: cd_children}
            ptuid_map = {ptuid: [cd_child for cd_child in cd_children
                                 if cd_child["ptuid"] == ptuid]
                         for ptuid in np.unique([cd_child["ptuid"] for cd_child in cd_children])}
            self._export_hdf5(dauid_map, cduid_map, ptuid_map, output_hdf5_path, compress_wkb_lz4)

    def plot_da_stats(self):
        """Computes, prints & plots useful statistics about dissemination areas."""
        import plotly.express as px
        df = pd.DataFrame(list(self.dauid_map.values()))
        fig = px.density_contour(
            df,
            x="area",
            y="pop",
            marginal_x="histogram",
            marginal_y="histogram",
            range_x=[0, 2000],
            range_y=[0, 2000],
        )
        fig.data[0].update(contours_coloring="fill", contours_showlabels=True)
        fig.show()
        below_100_pops = [m for m in self.dauid_map.values() if m["pop"] <= 100]
        below_200_pops = [m for m in self.dauid_map.values() if m["pop"] <= 200]
        print(f"DA count: {len(self.dauid_map)}")
        print(f"\tDAs below 100 pop: {len(below_100_pops)} (population sum = {sum([int(m['pop']) for m in below_100_pops])})")
        print(f"\tDAs below 200 pop: {len(below_200_pops)} (population sum = {sum([int(m['pop']) for m in below_200_pops])})")


    @staticmethod
    def _export_hdf5(
            dauid_map: typing.Dict[str, typing.Dict],
            cduid_map: typing.Dict[str, typing.List[typing.Dict]],
            ptuid_map: typing.Dict[str, typing.List[typing.Dict]],
            output_hdf5_path: typing.AnyStr,
            compress_wkb_lz4: bool = False,
    ):
        logger.debug(f"exporting data to: {output_hdf5_path}")
        with h5py.File(output_hdf5_path, "w") as fd:
            CensusDataExtractor._export_hdf5_header_attribs(dauid_map, cduid_map, ptuid_map, fd)
            fd.create_dataset(
                "ptuid",
                data=np.asarray([np.uint32(v["ptuid"]) for v in dauid_map.values()]),
            )
            fd.create_dataset(
                "cduid",
                data=np.asarray([np.uint32(v["cduid"]) for v in dauid_map.values()]),
            )
            fd.create_dataset(
                "dauid",
                data=np.asarray([np.uint32(v["dauid"]) for v in dauid_map.values()]),
            )
            stats = np.asarray([
                [np.float32(v[mid]) for mid in census_member_ids_of_interest.values()]
                for v in dauid_map.values()])
            fd.create_dataset("da_stats", data=stats)
            wkb_dataset = fd.create_dataset(
                "da_wkb",
                shape=(len(dauid_map),),
                maxshape=(len(dauid_map),),
                dtype=h5py.special_dtype(vlen=np.uint8),
            )
            assert not compress_wkb_lz4, "missing lz4 compression impl @@@"
            total_wkb_bytes = 0
            for da_idx, da_dict in tqdm.tqdm(enumerate(dauid_map.values()),
                                             total=len(dauid_map)):
                wkb_dataset[da_idx] = np.frombuffer(da_dict["wkb"], dtype=np.uint8)
                total_wkb_bytes += len(da_dict["wkb"])
            wkb_dataset.attrs["byte_count"] = total_wkb_bytes

    @staticmethod
    def _export_hdf5_header_attribs(dauid_map, cduid_map, ptuid_map, fd):
        fd.attrs["da_count"] = len(dauid_map)
        fd.attrs["cd_count"] = len(cduid_map)
        fd.attrs["pt_count"] = len(ptuid_map)
        fd.attrs["stat_cols"] = list(census_member_ids_of_interest.values())

    @staticmethod
    def _parse_records(
            census_records_dir_path: typing.AnyStr,
            expected_da_count: typing.Optional[int] = None,
    ):
        census_records_file_path = covid19geo.utils.get_unique_file_by_extension(
            directory_path=census_records_dir_path,
            extension="_CSV_data.csv",
        )
        logger.debug(f"parsing records CSV file: {census_records_file_path}")
        chunksize = 2 ** 16
        reader = pd.read_csv(
            filepath_or_buffer=census_records_file_path,
            chunksize=chunksize,
            usecols=census_cols_of_interest.keys(),
            dtype=census_cols_of_interest,
            na_values=["..", "...", "F", "x"],
        )
        dauid_map = {}
        if expected_da_count is not None:
            expected_rows = expected_da_count * 2247  # nb attribs per DA
            expected_iters = (expected_rows // chunksize) + 1
            reader = tqdm.tqdm(reader, desc="parsing CSV rows", total=expected_iters)
        else:
            reader = tqdm.tqdm(reader, desc="parsing CSV rows")
        for rows in reader:
            rows = rows.loc[rows[census_da_member_id].isin(census_member_ids_of_interest)]
            rows = rows.loc[rows["GEO_LEVEL"] == 4]  # keep only da-level information
            if not rows.empty:
                for idx, row in rows.iterrows():
                    if row[census_geo_code] not in dauid_map:
                        dauid_map[row[census_geo_code]] = {}
                    key = census_member_ids_of_interest[row[census_da_member_id]]
                    dauid_map[row[census_geo_code]][key] = row[census_da_attrib_value]
        return dauid_map

    @staticmethod
    def _parse_geometries(
            census_boundaries_dir_path: typing.AnyStr,
            project_epsg3347_to_epsg4326: bool = True,
    ):
        census_boundaries_file_path = covid19geo.utils.get_unique_file_by_extension(
            directory_path=census_boundaries_dir_path,
            extension=".shp",
        )
        logger.debug(f"parsing boundaries shapefile: {census_boundaries_file_path}")
        import ogr  # we import gdal components just before use in case someone can't install
        ogr_driver = ogr.GetDriverByName("ESRI Shapefile")
        boundaries_fd = ogr_driver.Open(census_boundaries_file_path, 0)
        boundaries_layer = boundaries_fd.GetLayer()
        boundaries_count = boundaries_layer.GetFeatureCount()
        logger.debug(f"will load {boundaries_count} features")
        boundaries_layer_defs = boundaries_layer.GetLayerDefn()
        boundaries_attrib_count = boundaries_layer_defs.GetFieldCount()
        attribs_full_str = f"each feature has {boundaries_attrib_count} attributes:"
        found_dauid = False
        for i in range(boundaries_layer_defs.GetFieldCount()):
            field_name = boundaries_layer_defs.GetFieldDefn(i).GetName()
            field_type_code = boundaries_layer_defs.GetFieldDefn(i).GetType()
            field_type = boundaries_layer_defs.GetFieldDefn(i).GetFieldTypeName(field_type_code)
            field_width = boundaries_layer_defs.GetFieldDefn(i).GetWidth()
            field_precision = boundaries_layer_defs.GetFieldDefn(i).GetPrecision()
            field_str = \
                f"{field_name}: {field_type} (width={field_width}, precision={field_precision})"
            attribs_full_str += f"\n\t{field_str}"
            if field_name == "DAUID":
                assert field_width == 8 and field_type == "String", \
                    "unexpected DAUID field width/type"
                found_dauid = True
        logger.debug(attribs_full_str)
        assert found_dauid, "features missing expected 'DAUID' field"
        # will load all dissemination areas into memory while keeping track of 3-level hierarchy
        dauid_map, cduid_map, ptuid_map = {}, {}, {}
        transform = None
        if project_epsg3347_to_epsg4326:
            import osr  # we import gdal components just before use in case someone can't install
            source_ref, target_ref = osr.SpatialReference(), osr.SpatialReference()
            source_ref.ImportFromEPSG(3347)
            assert source_ref.IsSame(boundaries_layer.GetSpatialRef())
            target_ref.ImportFromEPSG(4326)
            transform = osr.CoordinateTransformation(source_ref, target_ref)
        for feature in tqdm.tqdm(boundaries_layer, desc="parsing features"):
            dauid = feature.GetField("DAUID")
            cduid, ptuid = dauid[0:4], dauid[0:2]
            geom = feature.GetGeometryRef()
            if transform is not None:
                geom.Transform(transform)
            geometry_ref = dict(
                ptuid=ptuid,
                cduid=cduid,
                dauid=dauid,
                wkb=geom.ExportToWkb(),
            )
            assert dauid not in dauid_map, "unexpected duplicate DA geometry"
            dauid_map[dauid] = geometry_ref
            if cduid not in cduid_map:
                cduid_map[cduid] = []
            cduid_map[cduid].append(geometry_ref)
            if ptuid not in ptuid_map:
                ptuid_map[ptuid] = []
            ptuid_map[ptuid].append(geometry_ref)
        return dauid_map, cduid_map, ptuid_map


def validate_exported_hdf5(hdf5_path: typing.AnyStr):
    """Runs simple shape/count tests on exported HDF5 dataset attributes."""
    logger.debug(f"parse-testing hdf5 at: {hdf5_path}")
    with h5py.File(hdf5_path, "r") as fd:
        assert np.array_equal(fd.attrs["stat_cols"], list(census_member_ids_of_interest.values()))
        assert fd.attrs["da_count"] >= fd.attrs["cd_count"]
        assert fd.attrs["cd_count"] >= fd.attrs["pt_count"]
        assert fd["ptuid"].shape[0] == fd.attrs["da_count"]
        assert fd["cduid"].shape[0] == fd.attrs["da_count"]
        assert fd["dauid"].shape[0] == fd.attrs["da_count"]
        assert fd["da_stats"].shape[0] == fd.attrs["da_count"]
        assert fd["da_stats"].shape[1] == len(census_member_ids_of_interest)
        assert fd["da_wkb"].shape[0] == fd.attrs["da_count"]


def _download_and_extract_census_records(
        download_dir: typing.AnyStr,
        checksum: bool,
):
    census_records_path = os.path.join(download_dir, census_2016_csv_da_file_name)
    if not os.path.exists(census_records_path):
        logger.debug("downloading census 2016 csv da records (1.7GB total)")
        covid19geo.utils.download_file_with_progress_bar(
            file_url=census_2016_csv_da_full_url,
            file_path=census_records_path,
        )
    if checksum:
        census_records_md5 = covid19geo.utils.compute_md5_hash(census_records_path)
        assert census_records_md5 == census_2016_csv_da_file_md5, \
            "census 2016 CSV records failed checksum:\n" \
            f"\t{census_records_md5}\n\t\tvs\n\t{census_2016_csv_da_file_md5}"
    census_records_dir_name = os.path.splitext(census_2016_csv_da_file_name)[0]
    census_records_dir = os.path.join(download_dir, census_records_dir_name)
    if not os.path.exists(census_records_dir):
        logger.debug("extracting census 2016 csv da records (14.6GB total)")
        with zipfile.ZipFile(census_records_path, "r") as fd:
            fd.extractall(census_records_dir)
    return census_records_dir


def _download_and_extract_census_boundaries(
        download_dir: typing.AnyStr,
        checksum: bool,
):
    census_boundaries_path = os.path.join(download_dir, census_2016_boundaries_file_name)
    if not os.path.exists(census_boundaries_path):
        logger.debug("downloading census 2016 boundaries (67MB total)")
        covid19geo.utils.download_file_with_progress_bar(
            file_url=census_2016_boundaries_url,
            file_path=census_boundaries_path,
        )
    if checksum:
        census_boundaries_md5 = covid19geo.utils.compute_md5_hash(census_boundaries_path)
        assert census_boundaries_md5 == census_2016_boundaries_md5, \
            "census 2016 boundaries failed checksum:\n" \
            f"\t{census_boundaries_md5}\n\t\tvs\n\t{census_2016_boundaries_md5}"
    census_boundaries_dir_name = os.path.splitext(census_2016_boundaries_file_name)[0]
    census_boundaries_dir = os.path.join(download_dir, census_boundaries_dir_name)
    if not os.path.exists(census_boundaries_dir):
        logger.debug("extracting census 2016 boundaries (183MB total)")
        with zipfile.ZipFile(census_boundaries_path, "r") as fd:
            fd.extractall(census_boundaries_dir)
    return census_boundaries_dir


def download_and_extract_census_data(
        download_dir: typing.AnyStr,
        checksum: bool = True,
):
    """Downloads and extracts census records and boundary files for HDF5 exportation."""
    records_dir = _download_and_extract_census_records(download_dir, checksum)
    boundaries_dir = _download_and_extract_census_boundaries(download_dir, checksum)
    return records_dir, boundaries_dir
