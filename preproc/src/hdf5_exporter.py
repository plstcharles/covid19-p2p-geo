import logging
import os
import typing

import numpy as np
import ogr
import pandas as pd
import tqdm

import utils

logger = logging.getLogger(__name__)

census_geo_code = "GEO_CODE (POR)"
census_da_member_id = "Member ID: Profile of Dissemination Areas (2247)"
census_da_attrib_value = "Dim: Sex (3): Member ID: [1]: Total - Sex"

census_cols_of_interest = {
    census_geo_code: str,
    "GEO_LEVEL": np.int32,
    census_da_member_id: np.int32,
    census_da_attrib_value: np.float32,
}

census_member_ids_of_interest = {  # these should hopefully rarely have NaNs
    1: "pop",  # Population, 2016
    4: "dwellings",  # Total private dwellings
    5: "res-dwellings",  # Private dwellings occupied by usual residents
    7: "area",  # Land area in square kilometres
}


class CensusCompactor:
    """Compresses useful census data inside a HDF5 archive."""

    def __init__(
            self,
            census_records_dir_path: typing.AnyStr,
            census_boundaries_dir_path: typing.AnyStr,
    ):
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
        # @@@ todo next: merge stats + geometries
        # ... then, export to hdf5


    @staticmethod
    def _parse_records(
            census_records_dir_path: typing.AnyStr,
            expected_da_count: typing.Optional[int] = None,
    ):
        census_records_file_path = utils.get_unique_file_by_extension(
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
            census_boundaries_dir_path: typing.AnyStr
    ):
        census_boundaries_file_path = utils.get_unique_file_by_extension(
            directory_path=census_boundaries_dir_path,
            extension=".shp",
        )
        logger.debug(f"parsing boundaries shapefile: {census_boundaries_file_path}")
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
        for feature in tqdm.tqdm(boundaries_layer, desc="parsing features"):
            dauid = feature.GetField("DAUID")
            cduid, ptuid = dauid[0:4], dauid[0:2]
            geometry_ref = dict(
                ptuid=ptuid,
                cduid=cduid,
                dauid=dauid,
                wkb=feature.GetGeometryRef().ExportToWkb(),
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



if __name__ == "__main__":
    # @@@@ TODO: CONVERT TO PROPER TEST
    logging.basicConfig()
    logging.getLogger().setLevel(logging.NOTSET)
    root_data_path = "/home/perf6/dev/covid19-p2p-geo/data/"
    _census_records_dir_path = \
        os.path.join(root_data_path, "low_level_data_98-401-X2016044_eng_CSV")
    _census_boundaries_dir_path = \
        os.path.join(root_data_path, "boundaries_digital_dissemination_areas_lada000b16a_e")
    CensusCompactor(
        census_records_dir_path=_census_records_dir_path,
        census_boundaries_dir_path=_census_boundaries_dir_path,
    )