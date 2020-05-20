import copy
import os
import pathlib
import typing

import shapely.geometry
import shapely.strtree
import shapely.wkb

import covid19geo.census
import covid19geo.utils
from covid19geo.utils import LatLonPairType, BoundingBoxType


class ZoneBase:
    """Builds a hierarchy of geographic regions for a given zone bounding box.

    This base class implements functions common to the real & dummy zone implementations.
    """
    def __init__(self):
        """Initializes common attribs for derived classes."""
        self.dauid_map: typing.Dict[int, typing.Dict] = {}
        self.cduid_map: typing.Dict[int, typing.List[typing.Dict]] = {}
        self.rtree: typing.Optional[shapely.strtree.STRtree] = None

    def get_areas_map(self) -> typing.Dict[int, typing.Dict]:
        """Returns the map of all loaded dissemination areas (with features, if available)."""
        return self.dauid_map

    def get_divisions_map(self) -> typing.Dict[int, typing.List[typing.Dict]]:
        """Returns the map of all census divisions, where each division contains a list of areas."""
        return self.cduid_map

    def get_coord_dauid(self, coords: LatLonPairType) -> int:
        """Returns the dissemination area UID (dauid) at the specified lat/lon coordinates.

        If no match is found, the function will return -1.
        """
        query_pt = shapely.geometry.Point(*coords)
        query_res = self.rtree.query(query_pt)
        for hit_poly in query_res:
            if hit_poly.contains(query_pt):
                assert hasattr(hit_poly, "_dauid"), "attrib should have been added..."
                return hit_poly._dauid
        return -1


class DummyZone(ZoneBase):
    """Builds a hierarchy of geographic regions for a given zone bounding box.

    This implementation will construct dummy census divisions & dissemination areas based
    on the requested division/area counts. Census divisions correspond to high-level regions
    and dissemination areas correspond to low-level regions inside the former. The 'zone' here
    is equivalent to a province/territory in terms of region encoding, but it can actually be
    much smaller than one in practice.

    This class is pretty much only destined for testing/debug purposes inside the simulator.
    """

    def __init__(
            self,
            bounding_box: BoundingBoxType,
            division_rows: int,
            division_cols: int,
            area_rows: int,
            area_cols: int,
            zone_prefix: int = 99,
            divisions_offset: int = 0,
            areas_offset: int = 0,
            fill_dummy_census_stats: bool = True,
    ):
        """Parses and prepares divisions/areas based on the provided params.

        Args:
            bounding_box: the bounding box of the zone for which to build divisions/areas.
            division_rows: the number of census division rows to slice the zone into.
            division_cols: the number of census division columns to slice the zone into.
            area_rows: the number of dissemination area rows to slice each division into.
            area_cols: the number of dissemination area columns to slice each division into.
            zone_prefix: the integer used in region UIDs that identifies the province/territory.
            divisions_offset: the integer offset to start counting division UIDs from.
            areas_offset: the integer offset to start counting area UIDs from.
        """
        assert 0 <= zone_prefix < 100, "invalid zone prefix (should be in [0, 99])"
        assert 0 < division_rows * division_cols and 0 < area_rows * area_cols, \
            "zone divisions/areas counts must not be null"
        assert division_rows * division_cols + divisions_offset <= 100, \
            f"division UIDs must all fit inside the [{divisions_offset}, 99] range"
        assert area_rows * area_cols + areas_offset <= 10000, \
            f"area UIDs must all fit inside the [{areas_offset}, 9999] range"
        super().__init__()
        self.zone_top_left = bounding_box[0]  # lat, lon
        self.zone_bottom_right = bounding_box[1]  # lat, lon
        self.divisions = division_rows, division_cols
        self.areas = area_rows, area_cols
        self.zone_prefix = zone_prefix
        self.dauid_map, self.cduid_map = self._slice_regions(divisions_offset, areas_offset)
        if fill_dummy_census_stats:
            for dauid, da in self.dauid_map.items():
                for census_stat in covid19geo.census.census_member_ids_of_interest.values():
                    assert census_stat not in da
                    da[census_stat] = None
        da_polygons = [da["polygon"] for da in self.dauid_map.values()]
        self.rtree = shapely.strtree.STRtree(da_polygons)

    def _slice_regions(
            self,
            divisions_offset,
            areas_offset
    ) -> typing.Tuple[typing.Dict[int, typing.Dict], typing.Dict[int, typing.List[typing.Dict]]]:
        # to simplify the slicing logic, take the absolute min/max from the TL/BR pairs...
        zone_tl = (min(self.zone_top_left[0], self.zone_bottom_right[0]),
                   min(self.zone_top_left[1], self.zone_bottom_right[1]))
        zone_br = (max(self.zone_top_left[0], self.zone_bottom_right[0]),
                   max(self.zone_top_left[1], self.zone_bottom_right[1]))
        zone_height = zone_br[0] - zone_tl[0]
        zone_width = zone_br[1] - zone_tl[1]
        div_height = zone_height / self.divisions[0]
        div_width = zone_width / self.divisions[1]
        area_height = div_height / self.areas[0]
        area_width = div_width / self.areas[1]
        cduid_map, dauid_map = {}, {}
        for div_row in range(self.divisions[0]):
            for div_col in range(self.divisions[1]):
                div_top_left = (
                    zone_tl[0] + div_row * div_height,
                    zone_tl[1] + div_col * div_width,
                )
                div_index = divisions_offset
                cduid = covid19geo.utils.get_cduid_from_partial(self.zone_prefix, div_index)
                assert cduid not in cduid_map
                cduid_map[cduid] = []
                internal_areas_offset = areas_offset
                for area_row in range(self.areas[0]):
                    for area_col in range(self.areas[1]):
                        area_top_left = (
                            div_top_left[0] + area_row * area_height,
                            div_top_left[1] + area_col * area_width,
                        )
                        area_bottom_right = (
                            area_top_left[0] + area_height,
                            area_top_left[1] + area_width,
                        )
                        area_index = internal_areas_offset
                        dauid = covid19geo.utils.get_dauid_from_partial(
                            self.zone_prefix, div_index, area_index)
                        da_poly = covid19geo.utils.get_poly_from_bbox(
                            (area_top_left, area_bottom_right))
                        assert dauid not in dauid_map
                        dauid_map[dauid] = {
                            "ptuid": self.zone_prefix,
                            "cduid": cduid,
                            "dauid": dauid,
                            "polygon": da_poly,
                            "wkb": da_poly.wkb,
                        }
                        assert not hasattr(da_poly, "_dauid")
                        setattr(da_poly, "_dauid", dauid)  # for ref keeping in rtree
                        cduid_map[cduid].append(dauid_map[dauid])
                        internal_areas_offset += 1
                divisions_offset += 1
        return dauid_map, cduid_map


class Zone(ZoneBase):
    """Builds a hierarchy of geographic regions for a given zone bounding box.

    This implementation will use real census divisions & dissemination areas based on already-
    parsed census data given as input. The 'zone' here is equivalent to a province/territory
    in terms of region encoding, but it can actually be much smaller than one in practice.

    This class can be used for testing/debug or for zone-level processing inside the simulator.
    """

    def __init__(
            self,
            bounding_box: BoundingBoxType,
            # TODO: swap data extractor for a compatible HDF5 data reader
            census_data: covid19geo.census.CensusDataExtractor,
    ):
        """Parses and prepares divisions/areas based on the provided params.

        Args:
            bounding_box: the bounding box of the zone for which to get divisions/areas.
            census_data: the object that contains the necessary census data.
        """
        super().__init__()
        self.zone_top_left = bounding_box[0]  # lat, lon
        self.zone_bottom_right = bounding_box[1]  # lat, lon
        self.dauid_map, self.cduid_map = self._fetch_regions(bounding_box, census_data)
        da_polygons = [da["polygon"] for da in self.dauid_map.values()]
        self.rtree = shapely.strtree.STRtree(da_polygons)

    @staticmethod
    def _fetch_regions(
            zone_bbox: BoundingBoxType,
            census_data: covid19geo.census.CensusDataExtractor,
    ) -> typing.Tuple[typing.Dict[int, typing.Dict], typing.Dict[int, typing.List[typing.Dict]]]:
        # here, we will actually build a temporary R-Tree and use it to select the DAs/CDs to keep
        # ...yes, this is lazy & inefficient if we create several zones, but it will do for now
        if isinstance(census_data, covid19geo.census.CensusDataExtractor):
            # the census data extractor maps have string-based keys; here, we will convert to ints
            dauid_map = copy.deepcopy(census_data.dauid_map)
            for dauid, da in dauid_map.items():
                da["polygon"] = shapely.wkb.loads(da["wkb"])
                # convert all strings into integers at the same time
                for key in ["ptuid", "cduid", "dauid"]:
                    da[key] = int(da[key])
                # add dauid attrib to poly now (for ref keeping in rtree)
                assert not hasattr(da["polygon"], "_dauid")
                setattr(da["polygon"], "_dauid", int(dauid))
        else:
            raise NotImplementedError
        rtree = shapely.strtree.STRtree([da["polygon"] for da in dauid_map.values()])
        zone_poly = \
            covid19geo.utils.get_poly_from_bbox(zone_bbox)
        query_res = rtree.query(zone_poly)
        zone_dauids = []
        for hit_poly in query_res:
            if hit_poly.intersects(zone_poly):
                zone_dauids.append(hit_poly._dauid)
        dauid_map = {int(dauid): da for dauid, da in dauid_map.items() if int(dauid) in zone_dauids}
        cduid_map = {}
        for dauid, da in dauid_map.items():
            if da["cduid"] not in cduid_map:
                cduid_map[da["cduid"]] = []
            cduid_map[da["cduid"]].append(da)
        return dauid_map, cduid_map


if __name__ == "__main__":
    # TODO: convert to proper tests @@@@
    montreal_dummy_zone = DummyZone(
        bounding_box=covid19geo.utils.get_montreal_island_bbox(),
        division_rows=3, division_cols=3,
        area_rows=10, area_cols=10,
    )
    # 3x3 divisions with 10x10 areas in each = 900 DAUIDs
    assert len(montreal_dummy_zone.get_divisions_map()) == 9
    assert len(montreal_dummy_zone.get_areas_map()) == 900
    mila_dauid = montreal_dummy_zone.get_coord_dauid(covid19geo.utils.get_mila_coords())
    assert mila_dauid != -1

    # by default, we will go up two directories to find the C++ project root;
    # this might not be the correct default place to download stuff if you
    # installed the Python package directly. In that case, make sure you provide
    # a default argument!
    # TODO: @@@@ REPLACE BY HDF5 READER, CSV PARSING IS TOO SLOW
    download_dir = pathlib.Path(__file__).parent.absolute().parents[2] / "data"
    os.makedirs(download_dir, exist_ok=True)
    records_dir_path, boundaries_dir_path = \
        covid19geo.census.download_and_extract_census_data(download_dir)
    # extractor = covid19geo.census.CensusDataExtractor(
    #     census_records_dir_path=records_dir_path,
    #     census_boundaries_dir_path=boundaries_dir_path,
    # )
    import pickle
    # with open("/tmp/extractor.pkl", "wb") as fd:
    #     pickle.dump(extractor, fd)
    with open("/tmp/extractor.pkl", "rb") as fd:
        extractor = pickle.load(fd)
    montreal_real_zone = Zone(
        bounding_box=covid19geo.utils.get_montreal_island_bbox(),
        census_data=extractor,
    )
    assert len(montreal_real_zone.get_divisions_map()) > 0
    assert len(montreal_real_zone.get_areas_map()) > 0
    mila_dauid = montreal_real_zone.get_coord_dauid(covid19geo.utils.get_mila_coords())
    # check if we can match the data obtained manually from censusmapper.ca
    assert mila_dauid == 24661626
