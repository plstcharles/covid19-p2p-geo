import hashlib
import os
import typing

import requests
import shapely.geometry
import tqdm

LatLonPairType = typing.Tuple[float, float]
BoundingBoxType = typing.Tuple[LatLonPairType, LatLonPairType]  # TL, BR


def get_unique_file_by_extension(
        directory_path: typing.AnyStr,
        extension: typing.AnyStr,
):
    """Returns a unique file from a directory given a query extension.

    This function will throw if no file is found with the given extension or if
    multiple files are found with this extension.

    Args:
        directory_path: path to the directory in which the files should be searched.
        extension: the suffix to search the files with.

    Returns:
        A string that corresponds to the unique file found in the directory with the
        specified extension.
    """
    file_paths = os.listdir(directory_path)
    assert sum([p.endswith(extension) for p in file_paths]) == 1, \
        f"unexpected content in directory: {directory_path}\n" \
        f"\t(should contain a single '*{extension}' file)"
    target_file_name = next(p for p in file_paths if p.endswith(extension))
    target_file_path = os.path.join(directory_path, target_file_name)
    return target_file_path


def download_file_with_progress_bar(
        file_url: typing.AnyStr,
        file_path: typing.AnyStr,
        block_size: int = 1024,
):
    """Downloads a file from a given URL with a tqdm-based progress bar.

    Will throw if the number of downloaded bytes does not match the expected
    file size.

    Args:
        file_url: the URL to GET using requests.
        file_path: the path where the file should be saved.
        block_size: block size to download/update using.

    Returns:
        The number of downloaded bytes.
    """
    r = requests.get(file_url, stream=True)
    total_size = int(r.headers.get('content-length', 0))
    progr = tqdm.tqdm(total=total_size, unit='iB', unit_scale=True)
    with open(file_path, "wb") as fd:
        for data in r.iter_content(block_size):
            progr.update(len(data))
            fd.write(data)
    progr.close()
    assert total_size == 0 or progr.n != total_size, "unexpected file size"


def compute_md5_hash(
        file_path: typing.AnyStr,
        chunk_size: int = 4096,
):
    """Computes the md5 hash of a file for checksum purposes."""
    hash_md5 = hashlib.md5()
    with open(file_path, "rb") as f:
        for chunk in iter(lambda: f.read(chunk_size), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()


def get_cduid_from_partial(
        ptuid: int,
        cd_code: int,
) -> int:
    """Returns a census division UID (4-digit number) based off its parent UID."""
    assert 0 <= ptuid < 100, "invalid province/territory UID"
    assert 0 <= cd_code < 100, "invalid census division UID"
    cduid = ptuid * 100 + cd_code
    assert 0 <= cduid < 10000, "messed up logic, should be 4-digit integer"
    return cduid


def get_dauid_from_partial(
        ptuid: int,
        cd_code: int,
        da_code: int,
) -> int:
    """Returns a dissemination area UID (8-digit number) based off parent UIDs."""
    cduid = get_cduid_from_partial(ptuid, cd_code)
    assert 0 <= da_code < 10000, "invalid dissemination area code"
    dauid = cduid * 10000 + da_code
    assert 0 <= dauid < 100000000, "messed up logic, should be 8-digit integer"
    return dauid


def get_poly_from_bbox(
        bbox: BoundingBoxType,
) -> shapely.geometry.Polygon:
    """Returns a shapely polygon for a set of bounding box coordinates."""
    return shapely.geometry.Polygon([
        bbox[0],
        (bbox[0][0], bbox[1][1]),
        bbox[1],
        (bbox[1][0], bbox[0][1]),
    ])


def get_montreal_island_bbox() -> BoundingBoxType:
    """Utility function that returns the geo extent (TL-BR bbox) of the island of Montreal.

    The result can be useful for debugging & for generating a dummy grid of zones for the
    simulator.

    The returned bounding box coordinates are in WGS84 (EPSG:4326), and given as a pair
    of lat-lon pairs.
    """
    # this is a really rough approximation of the island's bounding box, don't rely on it!
    top_left = (45.682214, -74.003708)
    bottom_right = (45.404200, -73.491887)
    return top_left, bottom_right


def get_mila_neighborhood_bbox() -> BoundingBoxType:
    """Utility function that returns the geo extent (TL-BR bbox) of the Mila's neighborhood.

    The result can be useful for debugging & for generating a dummy grid of zones for the
    simulator.

    The returned bounding box coordinates are in WGS84 (EPSG:4326), and given as a pair
    of lat-lon pairs.
    """
    top_left = (45.5343278, -73.6181879)
    bottom_right = (45.5257197, -73.6101165)
    return top_left, bottom_right


def get_mila_coords() -> LatLonPairType:
    """Utility function that returns the geospatial coordinates of Mila's 6666 building.

    The returned bounding box coordinates are in WGS84 (EPSG:4326), and given as a pair
    of lat-lon pairs.
    """
    return 45.530637, -73.613989
