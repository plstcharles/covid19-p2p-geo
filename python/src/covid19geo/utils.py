import hashlib
import os
import typing

import requests
import tqdm


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
