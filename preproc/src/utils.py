import os
import typing


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
