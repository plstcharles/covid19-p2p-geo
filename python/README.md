# COVI-Canada Geographic Data Preprocessing Package

To download, process, and export the census data into HDF5 files, you will first
have to setup a Python environment with the dependencies listed in ``conda-env.yml``.
Using a conda environment is recommended over a regular virtualenv due to the use of
GDAL, which is very often broken in PyPi, and works out-of-the-box in conda-forge.

If Anaconda/Miniconda is installed on your machine, simply run:
```
  conda env create --file conda-env.yml
```
This should create a conda environment named ``covid19-p2p-geo``. To activate it, use:
```
  conda activate covid19-p2p-geo
```
You will then be ready to install the python package and launch the preprocessing script:
```
  pip install -e . --no-deps
  python -m covid19geo.hdf5_exporter
```
All the necessary data will be automatically downloaded and packaged in ``../data``.

