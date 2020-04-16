# Geographic Data Preprocessing

To download, process, and export the census data into HDF5 files, you will first
have to setup a Python environment with the dependencies listed in ``conda-env.yml``.

If Anaconda/Miniconda is installed on your machine, simply run:
```
  conda env create --file conda-env.yml
```
This should create a conda environment named ``covid19-p2p-geo``. To activate it, use:
```
  conda activate covid19-p2p-geo
```
You will then be ready to launch the preprocessing script:
```
  python hdf5_exporter.py
```
All the necessary data will be automatically downloaded and packaged in ``../data``.
