from setuptools import setup, find_packages


with open('requirements.txt', 'r') as f:
    requirements = [line.strip() for line in f.readlines() if not line.startswith("#")]

setup(
    name                 = "covid19geo",
    version              = "0.0.0.dev0",
    url                  = "https://github.com/covi-canada/geospatial",
    description          = "Geospatial helper library for COVI-Canada projects.",
    long_description     = "Geospatial helper library for COVI-Canada projects.",
    classifiers          = [
        "Development Status :: 1 - Planning",
        "Environment :: Console",
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3 :: Only",
        "Topic :: Scientific/Engineering",
        "Topic :: Scientific/Engineering :: Artificial Intelligence",
        "Topic :: Software Development",
        "Topic :: Software Development :: Libraries",
        "Topic :: Software Development :: Libraries :: Python Modules",
    ],
    zip_safe             = False,
    python_requires      = '>=3.7.4',
    install_requires     = requirements,
    packages             = find_packages("src"),
    package_dir          = {'': 'src'},
)
