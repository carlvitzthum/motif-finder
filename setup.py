import io
from setuptools import setup
from os import path

this_directory = path.abspath(path.dirname(__file__))
with io.open(path.join(this_directory, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

# dependencies must be handled carefully for py2 and py3 compatibility
requires = [
    'biopython<=1.74',
    'numpy>=1.10.1,<=1.16.3',
    'importlib_resources'
]

this_version = io.open(path.join(this_directory, "motiffinder/_version.py")).readlines()[-1].split()[-1].strip("\"'")

setup(
    name = "motiffinder",
    version = this_version,
    packages = ['motiffinder'],
    package_dir = {'motiffinder': 'motiffinder'},
    package_data = {'motiffinder': ['static/*']},
    description = "A novel approach to motif identification in proteins",
    long_description=long_description,
    long_description_content_type='text/markdown',
    url = "https://github.com/carlvitzthum/motif-finder",
    download_url = "https://github.com/carlvitzthum/motif-finder/tarball/" + this_version,
    author = "Carl Vitzthum",
    author_email = "carl.vitzthum@gmail.com",
    license = "MIT",
    keywords = ["bioinformatics", "genomics", "proteomics", "motif", "protein"],
    install_requires = requires,
    setup_requires = requires,
    tests_require = requires,
    test_suite = "test",
    entry_points = {
        'console_scripts': [
             'motiffinder = motiffinder.__main__:main',
        ]
    },
    classifiers = [
        "Programming Language :: Python :: 2.7",
        "Programming Language :: Python :: 3",
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Topic :: Software Development :: Libraries :: Python Modules",
        "Topic :: Scientific/Engineering :: Bio-Informatics"
    ],
    zip_safe = True
)
