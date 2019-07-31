# motif-finder
A novel approach to motif identification in proteins

### Important

This code is in BETA and features may change without notice.

### Usage

Provide a text query file, which is required. If not using manual motifs, also provide one more motif files using the `-m` flag to generate motifs dynamically.

```
pip install motiffinder
# main program
motiffinder <query file> -m <optional motif file>

# for help
motiffind -h

# for version
motiffinder -v
```

Or, clone this repository and set up with:

```
pip install -r requirements.text
python -m motiffinder <query file> -m <optional motif file>

# alternately, use setup.py to get the motiffinder package
python setup.py install
```

### Contributors

Written by Carl Vitzthum and Prof. Andrea Tilden
