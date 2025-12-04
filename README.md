# lib-maptial

## INSTALL
```
python -m pip install maptial
```

## For developers

### Either create a virtual environment
```
python3 -m venv .env-maps
source .env-maps/bin/activate 
pip install --upgrade pip
pip install -r requirements.txt
```

### Or an editable install in a conda environment
Note conda
 could be also mamba or micromamba or other conda-like tool.

To install `lib-maptial` in editable mode using conda:

1. Create and activate a conda environment (if you don't have one):
	```bash
	conda create -n myenv python=3.10
	conda activate myenv
	```

2. Install pip (if not already available):

	conda install pip
	```

3. Install lib-maptial as editable:
	```bash
	python -m pip install -e .
	```
	Run this command from the `lib-maptial` directory (where `pyproject.toml` is located).

Any changes you make to the library will be reflected immediately in your environment.
