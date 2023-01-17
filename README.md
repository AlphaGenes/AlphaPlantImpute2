# AlphaPlantImpute2

## Requirements

* Python 3
* NumPy
* Numba

## Installation

    wget https://github.com/AlphaGenes/AlphaPlantImpute2/blob/master/alphaplantimpute2-1.5.3-py3-none-any.whl
    pip install alphaplantimpute2-1.5.3-py3-none-any.whl

## User guide

See `docs/source/index.rst` or `alphaplantimpute2.pdf` in this repository.

## Build instructions

Run the following to build the Python wheel and user guide. You will need an installation of [Sphinx](https://www.sphinx-doc.org) and [LaTex](https://www.latex-project.org/get) to build the user guide.

    git clone https://github.com/AlphaGenes/AlphaPlantImpute2.git
    cd AlphaPlantImpute2
    git submodule init
    git submodule update
    bash build_pipeline.sh
    pip install dist/AlphaPlantImpute2*.whl

The wheel can be found in `dist/` and the user guide at `docs/build/latex/alphaplantimpute2.pdf`
