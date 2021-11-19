# AlphaPlantImpute2

## Requirements

* Python 3
* NumPy
* Numba

## Installation

    pip install alphaplantimpute2-1.5.3-py3-none-any.whl

## User guide

See `alphaplantimpute2.pdf`

## Build instructions

Run the following to build the Python wheel and user guide. You will need an installation of Latex (https://www.latex-project.org/get) to build the user guide.

    git clone git@github.com:AlphaGenes/AlphaPlantImpute2.git
    git submodule init
    git submodule update
    bash build_pipeline.sh 

The wheel can be found in `dist/` and the user guide at `docs/build/latex/alphaplantimpute2.pdf`
