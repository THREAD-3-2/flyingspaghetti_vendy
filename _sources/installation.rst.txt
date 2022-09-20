.. _installation:

==============
 Installation
==============

Virtual environment
===================

If you have a supported Python installation on your computer, you can
install the package in a virtual environment like so:

.. code-block:: bash

    # create a virtual environment (called venv)
    python3 -m venv venv

    # activate virtual environment
    . ./venv/bin/activate
    
    # install packages listed in `requirements.txt`, e.g. sphinx
    pip install -r requirements.txt
    
    # to build the documentation locally, run:
    cd docs
    make doctest # optionally check if your examples work
    make html
