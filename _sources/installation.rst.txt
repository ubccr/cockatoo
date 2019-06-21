Installing
======================
.. highlight:: bash

**cockatoo** is in the `Python Package Index <http://pypi.python.org/pypi/cockatoo/>`_.

Installing with Conda
----------------------

First download and installed either `Anaconda <https://www.anaconda.com/distribution/>`_ or
`Miniconda <https://docs.conda.io/en/latest/miniconda.html>`_.

Create a new conda environment, checkout the cockatoo source code, and install
the dependencies::

    $ conda create -c keiserlab -c rdkit -c sdaxen --name cockatoo-dev e3fp
    $ conda activate cockatoo-dev

    # Note on Linux may want to install ncurses if using Miniconda
    $ conda install -c default ncurses

    $ git clone git://github.com/ubccr/cockatoo.git cockatoo
    $ cd cockatoo
    $ pip install -e .

Installing with pip
----------------------

To install cockatoo via `pip <http://pypi.python.org/pypi/pip>`_::

  $ pip install cockatoo

Requirements
----------------------

cockatoo requires `e3fp <https://github.com/keiserlab/e3fp>`_, `RDKit <http://rdkit.org/>`_ and `marshmallow <http://marshmallow.readthedocs.org>`_::

  # Note e3fp requires RDKit
  $ pip install e3fp
  $ pip install marshmallow

cockatoo's command line tool requires `click <http://click.pocoo.org/>`_::

  $ pip install Click

The hierarchical clustering module requires python packages: numpy,scipy,matplotlib,brewer2mpl::
    
  $ pip install numpy
  $ pip install scipy
  $ pip install brewer2mpl

Installing from source
-----------------------

To install cockatoo from source checkout the latest code from `github <https://github.com/ubccr/cockatoo>`_::

  $ git clone git://github.com/ubccr/cockatoo.git cockatoo
  $ cd cockatoo
  $ python setup.py install
