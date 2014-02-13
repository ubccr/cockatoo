Installing
======================
.. highlight:: bash

**cockatoo** is in the `Python Package Index <http://pypi.python.org/pypi/cockatoo/>`_.

Requirements
-------------------

cockatoo requires `RDKit <http://rdkit.org/>`_. Debian/Ubuntu users can install RDKit via apt-get::

  $ sudo apt-get install python-rdkit librdkit1 rdkit-data

For information on building RDKit from source see the `manual <https://github.com/rdkit/rdkit/raw/master/Docs/Book/RDKit.pdf>`_.

Installing with pip
-------------------

To install cockatoo via `pip <http://pypi.python.org/pypi/pip>`_::

  $ pip install cockatoo

Installing from source
-----------------------

To install cockatoo from source checkout the latest code from github::

  $ git clone git://github.com/ubccr/cockatoo.git cockatoo
  $ cd cockatoo
  $ python setup.py install
