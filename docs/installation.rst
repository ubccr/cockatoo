Installing
======================
.. highlight:: bash

**cockatoo** is in the `Python Package Index <http://pypi.python.org/pypi/cockatoo/>`_.

Requirements
-------------------

cockatoo requires `pinky <https://github.com/ubccr/pinky>`_.

  $ pip install pinky

The hierarchical clustering module requires python packages: numpy,scipy,matplotlib,brewer2mpl,ete2::
    
  $ pip install numpy
  $ pip install scipy
  $ pip install brewer2mpl
  $ easy_install -U ete2

Installing with pip
-------------------

To install cockatoo via `pip <http://pypi.python.org/pypi/pip>`_::

  $ pip install cockatoo

Installing from source
-----------------------

To install cockatoo from source checkout the latest code from `github <https://github.com/ubccr/cockatoo>`_::

  $ git clone git://github.com/ubccr/cockatoo.git cockatoo
  $ cd cockatoo
  $ python setup.py install
