# A similarity metric for macromolecular crystallization conditions

![alt text][logo]

## What is cockatoo?

cockatoo is an implementation of a similarity metric used in the comparison of
macromolecular crystallization conditions (or cocktails).

Documention is online at: http://ubccr.github.io/cockatoo/

## Install with Conda

First download and installed either [Anaconda](https://www.anaconda.com/distribution/) or
[Miniconda](https://docs.conda.io/en/latest/miniconda.html).

Create a new conda environment, checkout the cockatoo source code, and install
the dependencies:

```
    $ conda create -c keiserlab -c rdkit -c sdaxen --name cockatoo-dev e3fp
    $ conda activate cockatoo-dev

    # Note on Linux may want to install ncurses if using Miniconda
    $ conda install -c default ncurses

    $ git clone git://github.com/ubccr/cockatoo.git cockatoo
    $ cd cockatoo
    $ pip install -e .
```

## License

cockatoo is released under the GNU General Public License ("GPL") Version 3.0.
See the LICENSE file.

[logo]: docs/images/cockatoo-logo-lg.jpg "Cockatoo Logo"
