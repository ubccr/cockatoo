try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

config = {
    'name': 'cockatoo',
    'description': 'cockatoo - A similarity metric for macromolecular crystallization conditions',
    'author': 'Andrew E. Bruno',
    'url': '',
    'download_url': '',
    'author_email': 'aebruno2@buffalo.edu',
    'version': '0.0.1',
    'install_requires': ['nose','rdkit'],
    'packages': ['cockatoo'],
    'scripts': [],
    'classifiers': [
        'Development Status :: 3 - Alpha',
        'Environment :: Console',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'Programming Language :: Python',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
     ],
}

setup(**config)
