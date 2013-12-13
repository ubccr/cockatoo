from distutils.core import setup

VERSION = '0.0.1'

setup(
    name='cockatoo',
    description='cockatoo - A similarity metric for macromolecular crystallization conditions',
    long_description='cockatoo - A similarity metric for macromolecular crystallization conditions',
    author='Andrew E. Bruno',
    url='https://github.com/ubccr/cockatoo',
    license='GNU General Public License v3 (GPLv3)',
    download_url='https://github.com/ubccr/cockatoo',
    author_email='aebruno2@buffalo.edu',
    version=VERSION,
    requires=['nose'],
    packages=['cockatoo'],
    package_data={'cockatoo': ['data/*.csv']},
    scripts=['bin/cockatoo-convert', 'bin/cockatoo-hclust'],
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Environment :: Console',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'Programming Language :: Python',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
     ]
)
