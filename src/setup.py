"""
my docstring

"""
from distutils.core import setup

setup(
    name='pyrem',
    version='trunk',
    author='Quentin Geissmann',
    author_email= 'quentin.geissmann13@imperial.ac.uk',
    packages=['pyrem',
              ],
    url="https://github.com/gilestrolab/pyrem",
    license="GPL3",
    description='todo', #TODO
    long_description=open('README').read(),
    package_data={'pyrem': ['data/classifiers/*.pkl']},
#    extras_require={
#        'pipes': ['luigi>=1.0.13'],
#    },
    install_requires=[
        "numpy>=1.6.1",
        "matplotlib>=1.1.0",
        "scipy>=0.10.1",
        "pandas>=0.14.0",
        "pywavelets>=0.2.2",
        "joblib>=0.8.2",
        # "scikit-learn>=0.14.1",
        "scikits.samplerate>=0.3.3",
    ],
)
