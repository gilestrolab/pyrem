from distutils.core import setup

setup(
    name='pyrem',
    version='trunk',
    author='Quentin Geissmann',
    author_email= 'quentin.geissmann13@imperial.ac.uk',
    packages=['pyrem',
              'pyrem.signal',
              ],
    url="https://github.com/gilestrolab/pyrem",
    license="GPL3",
    description='todo', #TODO
    long_description=open('README').read(),
#    extras_require={
#        'pipes': ['luigi>=1.0.13'],
#    },
    install_requires=[
        "numpy>=1.6.1",
        "matplotlib>=1.1.0",
        "scipy>=0.10.1",
    ],
)
