
from distutils.core import setup

setup(name='gmatch',
      version='1.0',
      description='Match pairs of objects using Groth\'s algorithm',
      author='Sergio Pascual',
      author_email='sergiopr@fis.ucm.es',
      url='https://bitbucket.org/sergiopr/gmatch',
      license='GPLv3',
      packages=['gmatch'],
      requires=['scipy'],
      classifiers=[
       "Programming Language :: Python",
       'Development Status :: 3 - Alpha',
       "Environment :: Other Environment",
       "Intended Audience :: Science/Research",
       "License :: OSI Approved :: GNU General Public License (GPL)",
       "Operating System :: OS Independent",
       "Topic :: Scientific/Engineering :: Astronomy",
      ],
)
