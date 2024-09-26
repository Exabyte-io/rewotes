from setuptools import setup, find_packages

VERSION = '0.0.1'
DESCRIPTION = 'A user friendly package that facilitates the convergence of the plane-wave energy cutoff parameter w.r.t total energy.'

if __name__ == '__main__':
    setup(name='convergence_tracker',
          version=VERSION,
          description=DESCRIPTION,
          author='Brendan Smith',
          author_email='brendansmithphd@gmail.com',
          requires=['Python(>3.3)'],
          install_requires=['Python(>3.3)'],
          packages=find_packages()  
        )   
