from setuptools import setup, find_packages

package_dir = {'rnmcpy': 'src'}

setup(name='RNMNpy',
      version='0.1',
      description='Reaction Network Monte Carlo - Python',
      url='https://github.com/BlauGroup/RNMC',
      author='Daniel Barter',
      author_email='danielbarter@gmail.com',
      license='LBNL',
      packages=find_packages(),
      install_requires=[
        "setuptools",
        "pymatgen>=2022.3.7",
        "monty>=3.0.2",
        "numpy>=1.20.1",
        "matplotlib>=3.5.0",
        "HiPRGen",
      ],
)
