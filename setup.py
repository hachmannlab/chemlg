import setuptools
from os import path
import chemlg

here = path.abspath(path.dirname(__file__))

# Get the long description from the README file
with open(path.join(here, 'README.md')) as f:
    long_description = f.read()

if __name__ == "__main__":
    setuptools.setup(
        name='chemlg',
        version=chemlg.__version__,
        author='Mohammad Atif Faiz Afzal, Johannes Hachmann',
        author_email='m27@buffalo.edu, hachmann@buffalo.edu',
        # url='https://github.com/hachmannlab/chemml',
        project_urls={
            'Source': 'https://github.com/hachmannlab/chemlg',
            'url': 'https://hachmannlab.github.io/chemlg/'
        },
        description=
        'ChemLG is a smart and massive parallel molecular library generator for  chemical and materials sciences.',
        long_description=long_description,
        scripts=['lib/chemlgshell'],
        keywords=[
            'Library Generator', 'Molecular Library',
            'Materials Science', 'Drug Discovery'
        ],
        license='BSD-3C',
        packages=setuptools.find_packages(),

        install_requires=[
            'future', 'six',
            'numpy', 'pandas',
            'scipy'
        ],
        extras_require={
            'docs': [
                'sphinx',
                'sphinxcontrib-napoleon',
                'sphinx_rtd_theme',
                'numpydoc',
            ],
            'tests': [
                'pytest',
                'pytest-cov',
                'pytest-pep8',
                'tox',
            ],
        },
        tests_require=[
            'pytest',
            'pytest-cov',
            'pytest-pep8',
            'tox',
        ],
        classifiers=[
            'Development Status :: 4 - Beta',
            'Natural Language :: English',
            'Intended Audience :: Science/Research',
            'Programming Language :: Python :: 2.7',
            'Programming Language :: Python :: 3.5',
            'Programming Language :: Python :: 3.6',
        ],
        zip_safe=False,
    )
