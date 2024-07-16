from setuptools import setup, find_packages

with open(file='README.md', mode='r') as fh:
    long_description = fh.read()

    setup(
    name='scPAGE',
    version='0.1',
    author='Fabio Dorazio',
    description='Single cell Perturbations and Analysis of Gene Expression',
    long_description=long_description,
    long_description_content_type='text/markdown',
    packages=find_packages(exclude=('docs', 'tests')),
    include_package_data=True,
    python_requires='>=3.7',
    entry_points={
        'console_scripts': ['scPAGE=scPAGE.__main__:main'],
    }
)