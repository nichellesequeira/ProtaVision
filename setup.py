from setuptools import setup, find_packages

setup(
    name='protavision',
    version='0.1',
    packages=find_packages(where='src'),
    package_dir={'': 'src'},
    install_requires=[
        'biopython',
        'matplotlib',
        'pandas' ,
        'py3Dmol' ,
        # Ajoutez d'autres dépendances ici
    ],
    entry_points={
        'console_scripts': [
            # Vous pouvez définir des scripts CLI ici si nécessaire
        ],
    },
)
