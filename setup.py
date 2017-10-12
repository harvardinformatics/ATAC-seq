from setuptools import setup

setup(
    name = 'ATAC-seq',
    version = '0.1',
    description = 'ATAC-seq guidelines',
    author = 'John M. Gaspar',
    author_email = 'jgaspar@fas.harvard.edu',
    license = 'GPLv2',
    packages = ['atacseq'],
    zip_safe = False,
    entry_points = {
        'console_scripts' : [
            'SAMtoBED=atacseq.SAMtoBED:main',
            'removeChrom=atacseq.removeChrom:main',
        ]
    }
)
