from setuptools import setup

setup(
    name = 'ATAC-seq',
    version = '0.1',
    description = 'Attack your seqs',
    author = 'John Gaspar',
    author_email = 'john_gaspar@harvard.edu',
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
