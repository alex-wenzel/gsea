from os import walk

from setuptools import setup

NAME = 'gsea'

packages = []
for dp, dns, fns in walk(NAME):
    if dp.split(sep='/')[-1] not in ['.git', '__pycache__']:
        packages.append(dp)
setup(
    name=NAME,
    version='0.0.1',
    description='Gene Set Enrichment Analysis',
    long_description='Library for gene-set enrichment analysis',
    url='https://github.com/KwatME/gsea',
    author='(Kwat) Huwate Yeerna',
    author_email='kwatme8@gmail.com',
    license='LICENSE',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.6',
        'Natural Language :: English',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
    keywords='Bioinformatics, Genomics, Transcriptome, RNA, Statistics',
    packages=packages,
    python_requires='>=3.3',
    install_requires=[
        'numpy',
        'pandas',
        'matplotlib',
    ])
