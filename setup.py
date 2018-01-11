from os import walk

from setuptools import setup

NAME = 'gsea'
URL = 'https://github.com/KwatME/gsea'

packages = []
for dp, dns, fns in walk(NAME):
    if dp.split(sep='/')[-1] not in (
            '.git',
            '__pycache__', ):
        packages.append(dp)
setup(
    name=NAME,
    version='0.0.5',
    description='Library for gene-set enrichment analysis.',
    long_description='See {} to learn more.'.format(URL),
    url=URL,
    author='(Kwat) Huwate Yeerna',
    author_email='kwatme8@gmail.com',
    license='LICENSE',
    classifiers=(
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.6',
        'Natural Language :: English',
        'Topic :: Scientific/Engineering :: Bio-Informatics', ),
    keywords='Bioinformatics, Gene Set, Enrichment, GSEA',
    packages=packages,
    python_requires='>=3.3',
    install_requires=(
        'numpy',
        'pandas',
        'matplotlib', ))
