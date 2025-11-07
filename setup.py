from setuptools import setup, find_packages
import codecs
import os

# Leia o conteúdo do arquivo README.md
with open('README.md', 'r', encoding='utf-8') as arq:
    readme = arq.read()

install_requires = ["scanpy>=1.10.1", "pybiomart", "scikit-image"]

VERSION = '0.0.91'
DESCRIPTION = 'This comprehensive toolkit enables the analysis of multiple spatial transcriptomics datasets, offering a wide range of analytical capabilities. It supports various types of analyses, including detailed plotting and advanced image analysis, to help you gain deeper insights into your spatial transcriptomics data.'

# Configuração do setup
setup(
    name="spatools",
    version=VERSION,
    author="Pedro Videira Pinho, Mariana Boroni",
    author_email="pedrovp161@gmail.com",
    description=DESCRIPTION,
    long_description_content_type="text/markdown",
    long_description=readme,
    packages=find_packages(),
    install_requires=install_requires,  # Melhorar os requests conforme necessário
    python_requires=">=3.9",
    extras_require={  # Dependências opcionais para desenvolvimento
        "dev": ["twine>=5.1.1"]
    },
    keywords=['Python', 'Spatial transcriptomics', 'Spatial', 'transcriptomics',
              'Multi-sample', 'Colocatlization-analysis'
              'Bioinformatics'],
    classifiers=[
        "Development Status :: 4 - Beta",
        "Programming Language :: Python :: 3",
        "Operating System :: Unix",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: Microsoft :: Windows",
        "License :: OSI Approved :: MIT License"
    ]
)
