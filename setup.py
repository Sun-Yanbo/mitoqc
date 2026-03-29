from setuptools import setup, find_packages

setup(
    name="mitoqc",
    version="0.2.0",
    description="Systematic quality control for metazoan mitochondrial genomes",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    author="MitoQC Development Team",
    python_requires=">=3.9",
    packages=find_packages(),
    install_requires=[
        "biopython>=1.79",
        "pandas>=1.5",
        "numpy>=1.23",
    ],
    extras_require={
        "dev": ["pytest>=7.0", "pytest-cov"],
    },
    entry_points={
        "console_scripts": [
            "mitoqc=mitoqc.__main__:main",
        ],
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Intended Audience :: Science/Research",
    ],
)
