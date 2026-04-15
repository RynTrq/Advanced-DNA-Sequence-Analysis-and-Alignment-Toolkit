"""Compatibility packaging shim for older build frontends."""

from setuptools import find_packages, setup


setup(
    name="advanced-dna-toolkit",
    version="1.0.0",
    description="DNA sequence alignment, edit distance, and Burrows-Wheeler exact matching toolkit.",
    long_description=open("README.md", encoding="utf-8").read(),
    long_description_content_type="text/markdown",
    packages=find_packages("src"),
    package_dir={"": "src"},
    python_requires=">=3.9",
    entry_points={"console_scripts": ["dna-toolkit=dna_toolkit.cli:main"]},
)

