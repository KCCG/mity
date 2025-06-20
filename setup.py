import setuptools

with open("mitylib/_version.py", "r") as fh:
    version = fh.read().replace("__version__ = ", "").strip('""\n')

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="{NAME}",
    version=version,
    description="A sensitive Mitochondrial variant detection pipeline from WGS data",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/KCCG/mity",
    author="Mark Cowley",
    author_email="mcowley@ccia.org.au",
    license="MIT",
    classifiers=[
        "Development Status :: 4 - Beta",
        "Environment :: Console",
        "Intended Audience :: Science/Research",
        "Intended Audience :: Healthcare Industry",
        "License :: Free for non-commercial use",
        "Natural Language :: English",
        "Operating System :: MacOS",
        "Operating System :: POSIX",
        "Operating System :: Unix",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    keywords="mitochondrial DNA genomics variant SNV INDEL",
    project_urls={
        "Documentation": "https://github.com/KCCG/mity/",
        "Source": "https://github.com/KCCG/mity/",
        "Funding": "http://garvan.org.au/kccg",
    },
    packages=setuptools.find_packages(),
    install_requires=["pysam", "pandas", "xlsxwriter", "scipy", "pyyaml", "vcf2pandas"],
    python_requires=">=3.10",
    package_data={"mitylib": ["annot_mt/*", "annot_chrm/*", "reference/*", "config/*"]},
    include_package_data=True,
    scripts=["mitylib/mity"],
)
