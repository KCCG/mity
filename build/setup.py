import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

with open("mitylib/_version.py", "r") as fh:
    version = fh.read().replace("__version__ = ", "").strip('""\n')

setuptools.setup(
    name="mitywgs",
    version=version,
    description="A sensitive Mitochondrial variant detection pipeline from WGS data",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/KCCG/mity",
    author="Clare Puttick",
    author_email="clare.puttick@gmail.com",
    license="MIT",
    classifiers=[
        'Development Status :: 4 - Beta',
        'Environment :: Console',
        'Environment :: MacOS X',
        'Intended Audience :: Science/Research',
        'Intended Audience :: Healthcare Industry',
        'License :: Free for non-commercial use',
        'Natural Language :: English',
        "Operating System :: MacOS",
        "Operating System :: POSIX",
        "Operating System :: Unix",
        "Programming Language :: Python :: 3.5",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Topic :: Scientific/Engineering :: Bio-Informatics"
    ],
    keywords='mitochondrial DNA genomics variant SNV INDEL',
    project_urls={
        'Documentation': 'https://github.com/KCCG/mity/',
        'Source': 'https://github.com/KCCG/mity/',
        'Funding': 'http://garvan.org.au/kccg',
    },
    packages=setuptools.find_packages(),
    install_requires=[
        'pyvcf',
        'pysam',
        'pandas',
        'xlsxwriter',
        'pyfastx',
        'scipy'
    ],
    python_requires='>=3.5.4',
    package_data={'mitylib': ['verchew.ini', 'annot/*', 'reference/*']},
    include_package_data=True,
    scripts=["mity"]
)
