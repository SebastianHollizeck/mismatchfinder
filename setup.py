from setuptools import setup, find_packages

vsn = {}
with open("./mismatchfinder/__meta__.py") as fp:
    exec(fp.read(), vsn)
version = vsn["__version__"]
description = vsn["description"]

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="mismatchfinder",
    version=version,
    author="Sebastian Hollizeck",
    author_email="sebastian.hollizeck@petermac.org",
    license="GNU",
    description=description,
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/sebastianhollizeck/mismatchfinder",
    packages=find_packages(),
    install_requires=[
        "numpy",
        "ncls",
        "pysam",
        "matplotlib",
        "zarr",
        "pandas",
        "quadprog",
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Healthcare Industry",
        "Intended Audience :: Science/Research",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.6",
    entry_points={"console_scripts": ["mismatchfinder=mismatchfinder.__main__:main"]},
    package_data={"mismatchfinder": ["ext/*.csv"]},
)
