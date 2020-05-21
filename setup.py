import setuptools
from Cython.Build import cythonize

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="tcrmatch",
    version="0.0.1",
    author="Austin Crinklaw",
    author_email="acrinklaw@lji.org",
    description="TCRMatch",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/IEDB/tcrmatch",
    packages=setuptools.find_packages(),
    package_data={'TCRMatch': ['data/*']},
    include_package_data=True,
    install_requires=[
        'pandas',
        'Cython',
        'numpy',
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU Affero General Public License v3 or later (AGPLv3+)",
        "Operating System :: Unix",
    ],
    ext_modules=cythonize("TCRMatch/mait_match.pyx", annotate=True),
)
