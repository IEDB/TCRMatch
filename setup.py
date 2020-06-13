from distutils.core import setup, Extension

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name="tcrmatch",
    version="0.0.1",
    author="Austin Crinklaw",
    author_email="acrinklaw@lji.org",
    description="TCRMatch",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/IEDB/tcrmatch",
    packages=['TCRMatch'],
    package_data={'TCRMatch': ['data/*']},
    include_package_data=True,
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU Affero General Public License v3 or later (AGPLv3+)",
        "Operating System :: Unix",
    ],
    ext_modules=[Extension('tcrmatch_c', ['TCRMatch/tcrmatch.c'], extra_compile_args=['-std=c99'])],
)
