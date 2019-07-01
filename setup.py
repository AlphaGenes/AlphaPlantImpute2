from setuptools import setup, Extension, find_packages

setup(
    name="Program Name",
    version="0.0.1",
    author="Andrew Whalen",
    author_email="awhalen@roslin.ed.ac.uk",
    description="Description",
    long_description="Description",
    long_description_content_type="text/markdown",
    packages=['program'],
    package_dir={'': 'src'},
    classifiers=[
        "Programming Language :: Python :: 3",
    ],
    entry_points = {
    'console_scripts': [
        'Program=program.program:main',
        ],
    },
    install_requires=[
        'numpy',
        'numba',
    ]
)
