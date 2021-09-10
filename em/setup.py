import setuptools

with open("src/evoseg/README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="evoseg-tecdatalab",
    version="0.0.6",
    author="Manuel Zumbado Corrales",
    author_email="manzumbado@ic-itcr.ac.cr",
    description="EvoSeg is a tool for automatic subunit segmentation of EM Maps",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/tecdatalab/biostructure",
    project_urls={
        "Bug Tracker": "https://github.com/tecdatalab/biostructure/issues",
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: Apache Software License",
        "Operating System :: OS Independent",
    ],
    package_dir={"": "src"},
    packages=setuptools.find_packages(where="src"),
    python_requires=">=3.6",
)
