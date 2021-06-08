import setuptools

with open("src/evoseg/README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="evoseg-tecdatalab",
    version="0.0.1a",
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
    install_requires=['numpy>=1.19.5','pandas>=1.1.5', 'scikit-image>=0.17.2', 'scikit-learn>=0.24.2','scipy>=1.5.4', ' mrcfile>=1.3.0'],
)
