
with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="quantumnoise-pkg-christophe-valahu", # Replace with your own username
    version="1.0.0",
    author="Christophe Valahu",
    author_email="christophe.valahu@yahoo.com",
    description="A package to model decoherence due to noise",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/pypa/sampleproject",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)