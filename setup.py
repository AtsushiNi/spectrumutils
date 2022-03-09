import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="spectrumutils",
    version="0.0.1",
    author="Niiihama Atsushi",
    author_email="atsushi9731@gmail.com",
    description="A Python package for spectroscopy",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/AtsushiNi/spectrumutils",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    package_dir={"spectrumutils": "spectrumutils"},
    packages=setuptools.find_packages(where="spectrumutils"),
    python_requires=">=3.6",
)
