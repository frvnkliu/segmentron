from setuptools import setup, find_packages

project_name = "segmentron_regevlab"
project_version = "9.0.1"

authors = [
    {"name": "Frank Liu", "email": "frankliu@nyu.edu"},
    {"name": "Nicholas Wang", "email": "nyw2006@nyu.edu"},
    {"name": "Antonio Vela Gartner", "email": 'Antonio.VelaGartner@nyulangone.org'},
    {"name": "Oded Regev", "email": "regev@cims.nyu.edu"}
]

description = "Python Package for optimal segmentation of a DNA sequence"
readme = "README.md"
requires_python = ">=3.8"

classifiers = [
    "Programming Language :: Python :: 3",
    "License :: OSI Approved :: MIT License",
    "Operating System :: OS Independent",
]

urls = {
    "Homepage": "https://github.com/frvnkliu/segmentron",
    "Issues": "https://github.com/frvnkliu/segmentron/issues"
}

setup(
    name=project_name,
    version=project_version,
    author=", ".join([author["name"] for author in authors]),
    author_email=authors[0]["email"],  # Using the first author's email for simplicity
    description=description,
    #long_description=open(readme).read(),
    long_description_content_type="text/markdown",
    packages=find_packages(),
    python_requires=requires_python,
    classifiers=classifiers,
    project_urls=urls,
)