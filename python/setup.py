import os
from setuptools import setup

# Utility function to read the README file.
# Used for the long_description.  It's nice, because now 1) we have a top level
# README file and 2) it's easier to type in the README file than to put a raw
# string in below ...
def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(
    name = "ndlutil",
    version = "0.0.4",
    author = "Neil Lawrence, James Hensman",
    author_email = "james.hensman@gmail.com",
    description = ("A kernel toolbox"),
    license = "BSD",
    keywords = "machine-learning utilities",
    url = "http://TODO",
    packages=['ndlutil'],
    long_description=read('README'),
    classifiers=[
        "Development Status :: 1 - Alpha",
        "Topic :: Machine Learning",
        "License :: OSI Approved :: BSD License",
    ],
)
