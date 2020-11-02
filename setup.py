import os
from setuptools import find_packages, setup

HERE = os.path.dirname(__file__)
with open(os.path.join(HERE, "README.md")) as o:
    readme = o.read()

setup(
    name='oligoss',
    version='0.0.1',
    description='de novo sequencer for heterogeneous oligomer mixtures',
    long_description=readme,
    long_description_content_type="text/markdown",
    author='ALife Team, Cronin Group',
    author_email='daviddoran20161@gmail.com',
    packages=find_packages(),
    install_requires=["bson", "psutil", "mzmlripper"],
    entry_points={
        "console_scripts": [
            "oligoss=oligoss.__main__:main"
        ]
    },
    license="AGPLv3",
    classifiers=[
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
        "License :: OSI Approved :: GNU Affero General Public License v3 or later (AGPLv3+)",
        "Operating System :: Unix",
        "Development Status :: 3 - Alpha"
    ],
    include_package_data=True
)
