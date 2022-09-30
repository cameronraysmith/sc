from setuptools import setup, find_packages

setup(
    name='sc',
    version='22.09.1',
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        'Click',
    ],
    entry_points={
        'console_scripts': [
            'scpreprocess = sc.scpreprocess:cli',
        ],
    },
)