from setuptools import setup

setup (
    name = 'diskvert',
    version = '251208',
    author = 'Dominik Gronkiewicz',
    author_email = 'gronki@gmail.com',
    description = u"Calculate vertical structure of accretion disks",
    license = "MIT",
    packages = [ 'diskvert' ],
    scripts = [
        'scripts/col2python',
        'scripts/diskvert-cooling2D',
        'scripts/diskvert-plot',
    ],
    entry_points = {
        'console_scripts': [
            'diskvert-random=diskvert.random:main',
        ],
    },
    install_requires = [
        'numpy', 'matplotlib',
    ],
)
