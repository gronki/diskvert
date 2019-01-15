from setuptools import setup

setup (
    name = 'diskvert',
    version = '190115',
    author = 'Dominik Gronkiewicz',
    author_email = 'gronki@gmail.com',
    description = u"Calculate vertical structure of accretion disks",
    license = "MIT",
    packages = [ 'diskvert' ],
    scripts = [
        'scripts/col2python',
        'scripts/diskvert-cooling2D',
        'scripts/diskvert-random',
        'scripts/diskvert-plot',
    ],
    install_requires = [
        'numpy', 'matplotlib',
    ],
)
