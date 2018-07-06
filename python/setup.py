from setuptools import setup

setup (
    name = 'pydiskvert',
    version = '180314',
    author = 'Dominik Gronkiewicz',
    author_email = 'gronki@gmail.com',
    description = u"Calculate vertical structure of accretion disks",
    license = "MIT",
    packages = [ 'diskvert' ],
    scripts = [
        'scripts/col2python',
        'scripts/diskvert-cooling2D',
    ],
    entry_points = {
        'console_scripts': [
            'dv-plot-rx=diskvert.plotrx:main_plotmulti',
            'dv-plot-rx-cumul=diskvert.plotrx:main_plotcumul',
        ],
    },
    install_requires = [
        'numpy', 'matplotlib', 'sympy', 'ipython<6.0',
    ],
)
