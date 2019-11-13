from setuptools import setup


setup(
      name='ds9_colormap',    # This is the name of your PyPI-package.
      version='1.0',
      description='Python package to bypass ds9',
      author='Vittorio Ghirardini',
      author_email='vittorio@mpe.mpg.de',
      url="https://github.com/vittoriodx7/ds9_colormap",
      packages=['ds9_colomap'],
      install_requires=[
            'astropy','matplotlib','numpy'
      ],
)
