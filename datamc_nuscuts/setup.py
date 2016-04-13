from setuptools import setup

setup(name='datamc_nuscuts',
      version='0.1',
      description='analyze NusCuts',
      url='http://github.com/aa158/DataMC_NusCuts',
      author='Aaron Markowitz',
      author_email='markowitz@college.harvard.edu',
      license='GNU GPL',
      packages=['datamc_nuscuts'],
      install_requires=[
          'markdown',
      ],
      zip_safe=False,
      test_suite='nose.collector',
      tests_require=['nose'])
