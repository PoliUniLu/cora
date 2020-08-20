from setuptools import setup


def readme():
    with open('README.rst') as f:
        return f.read()

setup(name='cora',
      version='0.1',
      description='Boolean minimization alghorithms',
      url='',
      author='Zuzana Sebechlebska',
      author_email='zuzanasebechlebska@gmail.com',
      license='',
      packages=['cora'],
      install_requires=[
	'numpy',
	'pandas'
    ]
      )
