import setuptools

with open("README.md", 'r') as f:
    long_description = f.read()

setuptools.setup(name='distributed-control-fenics',
                 version='1.0.0',
                 description='A module of helper functions that assembles' +
                             '(distributed) control and observation' +
                             'operators using FEniCS.',
                 license="MIT",
                 long_description=long_description,
                 long_description_content_type='text/markdown',
                 url="https://github.com/highlando/distr_control_fenics",
                 author='Jan Heiland',
                 author_email='jnhlnd@gmail.com',
                 packages=['distributed_control_fenics'],  # same as name
                 install_requires=['numpy', 'scipy'],  # extrnl pckgs
                 classifiers=[
                     "Programming Language :: Python :: 3",
                     "Programming Language :: Python :: 2",
                     "License :: OSI Approved :: MIT License",
                     ]
                 )
