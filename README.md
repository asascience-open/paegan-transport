Paegan-Transport
================

The parallelized Lagrangian transport model for NetCDF/OPeNDAP data, written on top of [Paegan](https://github.com/asascience-open/paegan)


Setup
------------------
You are using `virtualenv`, right?

1. Install [virtualenv-burrito](https://github.com/brainsik/virtualenv-burrito)
2. Create virtualenv named "paegan-transport-dev": `mkvirtualenv -p your_python_binary paegan-transport-dev`
3. Start using your new virtualenv: `workon paegan-transport-dev`


Installation
-------------
Paegan-Transport requires python 2.7.x and is available on PyPI.

The best way to install Paegan-Transport is through pip:

```bash
pip install paegan-transport
```

Paegan-Transport requires the following python libraries which will be downloaded and installed through `pip`:

* GDAL==1.9.1
* Fiona==0.8
* paegan==0.9.3
* requests==1.2.3

See the [Paegan](https://github.com/asascience-open/paegan) documentation for installing other dependencies related to Paegan.


Roadmap
--------
* More behaviors


Troubleshooting
---------------
If you are having trouble getting any of the paegan functionality to work, try running the tests:

```bash
git clone git@github.com:asascience-open/paegan-transport.git`
cd paegan-transport
pip install pytest
python -m pytest -s
```

If you want to run the model_controller or shoreline tests, you will need to edit the test files with paths appropriate for your system.

Some tests requires large files that are not in source control.  You can get them here:


Contributors
----------------
* Kyle Wilcox <kwilcox@asascience.com>
* Alex Crosby <acrosby@asascience.com>
* Dave Foster <dfoster@asascience.com>
