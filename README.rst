vaccontrib
==========

Code to quantify the contributions unvaccinated and vaccinated
subpopulations make towards the effectie reproduction number and new
infections.

-  repository: https://github.com/benmaier/vaccontrib/
-  documentation: http://vaccontrib.readthedocs.io/

.. code:: python

   >>> from vaccontrib.covid import get_reduced_vaccinated_susceptible_contribution_matrix_covid
   >>> R0 = 6
   >>> C = get_reduced_vaccinated_susceptible_contribution_matrix_covid(R0,variant='delta')
   >>> C/C.sum()
   array([[0.38159051, 0.17360365],
          [0.28493897, 0.15986686]])

Also, check out the `tutorial
notebook <https://github.com/benmaier/vaccontrib/blob/main/cookbook/notebooks/covid_examples.ipynb>`__
and an `explanatory notebook including stochastic
simulations <https://github.com/benmaier/vaccontrib/blob/main/cookbook/notebooks/first_examples.ipynb>`__.

Install
-------

.. code:: bash

   git clone git@github.com:benmaier/vaccontrib.git
   pip install ./vaccontrib

``vaccontrib`` was developed and tested for

-  Python 3.6
-  Python 3.7
-  Python 3.8

So far, the package's functionality was tested on Mac OS X and CentOS
only.

Dependencies
------------

``vaccontrib`` directly depends on the following packages which will be
installed by ``pip`` during the installation process

-  ``numpy>=1.17``
-  ``matplotlib>=3.3``
-  ``bfmplot>=0.0.11``

Manuscript
----------

Results found using this software package were reported in a
`preprint <https://medrxiv.org/cgi/content/short/2021.11.24.21266831v1>`__.
To replicate the results, use the scripts in the
`paper_analysis <https://github.com/benmaier/vaccontrib/tree/main/paper_analysis>`__
directory. E.g. to get a comprehensive overview of results for a
parameterset, run

.. code:: bash

   python compute.py DIRNAME1 DIRNAME2

e.g.

.. code:: bash

   python compute.py 00_lower 01_upper

See the help text:

.. code:: bash

   usage: compute.py [-h] [-u RU] [-v RV] [-f] [-R R0] dirs [dirs ...]

   Compute contribution matrices.

   positional arguments:
     dirs                directories for which contributions matrices should be computed

   optional arguments:
     -h, --help          show this help message and exit
     -u RU, --Ru RU      Base R-value of unvaccinated
     -v RV, --Rv RV      Base R-value of vaccinated
     -f, --save-figures  create, show, and save illustrations
     -R R0, --R0 R0      Base R0 to which the contribution will be scaled

Make sure you have
`numpyarray_to_latex <https://github.com/benmaier/numpyarray_to_latex>`__
installed.

Changelog
---------

Changes are logged in a `separate
file <https://github.com/benmaier/vaccontrib/blob/main/CHANGELOG.md>`__.

License
-------

This project is licensed under the `MIT
License <https://github.com/benmaier/vaccontrib/blob/main/LICENSE>`__.
Note that this excludes any images/pictures/figures shown here or in the
documentation.

Contributing
------------

If you want to contribute to this project, please make sure to read the
`code of
conduct <https://github.com/benmaier/vaccontrib/blob/main/CODE_OF_CONDUCT.md>`__
and the `contributing
guidelines <https://github.com/benmaier/vaccontrib/blob/main/CONTRIBUTING.md>`__.
In case you're wondering about what to contribute, we're always
collecting ideas of what we want to implement next in the `outlook
notes <https://github.com/benmaier/vaccontrib/blob/main/OUTLOOK.md>`__.

|Contributor Covenant|

Dev notes
---------

Fork this repository, clone it, and install it in dev mode.

.. code:: bash

   git clone git@github.com:YOURUSERNAME/vaccontrib.git
   make

If you want to upload to PyPI, first convert the new ``README.md`` to
``README.rst``

.. code:: bash

   make readme

It will give you warnings about bad ``.rst``-syntax. Fix those errors in
``README.rst``. Then wrap the whole thing

.. code:: bash

   make pypi

It will probably give you more warnings about ``.rst``-syntax. Fix those
until the warnings disappear. Then do

.. code:: bash

   make upload

.. |Contributor Covenant| image:: https://img.shields.io/badge/Contributor%20Covenant-v1.4%20adopted-ff69b4.svg
   :target: code-of-conduct.md
