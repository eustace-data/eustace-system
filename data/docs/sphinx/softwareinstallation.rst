Installation
============

The EUSTACE system is written in Python_ 2.7 for use in a distributed memory parallel processing environment running a Linux operating system. 
Development and testing has been carried out on the JASMIN_ high-performance computing facility.

Prerequisites
-------------

The `JASMIN Analysis Platform`_ is required. This is installed and maintained on JASMIN by STFC_.

Dedicated queues and machines within JASMIN were used to efficiently run the system.

In addition the following Python packages are needed:

    pylint_

    sphinx_

    coverage_

These can be installed by creating a virtual environment as described on the `JASMIN help page`_:

.. code-block:: bash    

    virtualenv --system-site-packages [my python path]
    source [my python path]/bin/activate
    pip install pylint
    pip install sphinx
    pip install coverage

Here :code:`[my python path]` should be replaced with a location to store the packages such as in group workspace :code:`/gws/nopw/j04/eustace/users/` `myusername` :code:`/pythonenv`.

This environment must always be active when using the EUSTACE System.  One way to achieve this is to edit :code:`$HOME/.bashrc` and add the following lines:

.. code-block:: bash

    source [my python path]/bin/activate

The Met Office rose tool must be available on the Python and shell paths.  To achieve this on JASMIN append the following lines to :code:`$HOME/.bashrc` :

.. code-block:: bash

    export PATH=$PATH:/apps/contrib/metomi/bin
    export PYTHONPATH=$PYTHONPATH:/apps/contrib/metomi/rose/lib/python

Obtaining the Code
------------------

The EUSTACE system is available in an SVN_ repository hosted by STFC_. To obtain a local copy in a folder called :code:`$HOME/code/svn-eustace/system`:

.. code-block:: bash

    mkdir -p $HOME/code/svn-eustace/system
    svn checkout http://proj.badc.rl.ac.uk/svn/eustace/platform/trunk/system $HOME/code/svn-eustace/system


Configuration
-------------

The system must be available on Python and shell paths, which can be done by adding the following lines to :code:`$HOME/.bashrc`:

.. code-block:: bash

    export PATH=$PATH:$HOME/code/svn-eustace/system/bin
    export PYTHONPATH=$PYTHONPATH:$HOME/code/svn-eustace/system/src

Input and output directories are configured in :code:`$HOME/code/svn-eustace/system/src/eustaceconfig.py`.  The following are set as defaults and can be modified if installation on a system other than JASMIN_ is required:

+------------------------+----------------------------------+---------------------------------------------------------------+
| **Variable name**      | **Meaning**                      | **Default value**                                             |
+------------------------+----------------------------------+---------------------------------------------------------------+
| :code:`WORKSPACE_PATH` | EUSTACE workspace for retrieving | :code:`/gws/nopw/j04/eustace`                                 |
|                        | and storing data                 |                                                               |
+------------------------+----------------------------------+---------------------------------------------------------------+
| :code:`DOCS_PATH`      | Output of developer guide HTML   | :code:`/gws/nopw/j04/eustace/public/developerguide`           |
+------------------------+----------------------------------+---------------------------------------------------------------+
| :code:`CODE_PATH`      | EUSTACE src folder               | Set automatically by examining the file system so should not  |
|                        |                                  | usually require modification.                                 |
|                        |                                  | Example value: :code:`$HOME/code/svn-eustace/system/src`      |
+------------------------+----------------------------------+---------------------------------------------------------------+
| :code:`SYSTEM_PATH`    | EUSTACE system folder            | Set automatically by examining the file system so should not  |
|                        |                                  | usually require modification.                                 |
|                        |                                  | Example value: :code:`$HOME/code/svn-eustace/system/`         |
+------------------------+----------------------------------+---------------------------------------------------------------+

The current values of these configuration variables can be checked from a command prompt by querying using Python, for example:

.. code-block:: bash

    python -c "import eustaceconfig;print eustaceconfig.SYSTEM_PATH"

The system contains native code that must be built before the system can be run.  To build:

.. code-block:: bash

    cd $HOME/code/svn-eustace/system/src
    python setup.py build_ext --inplace
    cd build
    cmake ../cpp
    make
   
System Test
-----------

To run a full system test:

.. code-block:: bash

    eustace_runtests.sh

This may take several minutes.  If successful the final lines should look like:

.. code-block:: bash

    ----------------------------------------------------------------------
    Ran [number of tests] tests in [time]s

    OK

At the time of archiving, there remain a number of known failing tests. These do not impact the working system.

.. _Python: http://www.python.org/
.. _JASMIN: http://www.ceda.ac.uk/projects/jasmin/
.. _`JASMIN Analysis Platform`: http://www.jasmin.ac.uk/services/jasmin-analysis-platform/
.. _STFC: http://www.stfc.ac.uk/
.. _JASMIN help page: http://www.jasmin.ac.uk/faq/#python
.. _pylint: https://pylint.org/
.. _sphinx: http://www.sphinx-doc.org/
.. _coverage: https://coverage.readthedocs.io/en/v4.5.x/
.. _SVN: https://subversion.apache.org/
