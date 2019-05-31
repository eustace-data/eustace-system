Documentation
=============

Sphinx
------

Documentation is stored and built using Sphinx

Files are stored under:

``/home/users/$USER/code/svn-eustace/system/data/docs/sphinx``


Process
-------


Build the documentation

``cd /home/users/$USER/code/svn-eustace/system/bin``

``sh eustace_docbuild.sh``

A number of errors will be reported during documentation build.
These are expected, and relate to the cholmod module that is not
available in the current build.


