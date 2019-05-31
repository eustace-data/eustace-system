Coding standards
================

Languages
----------
The programming language for the end-to-end system is `Python 2.7.3`_.

Modules written in other languages may be called from python modules. 

`BASH`_ script may be used to coordinate processes as required for submission to the host computer system.


Directory structure
-------------------

The directory structure for the system will adhere to an accepted layout for Python programs.
It is structured as follows:

+---------------------------------------+-----------------------------------------------------+
| **Directory**	                        | **Contents**                                        |
+---------------------------------------+-----------------------------------------------------+
| ``bin``                               | Executable scripts (Python and BASH)                |
+---------------------------------------+-----------------------------------------------------+
| ``data``                              | Auxiliary files required to run the system, such as |
|                                       | configuration parameter files and documentation.    |
+---------------------------------------+-----------------------------------------------------+
| ``data/docs``	                        | System documentation                                |
+---------------------------------------+-----------------------------------------------------+
| ``src``                               | Top-level directory for inclusion in Python path    |
+---------------------------------------+-----------------------------------------------------+
| ``src/setup.py``                      | Installation script                                 |
+---------------------------------------+-----------------------------------------------------+
| ``src/eustace``                       | Primary directory for Python source code            |
+---------------------------------------+-----------------------------------------------------+
| ``src/eustace/`` *module*             | Python module directory                             |
+---------------------------------------+-----------------------------------------------------+
| ``src/eustace/`` *module* ``/test``   | Tests corresponding to Python module                |
+---------------------------------------+-----------------------------------------------------+
| ``src/roses``                         | Primary directory for rose suite and cylc code      |
|                                       | relating to Fullstace and Satstace components       |
+---------------------------------------+-----------------------------------------------------+

In addition there is a ``src/eumopps`` directory structure containing EUMOPPS
(the EUSTACE Met Office Processing and Provenance System).

Style guide
-----------

The `PEP8 style guide for python code`_ will be used with some modifications, extensions, and exceptions.

Verification is done using the `pylint`_ package with appropriate customisation.


Modifications and extensions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Whitespace**

The PEP8 guideline entitled `Whitespace in Expressions and Statements`_ applies except where a tabular form of code is helpful for readability, for example when writing a matrix of numerical constants.  Since these are a common occurence in EUSTACE code, this guideline is not verified automatically.

**String representation**

We represent strings using single quotation marks: ``'some text here'``

**Data classes**

Where data must be grouped together, classes are used rather than named tuples. So where we have a method ``some_calculation(complicated_data)`` which expects attributes of ``complicated_data.some_stuff``, ``complicated_data.other_things``, ``complicated_data.one_more_thing`` then we would define a class like:

.. code-block:: python

    class ComplicatedData(object):
        """Represent complicated data required for calculation."""

        def __init__(self, some_stuff, other_things, one_more_thing):
            self.some_stuff = some_stuff
            self.other_things = other_things
	    self.one_more_thing = one_more_thing

**Test names and layout**

The names of a test class should be of the form Test\ *ClassUnderTest*, for example a test of a class called ``FileReader`` should be called ``TestFileReader``.  If the test suite does not apply to a class then *ClassUnderTest* should be replaced with an expression summarising the main purpose.

Test modules should be placed in the ``test`` of the package under test.  For example if EUSTACE output formats are in package ``eustace.outputformats`` then the corresponding tests will be in ``eustace.outputformats.test``.

Documentation strings are not required in test classes (classes derived from unittest.TestCase) if the class name and method names are meaningful enough to allow the purpose of the tests to be inferred.

**Considered use of single character variable names**

Local variable names of a single character are not generally recommended as they can reduce readability, however they may be used in the following circumstances:

  - Where the single character helps scientific understanding because it corresponds to a convention used in corresponding publications, for example using ``t`` for time.

  - In test cases where a single character helps to keep the test code clear and concise.

Exceptions
~~~~~~~~~~

The following PEP8 guidelines do *not* apply:

+-----------------------------+-----------------------------------------------------+
| **Guideline ignored**       | **Reason**                                          |
+-----------------------------+-----------------------------------------------------+
| Maximum line length         | Reduces readability for some expressions,           |
|                             | particularly in test suites.                        |
+-----------------------------+-----------------------------------------------------+
| Relative imports            | Relative imports are helpful particularly           |
|                             | if we may consider using a subcomponent of the      |
|                             | project in a stand-alone setting.                   |
+-----------------------------+-----------------------------------------------------+


Comments in source code
-----------------------

All source code files should have a header indicating the main purpose of the file.

Python packages must have an ``__init__.py`` file with a comment at the top indicating the purpose of the package.

All important classes, methods and global variables should have a comment.

Where possible code should be written to be self-documenting by use of meaningful variable names.  Comments must be used if the variable names alone would be ambiguous.  Additional comments to clarify design choices are encouraged.

Comments on packages, classes, methods, and global variables should be compatible with the `sphinx`_ auto-documentation package.


Testing
-------

There should be full test coverage of methods using the unittest package.

Test coverage is verified using:

``coverage run -m nose``

followed by

``coverage report`` 

The range of parameters tested depends on the application and must be selected appropriately by the person writing the code.  This choice is also subject to code review.


Code review
-----------

Platform architecture code should be reviewed by software engineers other than those who wrote it.

Underlying scientific methods are subject to peer review, and it is for the publishing scientist to provide prototype code whose output corresponds to peer-reviewed works, and to make a written assertion.


.. _Python 2.7.3: https://www.python.org/download/releases/2.7.3/
.. _PEP8 style guide for python code: https://www.python.org/dev/peps/pep-0008/
.. _pylint: https://pylint.org/
.. _sphinx: http://www.sphinx-doc.org/
.. _Whitespace in Expressions and Statements: https://www.python.org/dev/peps/pep-0008/#whitespace-in-expressions-and-statements
.. _BASH: https://www.gnu.org/software/bash/
