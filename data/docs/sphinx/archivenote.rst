Archive Note
============

This documentation relates to the development of research code for the EUSTACE project.

The code was developed on the JASMIN CEMS computing resource.

This documentation represents part of the archive of this project. 

The EUSTACE system comprises:

EUMOPPS

SatGrid

SatStace

FullStace

Output formats

Postprocessing

Preprocessing

Surface-air model

Analysis (Advanced standard analysis)


Elements of experimental and in development code and associated comments have 
been retained within this archive where they add context to the final work 
or may provide impetus for future work.

Note that several modules refer to the sksparse cholmod module.
This is not currently installed by the documented installation process
and was not used in this version of the analysis.

This release of the code contains a number of unit tests that currently fail.
These are known to the developers and largely relate to the cholmod module
that was only installed at certain stages of development on selected machines.
It is not installed as part of this published setup, therefore related tests 
will fail. This does not impact the current analysis.
