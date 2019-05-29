Development process
===================

System requirements document
----------------------------

A document detailing key system requirements is maintained within the version control system and should be cross-referenced by tickets where appropriate.


Product specification document
------------------------------

A document specifying technical product details (e.g. file formats) is maintained within the version control system.  This should be cross-referenced by tickets where appropriate.  In addition a program should be able to read the document and verify that specifications are met.


Developer workflow
------------------

All source code must be kept in a version control system, and a ticketing system is linked to this.  Workflow for addressing a new feature is as follows:

1. Create a ticket indicating the feature or bug to be addressed.  The name of the ticket should follow agile conventions for user story names, and therefore have the following form:

      As a *role* I can *task* so that *result*.
 
  A list of project milestones and tasks should be maintained and cross referenced from tickets.

  Where relevant, the ticket should cross reference documentation or prototype code which define the feature or bug being addressed.

  The ticket should convey a clear idea of how we determine acceptance.

2. Tests and code should be written to fulfill the ticket.

3. The `sphinx`_ documentation system should be run to generate updated documentation, read to verify that comments appear correctly, and comments revised if necessary.

4. Once tests are sufficient and pass, the code should be committed to the code repository, and the logged with a cross-reference to the ticket.

5. Comments should be placed on the ticket cross-referencing the repository version containing the changes.

6. Code should be submitted for code review where applicable.

7. Steps 2-7 should be repeated as required to adjust code, documentation and tests to achieve the desired outcome.

8. The ticket may be closed by the person who wrote the programs provided they are satisfied that compliance is met, and sufficient review undertaken.

.. _sphinx: http://www.sphinx-doc.org/
