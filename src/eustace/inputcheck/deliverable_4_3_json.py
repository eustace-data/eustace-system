"""Read EUSTACE Deliverable 4.3 from EUSTACE Wiki and output dataset descriptors.

Each row of output is of the form <dataset name>,<dataset path on CEMS>.
"""

__version__ = "$Revision: 136 $"
__author__ = "Joel R. Mitchelson"

import sys
import json
import deliverable_4_3

dataset_list = deliverable_4_3.get_dataset_list()
json.dump(obj=dataset_list, fp=sys.stdout, default=lambda o: o.__dict__)
