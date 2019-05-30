#!/bin/bash

key=ATSR2

echo -e "Counting ${key} data"

for index in `seq 5 9` ; do input=`ls /gws/nopw/j04/eustace/data/incoming/UoR_sst/l3c/${key}/199${index}/*/* | wc -l` ;  output=`ls /gws/nopw/j04/eustace/data/internal/D2.2/R001115/20180405/ocean/199${index}/* | wc -l`; real_output=`echo ${output}/2 | bc -l` ; echo -e "199${index} input=$input,  output=$real_output"; done

for index in `seq 0 3` ; do input=`ls /gws/nopw/j04/eustace/data/incoming/UoR_sst/l3c/${key}/200${index}/*/* | wc -l` ;  output=`ls /gws/nopw/j04/eustace/data/internal/D2.2/R001115/20180405/ocean/200${index}/* | wc -l`; real_output=`echo ${output}/2 | bc -l` ; echo -e "200${index} input=$input,  output=$real_output"; done


key=AATSR

echo -e "Counting ${key} data"

for index in `seq 2 9` ; do input=`ls /gws/nopw/j04/eustace/data/incoming/UoR_sst/l3c/${key}/200${index}/*/* | wc -l` ;  output=`ls /gws/nopw/j04/eustace/data/internal/D2.2/R001115/20180405/ocean/200${index}/* | wc -l`; real_output=`echo ${output}/2 | bc -l` ; echo -e "200${index} input=$input,  output=$real_output"; done

for index in `seq 10 12` ; do input=`ls /gws/nopw/j04/eustace/data/incoming/UoR_sst/l3c/${key}/20${index}/*/* | wc -l` ;  output=`ls /gws/nopw/j04/eustace/data/internal/D2.2/R001115/20180405/ocean/20${index}/* | wc -l`; real_output=`echo ${output}/2 | bc -l` ; echo -e "20${index} input=$input,  output=$real_output"; done