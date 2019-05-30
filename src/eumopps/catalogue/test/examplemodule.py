"""Example module for use with test_system."""

def run(title, time, output, inputdaily, inputfixed):
    """Run method (called by EUMOPPS)."""

    outputfile = open(output, 'w')
    outputfile.write(title + '\n')
    outputfile.write(time.strftime('%Y%m%d') + '\n')
    for inputpathname in [ inputdaily, inputfixed ]:
        outputfile.write(open(inputpathname, 'r').read())
    outputfile.close()

def run_region(title, region_index, measurementfilelist_read, output):
    """Run method (called by EUMOPPS)."""
    print region_index, time
    outputfile = open(output, 'w')
    outputfile.write(title + '\n')
    
    outputfile.write(str(region_index) + '\n')
    
    outputfile.close()
