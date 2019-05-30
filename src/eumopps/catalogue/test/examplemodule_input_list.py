"""Example module for use with test_system, combining list of multiple data sources."""

def run(title, time, output, inputlist, inputfixed):
    """Run method (called by EUMOPPS)."""

    outputfile = open(output, 'w')
    outputfile.write(title + '\n')
    outputfile.write(time.strftime('%Y%m%d') + '\n')
    for inputpathname in (inputlist + [ inputfixed ]):
        if inputpathname is None:
            outputfile.write('--\n')
        else:
            outputfile.write(open(inputpathname, 'r').read())
    outputfile.close()
