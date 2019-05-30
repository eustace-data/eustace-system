from eustace.outputformats.ensuredirectory import ensuredirectory

"""Generates json files describing observation inputs at each time step"""

def input_summary(time, output, inputsources):
    """Extract and store information required to produce a set of inputloaders at a time step
    
    Stored in a json file as [isotime, [inputs_0,...,input_n]]
    
    """
    
    from eumopps.catalogue.fileio.formatjson import CatalogueWriterJSON

    ensuredirectory(output)

    outputfile = open(output, 'w')
    
    CatalogueWriterJSON().save(output, [time.isoformat(), inputsources])
    
def merge_input_summaries(summary_files, output):
    """Concatenates lists of inputs into a single json file"""
    
    from eumopps.catalogue.fileio.formatjson import CatalogueReaderJSON, CatalogueWriterJSON
    
    list_of_inputs = []
    
    for filename in summary_files:
        input_date, input_definition = CatalogueReaderJSON().load(filename)
        list_of_inputs.append( [input_date, input_definition] )
    
    ensuredirectory(output)
    
    CatalogueWriterJSON().save(output, list_of_inputs)
    
def generate_input_descriptor(json_file):
    """Return a time ordered dict of input files loaded from a json file"""
    
    import collections
    from eumopps.catalogue.fileio.formatjson import CatalogueReaderJSON 
    
    inputs = CatalogueReaderJSON().load(json_file)
    input_descriptor = collections.OrderedDict(inputs)
    
    return input_descriptor
    
