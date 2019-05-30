"""
Example to read/write config files in JSON format to/from Python objects.

Output should look like:

Loaded the dining configuration: sushi.json
--Food:  norimaki
--Drink:  saki
--Cutlery list:  ['knife', 'fork']
--Comment:  Chopsticks would have been easier.

Loading a modified configuration:
--Food:  norimaki
--Drink:  saki
--Cutlery list:  ['chopstick', 'chopstick']
--Comment:  Ah now that's better.
"""

from eumopps.jsonobjects import jsonobjects

from tempfile import NamedTemporaryFile

class FineDining:
    """A class describing a fine dining experience."""

    def __init__(self, food, drink, cutlery):
        self.food = food
        self.drink = drink
        self.cutlery = cutlery

class Cutlery:
    """Class describing details of cutlery required."""

    def __init__(self, names, comment):
        self.names = names
        self.comment = comment

def main():
    """Entry point for this program."""

    # Load the sushi configuration and show parameters
    filename = 'sushi.json'
    dining = jsonobjects.load(open(filename, 'r'))
    print 'Loaded the dining configuration:', filename
    print '--Food: ', dining.food
    print '--Drink: ', dining.drink
    print '--Cutlery list: ', dining.cutlery.names
    print '--Comment: ', dining.cutlery.comment
    print ''

    # Modify the object
    dining.cutlery = Cutlery([ 'chopstick', 'chopstick' ], 'Ah now that\'s better.')

    # Save new version as a temporary text file
    modifiedfile = NamedTemporaryFile(
        prefix='example_json_config_', 
        suffix='.json')
    jsonobjects.save(open(modifiedfile.name, 'w'), dining)

    # Load it again 
    print 'Loading a modified configuration:'
    newexperience = jsonobjects.load(open(modifiedfile.name, 'r'))
    print '--Food: ', newexperience.food
    print '--Drink: ', newexperience.drink
    print '--Cutlery list: ', newexperience.cutlery.names
    print '--Comment: ', newexperience.cutlery.comment    
    print ''
    

# Call main entry point
if __name__ == '__main__':
    main()
