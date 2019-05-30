"""Test JSON configuration objects."""

import unittest
from eumopps.jsonobjects.jsonobjects import save_object
from eumopps.jsonobjects.jsonobjects import load_object
from eumopps.jsonobjects.jsonobjects import save
from eumopps.jsonobjects.jsonobjects import load
from datetime import datetime
from StringIO import StringIO

class SomeDetail(object):
    """Example class with some details."""

    def __init__(self, onething, anotherthing):
        """Construct from the given parameters."""
        self.onething = onething
        self.anotherthing = anotherthing

class SomeMultiLayerObject(object):
    """Example object which might hold sub-objects."""

    def __init__(self, a, b, c):
        self.a = a
        self.b = b
        self.c = c

class TestJSONObjects(unittest.TestCase):
    """Test JSON Input/Output."""

    def test_save_object(self):

        # Make test data
        testdata = SomeMultiLayerObject(
            a=72661, 
            b=SomeDetail(onething='some information', anotherthing=88), 
            c='more info')
        
        # Save as dictionary
        result = save_object(testdata)

        # Check dictionary as expected
        self.assertEqual(
            result,
            { 'a': 72661,
              'b': { 
                    'onething': 'some information',
                    'anotherthing': 88,
                    'python_class': 'eumopps.jsonobjects.test.test_jsonobjects.SomeDetail' },
              'c': 'more info',
              'python_class': 'eumopps.jsonobjects.test.test_jsonobjects.SomeMultiLayerObject' })


    def test_load_object(self):

        # Make test data (dictionary)
        testdata = {
            'python_class': 'eumopps.jsonobjects.test.test_jsonobjects.SomeMultiLayerObject',
            'a': 'apples',
            'b': 'bananas',
            'c': {
                'python_class': 'eumopps.jsonobjects.test.test_jsonobjects.SomeDetail',
                'onething': 'a fish',
                'anotherthing': 'not a fruit' }
            }

        # Convert to classes
        result = load_object(testdata)

        # Check it worked
        self.assertIsInstance(result, SomeMultiLayerObject)
        self.assertIsInstance(result.a, str)
        self.assertIsInstance(result.b, str)
        self.assertIsInstance(result.c, SomeDetail)
        self.assertEqual(result.c.onething, 'a fish')
        self.assertEqual(result.c.anotherthing, 'not a fruit')
        

    def test_save(self):
        
        # Make test data
        testdata = SomeMultiLayerObject(
            a=72661, 
            b=SomeDetail(onething='some information', anotherthing=88), 
            c='more info')
        
        # Save to memory stream
        outputstream = StringIO()
        save(outputstream, testdata)

        # Check as expected
        self.assertEqual(
            outputstream.getvalue(),
            '{"a": 72661, "c": "more info", "b": {"onething": "some information", "python_class": "eumopps.jsonobjects.test.test_jsonobjects.SomeDetail", "anotherthing": 88}, "python_class": "eumopps.jsonobjects.test.test_jsonobjects.SomeMultiLayerObject"}')

    def test_load(self):
        
        # Make test data string
        testdata = '{ "python_class": "eumopps.jsonobjects.test.test_jsonobjects.SomeMultiLayerObject", "a": "apples", "b": "bananas", "c": { "python_class": "eumopps.jsonobjects.test.test_jsonobjects.SomeDetail", "onething": "a fish", "anotherthing": "not a fruit" } }'

        # Convert to classes via input stream
        inputstream = StringIO(testdata)
        result = load(inputstream)

        # Check it worked
        self.assertIsInstance(result, SomeMultiLayerObject)
        self.assertIsInstance(result.a, str)
        self.assertIsInstance(result.b, str)
        self.assertIsInstance(result.c, SomeDetail)
        self.assertEqual(result.c.onething, 'a fish')
        self.assertEqual(result.c.anotherthing, 'not a fruit')

class TestDateTimeNumeric(unittest.TestCase):

    def test_save(self):

        testobject = SomeDetail(datetime(2012, 3, 7, 12, 05, 32), None)
        outputstream = StringIO()
        save(outputstream, testobject)
        self.assertEqual('{"onething": {"timestring": "20120307120532", "python_class": "datetime.datetime"}, "python_class": "eumopps.jsonobjects.test.test_jsonobjects.SomeDetail"}', outputstream.getvalue())

    def test_load(self):

        teststring = '{ "python_class": "eumopps.jsonobjects.test.test_jsonobjects.SomeDetail", "onething": {"python_class": "datetime.datetime", "timestring": "20120307120532"}, "anotherthing": 22 }'
        teststream = StringIO(teststring)
        result = load(teststream)
        self.assertIsInstance(result, SomeDetail)
        self.assertIsInstance(result.onething, datetime)
        self.assertEqual(datetime(2012, 3, 7, 12, 05, 32), result.onething)
        
