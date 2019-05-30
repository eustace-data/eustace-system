import unittest
import numpy
import apply_mask

class TestAnalysisSystem(unittest.TestCase):

    def test_compute_flags(self):

        # flag definitions
        flag_type = numpy.int16
        
        flag_0 = flag_type(1) << 0
        flag_1 = flag_type(1) << 1
        flag_2 = flag_type(1) << 2
        flag_3 = flag_type(1) << 3
        
        null_flag = flag_type(0)
        
        # the data flags that will be combined to produce a mask
        qc_flags = numpy.array([[flag_0,          flag_1,          flag_2,          flag_3,                              null_flag,],
                                [flag_0 | flag_1, flag_0 | flag_2, flag_0 | flag_3, flag_0 | flag_1 | flag_2 | flag_3,  flag_0 | flag_1 | flag_3]], flag_type)
        
        
        # mask on any flag_1 being set
        flags_to_mask = numpy.array([ [flag_1,    null_flag,],
                                     ], flag_type )
        
        mask = apply_mask.compute_mask(qc_flags, flags_to_mask)
        
        expected = numpy.array([[False, True, False, False, False,],
                                [True, False, False, True, True,]])
        
        numpy.testing.assert_array_equal(mask, expected)
        
        # mask on any flag being set
        flags_to_mask = numpy.array([ [flag_0,    null_flag,],
                                      [flag_1,    null_flag,],
                                      [flag_2,    null_flag,],
                                      [flag_3,    null_flag,],
                                     ], flag_type )
        
        mask = apply_mask.compute_mask(qc_flags, flags_to_mask)
        
        expected = numpy.array([[True, True, True, True, False,],
                                [True, True, True, True, True,]])
        
        numpy.testing.assert_array_equal(mask, expected)
        
        # mask on all of flag_0, flag_1 and flag_3 being set, but clear the mask if flag_2 is set
        flags_to_mask = numpy.array([ [flag_0 | flag_1 | flag_3,    flag_2,],
                                     ], flag_type )
        
        mask = apply_mask.compute_mask(qc_flags, flags_to_mask)
        
        expected = numpy.array([[False, False, False, False, False,],
                                [False, False, False, False, True,]])
                                
        numpy.testing.assert_array_equal(mask, expected)
        
        # mask on any of flag_0, flag_1 and flag_3 being set, but override each individually if flag_2 is set
        flags_to_mask = numpy.array([ [flag_0,    flag_2,],
                                      [flag_1,    flag_2,],
                                      [flag_3,    flag_2,],
                                     ], flag_type )
        
        mask = apply_mask.compute_mask(qc_flags, flags_to_mask)
        
        expected = numpy.array([[True, True, False, True, False,],
                                [True, False, True, False, True,]])
        
        numpy.testing.assert_array_equal(mask, expected)
        
        # mask on flag_1 and flag_3. Only mask on flag_0 if flag 2 is not set. Note that flag_2 being set only overrides flag_0. Other flags are not affected by the override.
        flags_to_mask = numpy.array([ [flag_0,    flag_2,],
                                      [flag_1,    null_flag,],
                                      [flag_3,    null_flag,],
                                     ], flag_type )
        
        mask = apply_mask.compute_mask(qc_flags, flags_to_mask)
        
        expected = numpy.array([[True, True, False, True, False,],
                                [True, False, True, True, True,]])
                                
        numpy.testing.assert_array_equal(mask, expected)
        
        # mask on flag_0 and flag_2. Override flag_0 if flag_2 is set. This isn't a sensible flagging setup as the flag_2 will still result in masking where both flag_0 and flag_2 are set, regardless of the override for flag_0.
        flags_to_mask = numpy.array([ [flag_0,    flag_2,],
                                      [flag_2,    null_flag,],
                                     ], flag_type )
        
        mask = apply_mask.compute_mask(qc_flags, flags_to_mask)
        
        expected = numpy.array([[True, False, True, False, False,],
                                [True, True, True, True, True,]])
                                
        numpy.testing.assert_array_equal(mask, expected)
