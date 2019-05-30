"""Processing pipeline for binning of satellite data."""

import numpy
import sys
from iris import load_cube
from iris import Constraint
from satellite import SatelliteCollection
from gridbin.gridfieldlist import GridFieldAxis
from gridbin.gridfieldlist import GridFieldDescriptor
from gridbin.gridbin import GridBins
from gridbin.operation.count import GridBinOperationCount
from gridbin.operation.sumcountmaxmin import GridBinOperationSumCountMaxMin
from gridbin.operation.sumsq import GridBinOperationSumSq
from gridbin.operation.sumsq import GridBinOperationSumConstSq
from gridbin.operation.sumsqdev import GridBinOperationSumSqDev

class GridBinPipelineLoadResult(object):
    """Represent result of loading a file."""

    STATUS_OK = 'ok'
    STATUS_MISSING = 'missing'
    STATUS_ERROR = 'error'

    def __init__(self, filename, status, message):
        self.filename = filename
        self.status = status
        self.message = message

class GridBinPipelineDiagnostics(object):
    """Represent results of running aggregation."""

    def __init__(self):
        self.load = [ ]

class PostProcessMean(object):
    """Compute mean from aggregated sum and count."""

    def __init__(self, fieldoutput, var_name_count, var_name_sum):
        self.fieldoutput = fieldoutput
        self.var_name_count = var_name_count
        self.var_name_sum = var_name_sum
        
    def run(self, bins):
        countvalid = bins.get_field(self.var_name_count).data
        meanfield = (bins.get_field(self.var_name_sum).data / countvalid).astype(self.fieldoutput.dtype)
        bins.create_field(self.fieldoutput, meanfield)

class PostProcessRandomError(object):
    """Compute mean random error as sum of squared errors divided by sum of squared count."""

    def __init__(self, fieldoutput, var_name_count, var_name_sumsq):
        self.fieldoutput = fieldoutput
        self.var_name_count = var_name_count
        self.var_name_sumsq = var_name_sumsq
        
    def run(self, bins):
        countvalid = bins.get_field(self.var_name_count).data
        errorfield = numpy.ma.sqrt(bins.get_field(self.var_name_sumsq).data / (countvalid * countvalid), dtype=self.fieldoutput.dtype)
        bins.create_field(self.fieldoutput, errorfield)

class PostProcessRMS(object):
    """Compute RMS error from sum of squared errors and count."""

    def __init__(self, fieldoutput, var_name_count, var_name_sumsq):
        self.fieldoutput = fieldoutput
        self.var_name_count = var_name_count
        self.var_name_sumsq = var_name_sumsq
        
    def run(self, bins):
        countvalid = bins.get_field(self.var_name_count).data
        rmsfield = numpy.ma.sqrt(bins.get_field(self.var_name_sumsq).data / countvalid, dtype=self.fieldoutput.dtype)
        bins.create_field(self.fieldoutput, rmsfield)

class PostProcessSamplingUncertainty(object):
    """Compute sampling uncertainty based on total possible samples, total valid samples, and sum of squared deviation from mean."""

    def __init__(self, fieldvariance, fielduncspl, var_name_countobs, var_name_countvalid, var_name_sumsqdev):
        self.fieldvariance = fieldvariance
        self.fielduncspl = fielduncspl
        self.var_name_countobs = var_name_countobs
        self.var_name_countvalid = var_name_countvalid
        self.var_name_sumsqdev = var_name_sumsqdev

    def run(self, bins):
        ntotal = bins.get_field(self.var_name_countobs).data
        nvalid = bins.get_field(self.var_name_countvalid).data
        sumsqdev = bins.get_field(self.var_name_sumsqdev).data
        variance = sumsqdev / (nvalid - 1)
        sampling_uncertainty = ((ntotal - nvalid) * numpy.ma.sqrt(variance)) / (ntotal - 1)
        bins.create_field(self.fieldvariance, variance.astype(self.fieldvariance.dtype))
        bins.create_field(self.fielduncspl, sampling_uncertainty.astype(self.fielduncspl.dtype))

class GridBinPipelineSatellite(object):
    """Processing pipeline for MODIS inputs."""

    @staticmethod
    def from_command_arguments(args):
        """Build pipeline from given arguments."""

        return GridBinPipelineSatellite(
            axis_lat=GridFieldAxis(args.y0-0.5*args.ys, args.ys, args.yn, 'latitude', args.wrapcoords,
                                 units='degrees_north', standard_name='latitude', long_name='Latitude (deg)'),
            axis_lon=GridFieldAxis(args.x0-0.5*args.xs, args.xs, args.xn, 'longitude', args.wrapcoords,
                                 units='degrees_east', standard_name='longitude', long_name='Longitude (deg)'),
            qc_mask_obs=args.qc_mask_obs, qc_filter_obs=args.qc_filter_obs,
            qc_mask_valid=args.qc_mask_valid, qc_filter_valid=args.qc_filter_valid)


    def __init__(self, axis_lat, axis_lon, qc_mask_obs, qc_filter_obs, qc_mask_valid, qc_filter_valid):
        """Build grid onto which observations can be binned."""

        self.bins = GridBins([axis_lat, axis_lon])
        self.qc_mask_obs = qc_mask_obs
        self.qc_filter_obs = qc_filter_obs
        self.qc_mask_valid = qc_mask_valid
        self.qc_filter_valid = qc_filter_valid

    def run(self, fieldnames, inputlist):
        """Run this pipeline on the given inputs."""

        # Diagnostics instance
        diagnostics = GridBinPipelineDiagnostics()

        # Initial operation to count valid observation locations
        # needed for sampling uncertainty
        obscounter = GridBinOperationCount() \
                         .set_input_descriptor(GridBinOperationCount.INPUT, GridFieldDescriptor(fieldnames.primary)) \
                         .set_output_descriptor(GridBinOperationCount.OUTPUT, GridFieldDescriptor('total_number_of_observations', dtype=numpy.int32, usemask=False)) 

        # Operations (first batch)
        opsA = [ ]

        # Post processing (first batch)
        postprocA = [ ]

        # Operations (second batch)
        opsB = [ ]

        # Post processing (second batch)
        postprocB = [ ]

        # Sum and count
        opsA.append(GridBinOperationSumCountMaxMin() \
                       .set_input_descriptor(GridBinOperationSumCountMaxMin.INPUT, GridFieldDescriptor(fieldnames.primary)) \
                       .set_output_descriptor(GridBinOperationSumCountMaxMin.OUTPUTTOTAL, GridFieldDescriptor(fieldnames.output_primary+'sum', dtype=numpy.float32)) \
                       .set_output_descriptor(GridBinOperationSumCountMaxMin.OUTPUTCOUNT, GridFieldDescriptor(fieldnames.output_primary+'_number_of_observations', dtype=numpy.int32, usemask=False)) \
                       .set_output_descriptor(GridBinOperationSumCountMaxMin.OUTPUTMAX, GridFieldDescriptor(fieldnames.output_primary+'max', dtype=numpy.float32), initialvalue=-1.0E20) \
                       .set_output_descriptor(GridBinOperationSumCountMaxMin.OUTPUTMIN, GridFieldDescriptor(fieldnames.output_primary+'min', dtype=numpy.float32), initialvalue=1.0E20))

        # Compute sum squares for each uncertainty field
        for index, inputname in enumerate(fieldnames.uncertainty_fields):

            # Variable name for this output
            var_name_sumsq = fieldnames.output_uncertainty[index] + '_sumsq'

            # Sum squares
            opsA.append(GridBinOperationSumSq() \
                           .set_input_descriptor(GridBinOperationSumSq.INPUT, GridFieldDescriptor(inputname)) \
                           .set_output_descriptor(GridBinOperationSumSq.OUTPUT, GridFieldDescriptor(var_name_sumsq, dtype=numpy.float32)))

        for index, inputname in enumerate(fieldnames.uncertainty_scalars):
                
            # Variable name for this output
            var_name_sumsq = fieldnames.output_uncertainty_binned_scalars[index] + '_sumsq'

            opsA.append(GridBinOperationSumConstSq() \
                            .set_input_descriptor(GridBinOperationSumConstSq.INPUT, GridFieldDescriptor(inputname)) \
                            .set_output_descriptor(GridBinOperationSumConstSq.OUTPUT, GridFieldDescriptor(var_name_sumsq, dtype=numpy.float32)))
            
        # Mean from sum and count
        postprocA.append(PostProcessMean(GridFieldDescriptor(fieldnames.output_primary+'mean', dtype=numpy.float32), var_name_count=fieldnames.output_primary+'_number_of_observations', var_name_sum=fieldnames.output_primary+'sum'))

        # RMS computation for each uncertainty component
        for index, outputname in enumerate(fieldnames.output_uncertainty):

            # Variable name for this output
            var_name_sumsq = outputname + '_sumsq'

            if index == 0:                
                # Random error means divide by square of count
                postprocA.append(PostProcessRandomError(GridFieldDescriptor(outputname, dtype=numpy.float32), var_name_count=fieldnames.output_primary+'_number_of_observations', var_name_sumsq=var_name_sumsq))
            else:
                # Use RMS
                postprocA.append(PostProcessRMS(GridFieldDescriptor(outputname, dtype=numpy.float32), var_name_count=fieldnames.output_primary+'_number_of_observations', var_name_sumsq=var_name_sumsq))

        # RMS computation for each (binned scalar) uncertainty component
        for outputname in fieldnames.output_uncertainty_binned_scalars:

            # Variable name for this output
            var_name_sumsq = outputname + '_sumsq'

            # RMS
            postprocA.append(PostProcessRMS(GridFieldDescriptor(outputname, dtype=numpy.float32), var_name_count=fieldnames.output_primary+'_number_of_observations', var_name_sumsq=var_name_sumsq))
            

        # Sum squared deviation from the mean computed above
        opsB.append(GridBinOperationSumSqDev() \
                       .set_input_descriptor(GridBinOperationSumSqDev.INPUTVALUE, GridFieldDescriptor(fieldnames.primary)) \
                       .set_output_descriptor(GridBinOperationSumSqDev.EXISTINGMEAN, GridFieldDescriptor(fieldnames.output_primary+'mean')) \
                       .set_output_descriptor(GridBinOperationSumSqDev.OUTPUT, GridFieldDescriptor(fieldnames.output_primary+'sumsqdev', dtype=numpy.float32)))

        # Post-processing is computation of sampling uncertainty from coverage information
        postprocB.append(PostProcessSamplingUncertainty(
                fieldvariance=GridFieldDescriptor(fieldnames.output_primary+'variance', dtype=numpy.float32),
                fielduncspl=GridFieldDescriptor(fieldnames.output_primary+'mean_unc_spl', dtype=numpy.float32),
                var_name_countobs='total_number_of_observations',
                var_name_countvalid=fieldnames.output_primary+'_number_of_observations',
                var_name_sumsqdev=fieldnames.output_primary+'sumsqdev'))

        # List of variables to remove afterwards (intermediate only)
        intermediate_only = \
            [ fieldnames.output_primary+'sum', fieldnames.output_primary+'sumsqdev' ] + \
            [ outputname + '_sumsq' for outputname in fieldnames.output_uncertainty ] + \
            [ outputname + '_sumsq' for outputname in fieldnames.output_uncertainty_binned_scalars ]

        # List of filenames to use
        sublist = inputlist.filenames

        # Build list of input sources
        # for filegroup in inputlist.filenames:
        for inputindex, filegroup in enumerate(sublist):

            # No available input yet
            satelliteinput = None

            # Check this file available
            if filegroup.check_readable():

                # diagnostic info
                # sys.stderr.write('satgrid pipeline operation A {0}/{1}\n'.format(inputindex, len(sublist)))
                # sys.stderr.write(filegroup[0] + '\n')
                # sys.stderr.write(filegroup[1] + '\n')

                try:

                    # Load this item
                    satelliteinput = SatelliteCollection.load(fieldnames, filegroup)
                                    
                    # Loaded ok if we get here without exception thrown
                    diagnostics.load.append(GridBinPipelineLoadResult(filegroup, GridBinPipelineLoadResult.STATUS_OK, ''))

                except Exception as result:
            
                    # Error during load - report now
                    diagnostics.load.append(GridBinPipelineLoadResult(filegroup, GridBinPipelineLoadResult.STATUS_ERROR, str(result)))
            else:
                
                # Report file(s) not found
                diagnostics.load.append(GridBinPipelineLoadResult(filegroup, GridBinPipelineLoadResult.STATUS_MISSING, ''))


            if satelliteinput:

		# Generate safety mask for excluding coordinates with invalid values
		satelliteinput.generate_safety_mask()

                # Apply filter so that we are looking only at valid observations
                satelliteinput.apply_filter(satelliteinput.get_filter_from_qc_flags(self.qc_mask_valid, self.qc_filter_valid))

                # Count field of observation locations (required for sampling uncertainty)
                obsmap = self.bins.compute_input_map(satelliteinput.cube, alternative_mask=satelliteinput.get_filter_from_qc_flags(self.qc_mask_obs, self.qc_filter_obs), input_coordinate_names=fieldnames.coordinate_fields)

                # Map for applying operation to valid observations
                inputmap = self.bins.compute_input_map(satelliteinput.cube, input_coordinate_names=fieldnames.coordinate_fields)

                # Diagnose a problem noticed in grid box number 527905
                #problemindices = numpy.nonzero(inputmap.mapto == 527905)
                #print 'indices'
                #print problemindices
                #print 'primary mask'
                #print satelliteinput.cube.data.mask.flat[inputmap.mapfrom[problemindices]]
                #print 'sfc mask'
                #sfc = load_cube(filegroup[1], constraint=Constraint(cube_func=lambda cube: cube.var_name=='LST_unc_loc_sfc'))
                #print sfc.data.mask.flat[inputmap.mapfrom[problemindices]]
                #print 'sfc data'
                #print sfc.data.data.flat[inputmap.mapfrom[problemindices]]
                #print 'atm mask'
                #atm = load_cube(filegroup[1], constraint=Constraint(cube_func=lambda cube: cube.var_name=='LST_unc_loc_atm'))
                #print atm.data.mask.flat[inputmap.mapfrom[problemindices]]
                #print 'atm data'
                #print atm.data.data.flat[inputmap.mapfrom[problemindices]]

                # Count observation locations (required for sampling uncertainty)
                obscounter.run(obsmap)

                # Apply operations
                for op in opsA:
                    op.run(inputmap)                


        # Finish now if nothing to load
        if not any([ loadresult.status == GridBinPipelineLoadResult.STATUS_OK for loadresult in diagnostics.load ]):
            return diagnostics

        # First batch of postprocessing
        for proc in postprocA:
            proc.run(self.bins)

        # Run second batch of operations
        for inputindex, loadresult in enumerate(diagnostics.load):

            if loadresult.status == GridBinPipelineLoadResult.STATUS_OK:

                # Re-load and prepare this file
                filegroup = sublist[inputindex]

                # diagnostic info
                # sys.stderr.write('satgrid pipeline operation B {0}/{1}\n'.format(inputindex, len(diagnostics.load)))
                # sys.stderr.write(filegroup[0] + '\n')
                # sys.stderr.write(filegroup[1] + '\n')

                # load and prepare
                satelliteinput = SatelliteCollection.load(fieldnames, filegroup)

		# Generate safety mask for excluding coordinates with invalid values
		satelliteinput.generate_safety_mask()

                satelliteinput.apply_filter(satelliteinput.get_filter_from_qc_flags(self.qc_mask_valid, self.qc_filter_valid))
                inputmap = self.bins.compute_input_map(satelliteinput.cube, input_coordinate_names=fieldnames.coordinate_fields)
                
                # Run second round of operations
                for op in opsB:
                    op.run(inputmap)

        # Second batch of postprocessing
        for proc in postprocB:
            proc.run(self.bins)

        # Remove intermediate items
        indices_to_remove = [ index for index, cube in enumerate(self.bins) if cube.var_name in intermediate_only ]
        for index in reversed(indices_to_remove):
            del self.bins[index]

        # Application-specific attributes to set
        attributes_to_set = {
            'ts_number_of_observations': { 'units': '1', 'long_name': 'Number of observations of surface temperature' },
            'total_number_of_observations': { 'units': '1', 'long_name': 'Total number of satellite observations in grid box' },
            'tsmean': { 'units': 'K', 'long_name': 'Mean surface temperature (K)', 'standard_name': 'surface_temperature' },
            'tsmax': { 'units': 'K', 'long_name': 'Maximum surface temperature (K)', 'standard_name': 'surface_temperature' },
            'tsmin': { 'units': 'K', 'long_name': 'Minimum surface temperature (K)', 'standard_name': 'surface_temperature' },
            'tsvariance': { 'units': 'K2', 'long_name': 'Surface temperature variance (K2)' },
            'tsmean_unc_ran': { 'units': 'K', 'long_name': 'Random uncertainty in mean surface temperature (K)' },
            'tsmean_unc_loc_atm': { 'units': 'K', 'long_name': 'Locally correlated uncertainty (atm component) in mean surface temperature (K)' },
            'tsmean_unc_loc_sfc': { 'units': 'K', 'long_name': 'Locally correlated uncertainty (sfc component) in mean surface temperature (K)' },
            'tsmean_unc_spl': { 'units': 'K', 'long_name': 'Sampling uncertainty in mean surface temperature (K)' },
            'tsmean_unc_sys': { 'units': 'K', 'long_name': 'Systematic uncertainty in mean surface temperature (K)' },
            }

        # Set the attributes
        for (var_name, attributes) in attributes_to_set.iteritems():
            field = self.bins.get_field(var_name)
            if field:
                for (attr_name, attr_value) in attributes.iteritems():
                    setattr(field, attr_name, attr_value)

        # return diagnostics object
        return diagnostics
