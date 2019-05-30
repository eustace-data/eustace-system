"""Tests for sphere mesh SPDE and original experimental display code."""

import unittest

from eustace.analysis.advanced_standard.stats.spde.spherical import SphereMeshSPDE
from eustace.analysis.advanced_standard.stats.spde.spherical_view import SphereMeshView, SphereMeshViewGlobal, SphereMeshViewSuperTriangle, SphereMeshViewLocal
import numpy

class TestSphereMeshView(unittest.TestCase):

    def test_derivative(self):
        """Check that numerical and analytical derivatives of Q match."""
    
        full_resolution_level = 4
        neighbourhood_level = 2
    
        full_spde = SphereMeshSPDE(level=full_resolution_level)
        
        active_triangles = full_spde.neighbours_at_level(neighbourhood_level, 0)
        
        spde = SphereMeshView(full_resolution_level, active_triangles)
        
        Q = spde.build_Q_stationary(numpy.log(1.0), numpy.log(1.0), 2)

        epsilon = 0.000001

        expected_dQ0 = \
            ((spde.build_Q_stationary(numpy.log(1.0)+epsilon, numpy.log(1.0), 2) - 
              spde.build_Q_stationary(numpy.log(1.0)-epsilon, numpy.log(1.0), 2)) / (2.0 * epsilon)).todense()

        expected_dQ1 = \
            ((spde.build_Q_stationary(numpy.log(1.0), numpy.log(1.0)+epsilon, 2) - 
              spde.build_Q_stationary(numpy.log(1.0), numpy.log(1.0)-epsilon, 2)) / (2.0 * epsilon)).todense()

        dQ0 = spde.build_dQdp_stationary(numpy.log(1.0), numpy.log(1.0), 2, 0).todense()
        dQ1 = spde.build_dQdp_stationary(numpy.log(1.0), numpy.log(1.0), 2, 1).todense()

        # print numpy.abs(dQ0 - expected_dQ0).ravel().max() / numpy.abs(expected_dQ0).ravel().max()
        # print numpy.abs(dQ1 - expected_dQ1).ravel().max() / numpy.abs(expected_dQ1).ravel().max()
        
        numpy.testing.assert_almost_equal(dQ0, expected_dQ0, decimal=7)
        numpy.testing.assert_almost_equal(dQ1, expected_dQ1, decimal=7)
        
    def test_non_stationary_merge(self):
        """Simple test of hyperparameter merging for building a nonstationary model"""
        
        full_resolution_level = 1
        neighbourhood_level = 0
    
        full_spde = SphereMeshViewGlobal(level=full_resolution_level)
        
        active_triangles = full_spde.neighbours_at_level(neighbourhood_level, 0)
        
        n_regions = full_spde.n_triangles_at_level(neighbourhood_level)
    
        sigma_accumulator = None
        rho_accumulator = None
        contribution_counter = None
        
        hyperparameters = numpy.log( numpy.concatenate( [numpy.random.uniform(1.0,5.0, 1),
                                                         numpy.random.uniform(10.0,45.0, 1) *numpy.pi/180.] ) )
        
        for region_index in range(n_regions):
            local_spde = SphereMeshViewSuperTriangle(full_resolution_level, neighbourhood_level, region_index)
            local_hyperparameters = hyperparameters
            
            accumulators = SphereMeshViewGlobal.accumulate_local_parameterisations(sigma_accumulator,
                                                                                    rho_accumulator,
                                                                                    contribution_counter,
                                                                                    local_spde,
                                                                                    local_hyperparameters)
                                                                                                                                            
            sigma_accumulator, rho_accumulator, contribution_counter = accumulators
            
    
        log_sigmas, log_rhos = SphereMeshViewGlobal.finalise_local_parameterisation_sigma_rho(sigma_accumulator,
                                                                                                rho_accumulator,
                                                                                                contribution_counter)
        
        numpy.testing.assert_almost_equal( log_sigmas, hyperparameters[0] * numpy.ones(full_spde.triangulation.points.shape[0]) )
        numpy.testing.assert_almost_equal( log_rhos, hyperparameters[1] * numpy.ones(full_spde.triangulation.points.shape[0]) )
        
        self.assertEqual( full_spde.n_latent_variables(), len(log_sigmas))
        self.assertEqual( full_spde.n_latent_variables(), len(log_rhos))
        
def demo_non_stationary():
    
    full_resolution_level = 5
    neighbourhood_level = 2

    full_spde = SphereMeshViewGlobal(level=full_resolution_level)
    
    active_triangles = full_spde.neighbours_at_level(neighbourhood_level, 0)
    
    n_regions = full_spde.n_triangles_at_level(neighbourhood_level)
    
    merge_method = 'new'
    if merge_method == 'old':
        local_spdes = []
        local_hyperparameters = []
        
        for region_index in range(n_regions):
            
            local_spdes.append( SphereMeshViewSuperTriangle(full_resolution_level, neighbourhood_level, region_index) )
            
            hyperparameters = numpy.array([numpy.float64(region_index), numpy.float64(region_index)])
            hyperparameters = numpy.log( numpy.concatenate( [numpy.random.uniform(1.0,3.0, 1), numpy.random.uniform(5.0,30.0, 1) *numpy.pi/180.] ) )
            #hyperparameters = numpy.log( numpy.concatenate( [numpy.ones(1), numpy.random.uniform(15.0,45.0, 1) *numpy.pi/180.] ) )
            #hyperparameters = numpy.array([2.0, 3.0])
            
            #hyperparameters = numpy.log([2.0, numpy.pi/4])
            
            local_hyperparameters.append( hyperparameters  )
            
        global_hyperparameters, global_sigma_design, global_rho_design = full_spde.merge_local_parameterisations( local_spdes, local_hyperparameters, merge_method = 'exp_average' )
        
        log_sigmas = global_sigma_design.dot(global_hyperparameters)
        log_rhos   = global_rho_design.dot(global_hyperparameters)
        
    elif merge_method == 'new':

        sigma_accumulator = None
        rho_accumulator = None
        contribution_counter = None
        
        for region_index in range(n_regions):
            local_spde = SphereMeshViewSuperTriangle(full_resolution_level, neighbourhood_level, region_index)
            local_hyperparameters = hyperparameters = numpy.log( numpy.concatenate( [numpy.random.uniform(1.0,5.0, 1),
                                                                                        numpy.random.uniform(10.0,45.0, 1) *numpy.pi/180.] ) )
            
            accumulators = SphereMeshViewGlobal.accumulate_local_parameterisations(sigma_accumulator,
                                                                                    rho_accumulator,
                                                                                    contribution_counter,
                                                                                    local_spde,
                                                                                    local_hyperparameters)
                                                                                                                                            
            sigma_accumulator, rho_accumulator, contribution_counter = accumulators
            
    
        log_sigmas, log_rhos = SphereMeshViewGlobal.finalise_local_parameterisation_sigma_rho(sigma_accumulator,
                                                                                                rho_accumulator,
                                                                                                contribution_counter)
        
    
    
    
    
    #print global_hyperparameters, global_sigma_design, global_rho_design
    
    import matplotlib.pyplot as plt
    from eustace.analysis.mesh.geometry import cartesian_to_polar2d
    polar_coords = cartesian_to_polar2d( full_spde.triangulation.points )         
    
    plt.figure()
    plt.scatter(polar_coords[:,1], polar_coords[:,0], c = 255.* log_sigmas / numpy.max(numpy.abs(log_sigmas)), linewidth = 0.0, s = 8.0 )
    
    plt.figure()
    plt.scatter(polar_coords[:,1], polar_coords[:,0], c = 255.* log_rhos / numpy.max(numpy.abs(log_rhos)), linewidth = 0.0, s = 8.0 )
    
    #plt.show()
    
    #numpy.testing.assert_almost_equal( log_sigmas, 2.0 * numpy.ones(full_spde.triangulation.points.shape[0]) )
    #numpy.testing.assert_almost_equal( log_rhos, 3.0 * numpy.ones(full_spde.triangulation.points.shape[0]) )
    
    from eustace.analysis.advanced_standard.components.storage_inmemory import ComponentStorage_InMemory
    from eustace.analysis.advanced_standard.components.storage_inmemory import SpatialComponentSolutionStorage_InMemory
    from eustace.analysis.advanced_standard.components.spatialdelayed import DelayedSpatialComponent
    from eustace.analysis.advanced_standard.elements.local_view import NonStationaryLocal, ExpandedLocalHyperparameters
    from eustace.analysis.advanced_standard.elements.local import LocalElement, LocalHyperparameters
    
    nonstationary_component = DelayedSpatialComponent(
            ComponentStorage_InMemory(NonStationaryLocal(full_resolution_level), ExpandedLocalHyperparameters(log_sigma = log_sigmas, log_rho = log_rhos)),
            SpatialComponentSolutionStorage_InMemory())
    
    #nonstationary_component = DelayedSpatialComponent(
            #ComponentStorage_InMemory(LocalElement(full_resolution_level), LocalHyperparameters(log_sigma = hyperparameters[0], log_rho = hyperparameters[1])),
            #SpatialComponentSolutionStorage_InMemory())
    
    #print log_sigmas, log_rhos
    
    #plt.figure()
    #plt.scatter(polar_coords[:,1], polar_coords[:,0], c = 255.* process_sample / numpy.max(numpy.abs(process_sample)), linewidth = 0.0, s = 8.0 )
    
    #plt.figure()
    #plt.imshow( numpy.asarray( Q.todense() ) )
    
    # setup an output grid
    out_lats = numpy.linspace(-89.5, 89.5, 180 )
    out_lons = numpy.linspace(-179.5, 179.5, 360 )
    out_lons, out_lats = numpy.meshgrid( out_lons, out_lats )
    out_coords = numpy.vstack( [out_lats.ravel(), out_lons.ravel()] ).T
    
    design_matrix = nonstationary_component.storage.element.spde.build_A( out_coords )
    
    # setup solver for sampling
    from eustace.analysis.advanced_standard.linalg.extendedcholmodwrapper import ExtendedCholmodWrapper
    Q = nonstationary_component.storage.element.element_prior(nonstationary_component.storage.hyperparameters).prior_precision()
    factor = ExtendedCholmodWrapper.cholesky( Q )
    
    # draw samples, project onto output grid and plot
    random_values = numpy.random.normal(0.0,1.0, (Q.shape[0],1))        
    process_sample = factor.solve_backward_substitution(random_values)        
    out_values = design_matrix.dot( process_sample )
    plt.figure()
    plt.scatter(out_coords[:,1], out_coords[:,0], c = 255.* out_values / numpy.max(numpy.abs(out_values)), linewidth = 0.0, s = 8.0 )
    
    random_values = numpy.random.normal(0.0,1.0, (Q.shape[0],1))        
    process_sample = factor.solve_backward_substitution(random_values)        
    out_values = design_matrix.dot( process_sample )
    plt.figure()
    plt.scatter(out_coords[:,1], out_coords[:,0], c = 255.* out_values / numpy.max(numpy.abs(out_values)), linewidth = 0.0, s = 8.0 )
    
    random_values = numpy.random.normal(0.0,1.0, (Q.shape[0],1))        
    process_sample = factor.solve_backward_substitution(random_values)        
    out_values = design_matrix.dot( process_sample )
    plt.figure()
    plt.scatter(out_coords[:,1], out_coords[:,0], c = 255.* out_values / numpy.max(numpy.abs(out_values)), linewidth = 0.0, s = 8.0 )
    
    
    plt.show()        
        