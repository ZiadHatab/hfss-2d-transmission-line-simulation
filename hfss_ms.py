"""
Author: @Ziad (https://github.com/ZiadHatab)

A code to run HFSS to solve a 2D full-wave microstrip line Z0 and gamma

NOTE:
    - All material parameters can be frequency-dependent
    - Frequency sweep is only solved at provided frequency points. 
    You need to interpolate the results yourself. Tip: use logspace for the frequency grid.
    
TODO:
    - (DONE!) Use H-symmetry boundary to half the simulation
    - (DOESN'T WORK!) Enable analytical derivative computation
    - (DONE!) Allow dielectric properties to be frequency-dependent
    - (DONE!) Add surface roughness through Surface Impedance
    - (DONE!) Add platting through Surface Impedance
    - (In progress) Clean up code
"""

import pyaedt
import numpy as np
import surfz  # my code for computing surface impedance. Should be in same folder as this script

class HFSSMS:
    """
    Run AEDT HFSS in headless (or graphic) mode to simulate a microstrip line 
    The script gives back the cross-section parameters, Z0 and gamma, 
    as well as (when enabled) their Jacobian with respect to input parameters (TODO).
    
    ## ===============================================             
    ##                    <-- w -->                                
    ##                 ++++++++++++++ t                 Signal line
    ## -----------------------------------------------             
    ## |                                             ^             
    ## |                                             h  Substrate  
    ## |                                             v                                                     
    ## -----------------------------------------------             
    ## -----------------------------------------------  GND         
    ##                  <-- wgnd -->                               
    ## ===============================================             

    Parameters
    ----------
    f : number or 1d-array
        frequency in Hz.
    w : number
        width of signal line in meters.
    h : number
        height of dielectric substrate in meters.
    t : number
        thickness of conductor in meters. Same used for GND (it doesn't really matter, as GND is at port boundary)
    wgnd : number
        width of GND plane in meters. Also, the width of the substrate!
    er : number or 1d-array
        dielectric relative permittivity
    etand : number or 1d-array
        dielectric loss tangent
    mur : number or 1d-array
        dielectric relative permeability
    mutand : number or 1d-array
        dielectric magnetic loss tangent
    conductor_gnd_top : list of dict or None
        surface impedance of the conductor at GND plane. If None, defaults to PEC.
    conductor_sig_bottom : list of dict or None
        surface impedance of the conductor at the bottom of signal line. If None, defaults to PEC.
    conductor_sig_side : list of dict or None
        surface impedance of the conductor at the side of signal line. If None, defaults to PEC.
    conductor_sig_top : list of dict or None
        surface impedance of the conductor at the top of signal line. If None, defaults to PEC.
    use_pec : boolean
        if you want use perfect electric conductor instead of specified values

    Attributes
    ----------
    Z0 : 1d-array
        characteristic impedance
    gamma : 1d-array
        propagation constant
    """
    def __init__(self, f, w, h, t, wgnd, er=1, etand=0, mur=1, mutand=0, 
                 conductor_gnd_top=None, 
                 conductor_sig_bottom=None,
                 conductor_sig_side=None,
                 conductor_sig_top=None,
                 use_pec=False):
        
        self.f  = np.atleast_1d(f)
        
        # geometry
        self.w    = w
        self.h    = h
        self.t    = t
        self.wgnd = wgnd

        # dielectric
        self.er    = er*np.ones_like(self.f)
        self.etand = etand*np.ones_like(self.f)
        self.mur    = mur*np.ones_like(self.f)
        self.mutand = mutand*np.ones_like(self.f)
        
        # conductor defined as through surface impedance
        # When not defined, i.e., None, then use PEC by default... Clean after DONE.
        """
        # Example:
        conductor_sig_bottom = [ {'sigma': 58e6, 'mur': 1-0j, 'Rrms': 0, 'boundary_loc': 0, 'distribution': 'norm'} ] # default dielectric on interface
        conductor_sig_side = [ {'sigma': 58e6, 'mur': 1-0j, 'Rrms': 0, 'boundary_loc': 0, 'distribution': 'norm'} ]   # default air on interface
        conductor_sig_top = [ {'sigma': 58e6, 'mur': 1-0j, 'Rrms': 0, 'boundary_loc': 0, 'distribution': 'norm'} ]    # default air on interface
        conductor_gnd_top = [ {'sigma': 58e6, 'mur': 1-0j, 'Rrms': 0, 'boundary_loc': 0, 'distribution': 'norm'} ]    # default dielectric on interface
        """
        self.conductor_gnd_top = conductor_gnd_top
        self.conductor_sig_bottom = conductor_sig_bottom
        self.conductor_sig_side = conductor_sig_side
        self.conductor_sig_top = conductor_sig_top
        self.use_pec = use_pec  # in case a PEC (perfect electric conductor) material is to be used 

        # constants
        self.mu0 = 1.25663706212e-6
        self.ep0 = 8.8541878128e-12
        self.c0  = 1/np.sqrt(self.mu0*self.ep0) # 299792458 # speed of light in vacuum (m/s)
        self.sigma_cu = 58e6
        
        # some useful inline functions
        self.mag2db = lambda x: 20*np.log10(abs(x))
        self.db2mag = lambda x: 10**(x/20)
        self.gamma2ereff = lambda x,f: -(self.c0/2/np.pi/f*x)**2
        self.ereff2gamma = lambda x,f: 2*np.pi*f/self.c0*np.sqrt(-(x-1j*np.finfo(float).eps))  # eps to ensure positive square-root
        self.gamma2dbmm  = lambda x: self.mag2db(np.exp(x.real*1e-3))  # losses dB/mm

    def fill_missing_vals(self, required_keys, default_values, material_properties):
        # fill in the missing properties not provided by the user with default values
        for material in material_properties:
            for key in required_keys:
              if key not in material:
                material[key] = default_values.get(key, None)
        return material_properties

    def run_simulation(self, closs_aedt_at_finish=True, solution_freq=None, port_accuracy=0.01, max_passes=20, headless=True):
        # solution freq if not provided
        solution_freq = self.f.max()*0.8 if solution_freq is None else solution_freq
        
        project_name = pyaedt.generate_unique_project_name(project_name='transmission_line')
        #aedt_version = "2024.2"  
        hfss = pyaedt.Hfss(project=project_name,
                #version=aedt_version,  # leave empty and automatically chosen!
                design='microstrip',
                non_graphical=headless,
                new_desktop=True,
                solution_type='Modal',
                close_on_exit=True,  # this trigger when I run the function to close AEDT
                )
        hfss.change_material_override(True)    # allows conductor to overlap with dielectric and give priorty to conductors
        hfss.change_automatically_use_causal_materials(False)  # disable HFSS non-sense with "Causal fit of dielectric materials"
        hfss.modeler.model_units = 'meter'  # let's stay metric!!
        hfss.mesh.assign_initial_mesh_from_slider(applycurvilinear=False)  # not sure if this needed?!
        
        try:
            freq = self.f.tolist()
            # create dielectric material
            hfss.create_dataset('er', x=freq, y=self.er.tolist(), 
                                is_project_dataset=True, x_unit='Hz', y_unit='')
            hfss.create_dataset('etand', x=freq, y=self.etand.tolist(), 
                                is_project_dataset=True, x_unit='Hz', y_unit='')
            hfss.create_dataset('mur', x=freq, y=self.mur.tolist(), 
                                is_project_dataset=True, x_unit='Hz', y_unit='')
            hfss.create_dataset('mutand', x=freq, y=self.mutand.tolist(), 
                                is_project_dataset=True, x_unit='Hz', y_unit='')
            dielectric = hfss.materials.add_material('my_dielectric')
            dielectric.conductivity = 0
            dielectric.permittivity = 'pwl($er, Freq)'
            dielectric.permeability = 'pwl($mur, Freq)'
            dielectric.dielectric_loss_tangent = 'pwl($etand, Freq)'
            dielectric.magnetic_loss_tangent = 'pwl($mutand, Freq)'
            dielectric.material_appearance = [64,128,128,0]
            dielectric.update()
            
            conductor_name = 'pec' #if self.use_pec else conductor.name
            
            # simulation parameters
            hfss['line_length'] = 1e-3  # this isn't important as we well simulate in 2D, it cannot be zero!
            hfss['sub_height']  = self.h
            hfss['trace_thickness'] = self.t
            hfss['sig_width'] = self.w
            hfss['sub_width'] = self.wgnd
            hfss['port_height'] = 'sub_width/2' # this is just an approximate... I probably should look into a better approach
    
            # draw the substrate 
            sub = hfss.modeler.create_box(origin=['0', '0', '0'], 
                                        sizes=['sub_width/2', 'sub_height', 'line_length'],
                                        material=dielectric.name,
                                        name='SUB')
            sub.transparency = 0
    
            # draw the signal line
            sig = hfss.modeler.create_box(origin=['0', 'sub_height', '0'], 
                                      sizes=['sig_width/2', 'trace_thickness', 'line_length'],
                                      material=conductor_name,
                                      name='SIG')
            sig.transparency = 0
            
            # draw the ground plane
            gnd = hfss.modeler.create_box(origin=['0', '0', '0'], 
                                        sizes=['sub_width/2', '-trace_thickness', 'line_length'],
                                        material=conductor_name,
                                        name='GND')
            gnd.transparency = 0
            
            if not self.use_pec:
                # common stuff for all boundaries
                required_keys = ['sigma', 'mur', 'er', 'Rrms', 'boundary_loc', 'distribution']
                default_values = {'sigma': 58e6, 'mur': 1-0j, 'er': None, 'Rrms': 0, 'boundary_loc': 0, 'distribution': 'norm'}
                x_datalist = freq
                mur_dielectric = self.mur*(1 - 1j*self.mutand)
                er_dielectric  = self.er*(1 - 1j*self.etand)
                mur_air = 1 - 0j
                er_air  = 1 - 0j
                # assign surface impedance to conductors
                if self.conductor_gnd_top is not None:
                    conductor_properties = self.fill_missing_vals(required_keys, default_values, self.conductor_gnd_top)
                    Rrms = [x['Rrms'] for x in conductor_properties]
                    boundary_loc = [x['boundary_loc'] for x in conductor_properties]
                    distribution = [x['distribution'] for x in conductor_properties]
                    material_properties = [{'mur': mur_dielectric, 'er': er_dielectric}] + conductor_properties
                    Zs_gnd_top = surfz.surface_impedance(self.f, material_properties, 
                                                         Rrms=Rrms, boundary_loc=boundary_loc, distribution=distribution)
                    # actually assign the boundary
                    y_datalist = Zs_gnd_top.real.tolist()
                    hfss.create_dataset('Zs_gnd_top_real', x_datalist, y_datalist, is_project_dataset=False, x_unit='Hz', y_unit='ohm')
                    y_datalist = Zs_gnd_top.imag.tolist()
                    hfss.create_dataset('Zs_gnd_top_imag', x_datalist, y_datalist, is_project_dataset=False, x_unit='Hz', y_unit='ohm')
                    hfss.assign_impedance_to_sheet(gnd.top_face_y, sourcename='GND-TOP', 
                                                   resistance='pwl(Zs_gnd_top_real, Freq)', reactance='pwl(Zs_gnd_top_imag, Freq)', is_infground=False)
                    
                if self.conductor_sig_top is not None:
                    conductor_properties = self.fill_missing_vals(required_keys, default_values, self.conductor_sig_top)
                    Rrms = [x['Rrms'] for x in conductor_properties]
                    boundary_loc = [x['boundary_loc'] for x in conductor_properties]
                    distribution = [x['distribution'] for x in conductor_properties]
                    material_properties = [{'mur': mur_air, 'er': er_air}] + conductor_properties
                    Zs_sig_top = surfz.surface_impedance(self.f, material_properties, 
                                                         Rrms=Rrms, boundary_loc=boundary_loc, distribution=distribution)
                    self.Zs_sig_top = Zs_sig_top
                    # actually assign the boundary
                    y_datalist = Zs_sig_top.real.tolist()
                    hfss.create_dataset('Zs_sig_top_real', x_datalist, y_datalist, is_project_dataset=False, x_unit='Hz', y_unit='ohm')
                    y_datalist = Zs_sig_top.imag.tolist()
                    hfss.create_dataset('Zs_sig_top_imag', x_datalist, y_datalist, is_project_dataset=False, x_unit='Hz', y_unit='ohm')
                    hfss.assign_impedance_to_sheet(sig.top_face_y, sourcename='SIG-TOP', 
                                                   resistance='pwl(Zs_sig_top_real, Freq)', reactance='pwl(Zs_sig_top_imag, Freq)', is_infground=False)
                        
                if self.conductor_sig_bottom is not None:
                    conductor_properties = self.fill_missing_vals(required_keys, default_values, self.conductor_sig_bottom)
                    Rrms = [x['Rrms'] for x in conductor_properties]
                    boundary_loc = [x['boundary_loc'] for x in conductor_properties]
                    distribution = [x['distribution'] for x in conductor_properties]
                    material_properties = [{'mur': mur_dielectric, 'er': er_dielectric}] + conductor_properties
                    Zs_sig_bottom = surfz.surface_impedance(self.f, material_properties, 
                                                         Rrms=Rrms, boundary_loc=boundary_loc, distribution=distribution)
                    # actually assign the boundary
                    y_datalist = Zs_sig_bottom.real.tolist()
                    hfss.create_dataset('Zs_sig_bottom_real', x_datalist, y_datalist, is_project_dataset=False, x_unit='Hz', y_unit='ohm')
                    y_datalist = Zs_sig_bottom.imag.tolist()
                    hfss.create_dataset('Zs_sig_bottom_imag', x_datalist, y_datalist, is_project_dataset=False, x_unit='Hz', y_unit='ohm')
                    hfss.assign_impedance_to_sheet(sig.bottom_face_y, sourcename='SIG-BOT', 
                                                   resistance='pwl(Zs_sig_bottom_real, Freq)', reactance='pwl(Zs_sig_bottom_imag, Freq)', is_infground=False)
                    
                if self.conductor_sig_side is not None:
                    conductor_properties = self.fill_missing_vals(required_keys, default_values, self.conductor_sig_side)
                    Rrms = [x['Rrms'] for x in conductor_properties]
                    boundary_loc = [x['boundary_loc'] for x in conductor_properties]
                    distribution = [x['distribution'] for x in conductor_properties]
                    material_properties = [{'mur': mur_air, 'er': er_air}] + conductor_properties
                    Zs_sig_side = surfz.surface_impedance(self.f, material_properties, 
                                                         Rrms=Rrms, boundary_loc=boundary_loc, distribution=distribution)
                    # actually assign the boundary
                    y_datalist = Zs_sig_side.real.tolist()
                    hfss.create_dataset('Zs_sig_side_real', x_datalist, y_datalist, is_project_dataset=False, x_unit='Hz', y_unit='ohm')
                    y_datalist = Zs_sig_side.imag.tolist()
                    hfss.create_dataset('Zs_sig_side_imag', x_datalist, y_datalist, is_project_dataset=False, x_unit='Hz', y_unit='ohm')
                    hfss.assign_impedance_to_sheet(sig.top_face_x, sourcename='SIG-SID', 
                                                   resistance='pwl(Zs_sig_side_real, Freq)', reactance='pwl(Zs_sig_side_imag, Freq)', is_infground=False)
                    
                    
            # draw the port and assign wave-port
            port = hfss.modeler.create_polyline(points=[['0', '0', '0'], ['sub_width/2', '0', '0']],
                                                name='PORT').sweep_along_vector(sweep_vector=['0', 'port_height', '0'])
            
            P1 = hfss.wave_port(port, integration_line=[[0, self.h, 0], [0, 0, 0]],
                            modes=1, renormalize=False, deembed=0, name='P1')
            P1['Modes/Mode1/CharImp'] = 'Zpi'   # define characteristic impedance type
            
            # create air region 
            air = hfss.modeler.create_region(pad_value=[0,0,0,0,0,0], pad_type='Percentage Offset',
                                            name='AIRBOX')
            # assign symmetry boundary (H-plane symmetry)
            hfss.assign_symmetry(assignment=[air.bottom_face_x], 
                                 name='Symmetry', is_perfect_e=False)
            hfss.set_impedance_multiplier(0.5)   # to compensate for the H-symmetry
            # assign radiation boundary
            hfss.assign_radiation_boundary_to_objects(assignment=[air.top_face_x, air.top_face_y, air.top_face_z], 
                                                      name='Radiation')
            
            # fit all in screen (probably only useful if you looking at the GUI)
            hfss.modeler.fit_all()
            
            # define simulation setup for establishing the mesh
            setup1 = hfss.create_setup('setup1')
            setup1.enable_adaptive_setup_single(freq=solution_freq*1e-9, # this is stupid to must give in GHz
                                                max_passes=max_passes)
            setup1.props['BasisOrder'] = 3
            setup1.props['SolveType'] = 'PortsOnly'
            setup1.props['PortAccuracy'] = port_accuracy
            setup1.analyze()   # runs to get optimal mesh based on given parameters 
    
            # second setup settings that contain the frequency sweep
            # the mesh is linked to first setup
            setup2 = hfss.create_setup('setup2')
            setup2.props['Frequency'] = f'{solution_freq*1e-9:.6f}GHz'
            setup2.props['BasisOrder'] = 3
            setup2.props['SolveType'] = 'PortsOnly'
            setup2.props['DoLambdaRefine'] = False
            setup2.add_mesh_link(design = hfss.design_name, adapt_port=False)
            setup2.props['MeshLink']['Soln'] = f'{setup1.name} : PortOnly'
    
            # define frequency sweep setup
            sweep = setup2.create_single_point_sweep(unit='Hz', freq=freq,
                                                     save_single_field=False, save_fields=False, name='Single')
    
            setup2.analyze()  # run the simulation
    
            sol = hfss.post.get_solution_data(expressions=['Zo(P1)', 'Gamma(P1)'], setup_sweep_name=f'{setup2.name} : {sweep.name}',
                                        domain='Sweep', primary_sweep_variable='Freq')
            self.gamma = np.array([complex(x,y) for x,y in zip(sol.data_real('Gamma(P1)'), sol.data_imag('Gamma(P1)'))])
            self.Z0 = np.array([complex(x,y) for x,y in zip(sol.data_real('Zo(P1)'), sol.data_imag('Zo(P1)'))])
            
            if closs_aedt_at_finish:
                hfss.close_project()
                hfss.delete_project()
                hfss.close_desktop()
                return True
            else:
                return hfss
        except:
            if closs_aedt_at_finish:
                hfss.close_project()
                hfss.delete_project()
                hfss.close_desktop()
                return False
            else:
                return hfss

if __name__ == '__main__':
    pass
    
# EOF
