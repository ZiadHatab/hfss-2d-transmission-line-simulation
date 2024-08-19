"""
Author: @Ziad (https://github.com/ZiadHatab)

A code to run HFSS to solve a 2D full-wave microstrip line Z0 and gamma

work in progress....

"""

import pyaedt
import numpy as np

class HFSSMS:
    """
    Run AEDT HFSS in headlesss (or graphic) mode to simulate a microstrip line 
    The script gives back the cross-section parameters, Z0 and gamma, 
    as well as (when enabled) their jacobian with respect to input parameters (TODO).
    
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

    ## Geometry parameters

    w : number
        width of signal line in meters.
    h : number
        height of dielectric substrate in meters.
    t : number
        thickness of conductor in meters. Same used for GND (it doesn't really matter, as GND is at port boundary)
    wgnd : number
        width of GND plane in meters. Also, the width of the substrate!

    ## material properties... dielectric

    er : number
        relative permittivity
    etand : number
        dielectric loss tangent

    ## material properties... conductor

    mur : number
        relative permeability
    mutand : number
        magnetic loss tangent
    sigmar : number
        relative conductivity reference to copper (58e6 S/m)
        
    use_pec : boolean
            if you want use perfect electric conductor instead specified values
    """
    def __init__(self, f, w, h, t, wgnd, er=1, etand=0, mur=1, mutand=0, sigmar=1, use_pec=False):
        
        self.f  = np.atleast_1d(f)
        
        # geometry
        self.w    = w
        self.h    = h
        self.t    = t
        self.wgnd = wgnd
        # dielectric
        self.er    = er
        self.etand = etand
        # conductor
        self.mur    = mur
        self.mutand = mutand
        self.sigmar = sigmar
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


    def run_simulation(self, closs_aedt_at_finish=True, solution_freq=None, port_accuracy=0.01, max_passes=20, headless=True):
        
        # solution freq if not provided
        solution_freq = self.f[-1]*0.8 if solution_freq is None else solution_freq
        
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
            # create dielectric material
            dielectric = hfss.materials.add_material('my_dielectric')
            dielectric.conductivity = 0
            dielectric.permittivity = self.er
            dielectric.permeability = 1
            dielectric.dielectric_loss_tangent = self.etand
            dielectric.material_appearance = [64,128,128,0]
            dielectric.update()
            
            # create conductor material 
            conductor = hfss.materials.add_material('my_conductor')
            conductor.conductivity = self.sigmar*self.sigma_cu
            conductor.permittivity = 1
            conductor.permeability = self.mur
            conductor.magnetic_loss_tangent = self.mutand
            conductor.dielectric_loss_tangent = 0
            conductor.material_appearance = [255,174,26,0]
            conductor.update()
            conductor_name = 'pec' if self.use_pec else conductor.name
    
            # simulation parameters
            hfss['line_length'] = 1e-3  # this isn't important as we well simulate in 2D, it cannot be zero!
            hfss['sub_height']  = self.h
            hfss['trace_thickness'] = self.t
            hfss['sig_width'] = self.w
            hfss['sub_width'] = self.wgnd
            hfss['port_height'] = 'sub_width/2' # this is just an approximate... I probably should look into a better approach
    
            # draw the substrate 
            sub = hfss.modeler.create_box(origin=['-sub_width/2', '0', '0'], 
                                        sizes=['sub_width', 'sub_height', 'line_length'],
                                        material=dielectric.name,
                                        name='SUB')
            sub.transparency = 0
    
            # draw the signal line
            sig = hfss.modeler.create_box(origin=['-sig_width/2', 'sub_height', '0'], 
                                      sizes=['sig_width', 'trace_thickness', 'line_length'],
                                      material=conductor_name,
                                      name='SIG')
            sig.transparency = 0
            
            # draw the ground plane
            gnd = hfss.modeler.create_box(origin=['-sub_width/2', '0', '0'], 
                                        sizes=['sub_width', '-trace_thickness', 'line_length'],
                                        material=conductor_name,
                                        name='GND')
            gnd.transparency = 0
            
            # draw the port and assign wave-port
            port = hfss.modeler.create_polyline(points=[['-sub_width/2', '0', '0'], ['sub_width/2', '0', '0']],
                                                name='PORT').sweep_along_vector(sweep_vector=['0', 'port_height', '0'])
            
            P1 = hfss.wave_port(port, integration_line=[[0, self.h, 0], [0, 0, 0]],
                            modes=1, renormalize=False, deembed=0, name='P1')
            P1['Modes/Mode1/CharImp'] = 'Zpi'   # define characteristic impedance type
    
            # create air region 
            air = hfss.modeler.create_region(pad_value=[0,0,0,0,0,0], pad_type='Percentage Offset',
                                            name='AIRBOX')
            hfss.assign_radiation_boundary_to_objects(air, name='Radiation')
    
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
            freq = self.f.tolist()
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
    h = 0.1e-3
    t = 0.02e-3
    w = 0.24e-3
    wgnd = 5e-3

    f = np.logspace(0.01, 2, 30, base=10)*1e9
    
    er = 3
    etand = 0.002
    
    ms = HFSSMS(f, w, h, t, wgnd, er=er, etand=etand)

    hfss = ms.run_simulation(closs_aedt_at_finish=True, headless=False, port_accuracy=0.001)

    import matplotlib.pyplot as plt
    ereff = ms.gamma2ereff(ms.gamma, ms.f)
        
    plt.figure()
    plt.plot(f*1e-9, ereff.real, lw=2)
    plt.xlabel('Frequency (GHz)')
    plt.ylabel('Effective Dk Real')
    plt.xlim([0, 140])
    
    plt.figure()
    plt.plot(f*1e-9, -ereff.imag/ereff.real, lw=2)
    plt.xlabel('Frequency (GHz)')
    plt.ylabel('Effective tand')
    plt.xlim([0, 140])
    
    plt.figure()
    plt.plot(f*1e-9, ms.Z0.real, lw=2, label='real')
    plt.xlabel('Frequency (GHz)')
    plt.ylabel('Impedance Z0 (ohm)')
    plt.legend()
    plt.xlim([0, 140])
    
    plt.figure()
    plt.plot(f*1e-9, ms.Z0.imag, lw=2, label='imag')
    plt.xlabel('Frequency (GHz)')
    plt.ylabel('Impedance Z0 (ohm)')
    plt.legend()
    plt.xlim([0, 140])
    
    
    plt.show()
    
    
# EOF
