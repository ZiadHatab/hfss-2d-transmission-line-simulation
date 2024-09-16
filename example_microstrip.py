"""
Author: @Ziad (https://github.com/ZiadHatab)

Example to simulate microstrip line using hfss_ms.py, which runs HFSS to solve a 2D full-wave microstrip line.

The material parameters I'm using are not really correct, I just retrofit them to the measurement.
Especially Nickel, I'm sure I have it wrong. Also the substrate uses a constant rel. permittivity (which is for sure wrong).
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os

from hfss_ms import HFSSMS   # my code, need to be in same folder as this script

def LLmodel(f, muri, murf, f0, gamma):
    """
    Landau-Lifshits model for complex relative permeability.
    Y. Shlepnev and S. McMorrow, "Nickel characterization for interconnect analysis," 
    2011 IEEE International Symposium on Electromagnetic Compatibility, Long Beach, CA, USA, 2011, 
    pp. 524-529, doi: https://doi.org/10.1109/ISEMC.2011.6038368.
    """
    return murf + (muri - murf)*(f0**2 + 1j*f*gamma)/( f0**2 + 2j*f*gamma - f**2)

if __name__ == '__main__':
    # microstrip measurement of effective permittivity and per unit length loss (these are just for comparison)
    script_folder = os.path.dirname(os.path.realpath(__file__))

    df = pd.read_csv(script_folder + '\\measurements\\Copper_Dkeff_losses.csv' )
    fmeas = df['Freq (GHz)'].values
    ereff_copper = df['Effective Dk (1)'].values
    ereff_copper_std = df['Effective Dk std (1)'].values
    loss_copper  = df['Losses (dB/mm)'].values
    loss_copper_std = df['Losses std (dB/mm)'].values

    df = pd.read_csv(script_folder + '\\measurements\\ENIG_Dkeff_losses.csv' )
    fmeas = df['Freq (GHz)'].values
    ereff_enig = df['Effective Dk (1)'].values
    ereff_enig_std = df['Effective Dk std (1)'].values
    loss_enig  = df['Losses (dB/mm)'].values
    loss_enig_std = df['Losses std (dB/mm)'].values

    # microstrip dimensions    
    h = 0.127e-3
    t = 0.02e-3
    w = 0.284e-3
    wgnd = 5e-3

    # frequency
    f = np.logspace(0, np.log10(150), 30, base=10)*1e9
    
    # ISOAL Astra® MT77
    # https://www.isola-group.com/pcb-laminates-prepreg/astra-mt77-laminate-and-prepreg/
    # @10GHz (this is not accurate, a proper model is needed)
    er = 2.91
    etand = 0.0017
    
    # Bare copper simulation
    # top gnd and bottom sig are rough
    conductor_gnd_top = [ {'sigma': 58e6, 'mur': 0.999991, 'Rrms': 1e-6, 'boundary_loc': 0, 'distribution': 'norm'} ]
    conductor_sig_bottom = [ {'sigma': 58e6, 'mur': 0.999991, 'Rrms': 1e-6, 'boundary_loc': 0, 'distribution': 'norm'} ]
    
    # top and side of sig are platted (here not, hence only copper listed)
    conductor_sig_top = [ {'sigma': 58e6, 'mur': 0.999991, 'Rrms': 50e-9, 'boundary_loc': 0, 'distribution': 'norm'} ]
    conductor_sig_side = conductor_sig_top
    
    ms_copper = HFSSMS(f, w, h, t, wgnd, er=er, etand=etand, conductor_gnd_top=conductor_gnd_top,
                conductor_sig_top=conductor_sig_top, conductor_sig_bottom=conductor_sig_bottom,
                conductor_sig_side=conductor_sig_side)
    hfss = ms_copper.run_simulation(closs_aedt_at_finish=True, headless=True, port_accuracy=0.001)
    
    # ENIG simulation
    # top gnd and bottom sig are rough
    conductor_gnd_top = [ {'sigma': 58e6, 'mur': 0.999991, 'Rrms': 1e-6, 'boundary_loc': 0, 'distribution': 'norm'} ]
    conductor_sig_bottom = [ {'sigma': 58e6, 'mur': 0.999991, 'Rrms': 1e-6, 'boundary_loc': 0, 'distribution': 'norm'} ]

    # top and side of sig are platted
    conductor_sig_top = [{'sigma': 41e6, 'mur': 0.99996, 'Rrms': 50e-9, 'boundary_loc': 0, 'distribution': 'norm'},  # gold
                         {'sigma': 14.5e6, 'mur': LLmodel(f, 20, 1, 120e9, 0.2*120e9), 'Rrms': 50e-9, 'boundary_loc': 0.05e-6, 'distribution': 'norm'},  # nickel
                         {'sigma': 58e6, 'mur': 0.999991, 'Rrms': 50e-9, 'boundary_loc': 5.2e-6+0.05e-6, 'distribution': 'norm'}  # copper
                         ]
    conductor_sig_side = conductor_sig_top  # same as top
    
    ms_enig = HFSSMS(f, w, h, t, wgnd, er=er, etand=etand, conductor_gnd_top=conductor_gnd_top,
                conductor_sig_top=conductor_sig_top, conductor_sig_bottom=conductor_sig_bottom,
                conductor_sig_side=conductor_sig_side)
    hfss = ms_enig.run_simulation(closs_aedt_at_finish=True, headless=True, port_accuracy=0.001)
    
    ereff_copper_sim = ms_copper.gamma2ereff(ms_copper.gamma, ms_copper.f)
    ereff_enig_sim = ms_enig.gamma2ereff(ms_enig.gamma, ms_enig.f)

    # PLOTS
    plt.figure()
    plt.plot(f/1e9, ereff_copper_sim.real, lw=2, label='Simulated Copper')
    plt.plot(f/1e9, ereff_enig_sim.real, lw=2, label='Simulated ENIG')
    plt.plot(fmeas, ereff_copper, lw=2, label='Measured Copper')
    plt.fill_between(fmeas, ereff_copper+2*ereff_copper_std, ereff_copper-2*ereff_copper_std, color='C2', alpha=0.3)
    plt.plot(fmeas, ereff_enig, lw=2, label='Measured ENIG')
    plt.fill_between(fmeas, ereff_enig+2*ereff_enig_std, ereff_enig-2*ereff_enig_std, color='C3', alpha=0.3)
    plt.xlabel('Frequency (GHz)')
    plt.ylabel('Effective Dk')
    plt.title('Real part of effective relative permittivity')
    plt.xlim([0, 150])
    plt.legend()
    
    plt.figure()
    plt.plot(f/1e9, ms_copper.gamma2dbmm(ms_copper.gamma), lw=2, label='Simulated Copper')
    plt.plot(f/1e9, ms_enig.gamma2dbmm(ms_enig.gamma), lw=2, label='Simulated ENIG')
    plt.plot(fmeas, loss_copper, lw=2, label='Measured Copper')
    plt.fill_between(fmeas, loss_copper+2*loss_copper_std, loss_copper-2*loss_copper_std, color='C2', alpha=0.3)
    plt.plot(fmeas, loss_enig, lw=2, label='Measured ENIG')
    plt.fill_between(fmeas, loss_enig+2*loss_enig_std, loss_enig-2*loss_enig_std, color='C3', alpha=0.3)
    plt.xlabel('Frequency (GHz)')
    plt.ylabel('Loss (dB/mm)')
    plt.title('Loss per unit length')
    plt.xlim([0, 150])
    plt.legend()
    
    plt.figure()
    plt.plot(f/1e9, ms_copper.Z0.real, lw=2, label='Simulated Copper')
    plt.plot(f/1e9, ms_enig.Z0.real, lw=2, label='Simulated ENIG')
    plt.xlabel('Frequency (GHz)')
    plt.ylabel('Z0 (ohm)')
    plt.title('Real part of characteristic impedance')
    plt.xlim([0, 150])
    plt.legend()
    
    plt.figure()
    plt.plot(f/1e9, ms_copper.Z0.imag, lw=2, label='Simulated Copper')
    plt.plot(f/1e9, ms_enig.Z0.imag, lw=2, label='Simulated ENIG')
    plt.xlabel('Frequency (GHz)')
    plt.ylabel('Z0 (ohm)')
    plt.title('Imaginary part of characteristic impedance')
    plt.xlim([0, 150])
    plt.legend()


    plt.show()

# EOF