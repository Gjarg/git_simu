from lib_opa import *
from lib_opa.new_refrac import get_material, material_info_dict
import numpy as np


def all_pulses(OPA):
    lp = OPA.SYSTEM['PumpBeam'].lambda_c
    ls = OPA.SYSTEM['SignalBeam'].lambda_c
    F1 = OPA.SYSTEM['OpaFramework']
    Idler = (Pulse(lambda_c=1 / (1 / lp - 1 / ls),
                   duration=None,
                   beta2=None,
                   delay=None,
                   energy=None,
                   radius=None,
                   Framework=F1))
    Pump_2 = (Pulse(lambda_c=lp/2,
                    duration=None,
                    beta2=None,
                    delay=None,
                    energy=None,
                    radius=None,
                    Framework=F1))
    Signal_2 = (Pulse(lambda_c=ls/2,
                      duration=None,
                      beta2=None,
                      delay=None,
                      energy=None,
                      radius=None,
                      Framework=F1))
    Idler_2 = (Pulse(lambda_c=(1 / (1 / lp - 1 / ls))/2,
                     duration=None,
                     beta2=None,
                     delay=None,
                     energy=None,
                     radius=None,
                     Framework=F1))
    return Idler, Pump_2, Signal_2, Idler_2


def allComponents(OPA):
    Idler, Pump_2, Signal_2, Idler_2 = all_pulses(OPA)
    OPA.add_component(Idler, "IdlerBeam")
    OPA.add_component(Pump_2, "PumpBeam_2")
    OPA.add_component(Signal_2, "SignalBeam_2")
    OPA.add_component(Idler_2, "IdlerBeam_2")


def generate_indices(OPA):
    crystal = get_material(OPA.SYSTEM['OpaCrystal'].name)
    if crystal.__bases__[0] == material_info_dict['BBO'].__bases__[0]:
        print('c est uniaxe')
        get_indices_uniaxe(OPA, crystal)
    elif crystal.__bases__[0] == material_info_dict['PPLN'].__bases__[0]:
        print(' c qpm')
        get_indices_qpm(OPA, crystal)
    elif crystal.__bases__[0] == material_info_dict['LBO'].__bases__[0]:
        print(' c biaxe')
        get_indices_biaxial(OPA, crystal)

# dev


def get_indices_biaxial(OPA, crystal):
    print("la fonction get_indices_biaxial est active")
    n_i, n_s, n_p = get_indice_name_biaxe(crystal, OPA)
    print('voici ni', n_i[0:15], 'voici ns', n_s[0:15], 'voici np', n_p[0:15])


def get_indice_name_biaxe(crystal, OPA):
    n = [0, 0, 0]
    print(OPA.SYSTEM['OpaCrystal'].plan)
    if OPA.SYSTEM['OpaCrystal'].plan == 'XY':
        print('je suis dans le plan XY **')
        print(OPA.SYSTEM['OpaCrystal'].phi)
        if OPA.Opa_PM[0] == 'e':
            print('je suis bien dans le cas 0 = e')
            # idler refractive index
            n[0] = crystal.n_angle(wavelength=OPA.SYSTEM['IdlerBeam'].lp,
                                   theta=OPA.SYSTEM['OpaCrystal'].theta,
                                   phi=OPA.SYSTEM['OpaCrystal'].phi,
                                   plan=OPA.SYSTEM['OpaCrystal'].plan)
        else:
            n[0] = crystal.nz(OPA.SYSTEM['IdlerBeam'].lp)

        if OPA.Opa_PM[1] == 'e':
            # signal refractive index
            n[1] = crystal.n_angle(OPA.SYSTEM['SignalBeam'].lp,
                                   theta=OPA.SYSTEM['OpaCrystal'].theta,
                                   phi=OPA.SYSTEM['OpaCrystal'].phi,
                                   plan=OPA.SYSTEM['OpaCrystal'].plan)
        else:
            n[1] = crystal.nz(OPA.SYSTEM['SignalBeam'].lp)

        if OPA.Opa_PM[2] == 'e':
            # pump refractive index
            n[2] = crystal.n_angle(OPA.SYSTEM['PumpBeam'].lp,
                                   theta=OPA.SYSTEM['OpaCrystal'].theta,
                                   phi=OPA.SYSTEM['OpaCrystal'].phi,
                                   plan=OPA.SYSTEM['OpaCrystal'].plan)
        else:
            n[2] = crystal.nz(OPA.SYSTEM['PumpBeam'].lp)
        return n

    elif OPA.SYSTEM['OpaCrystal'].plan == 'YZ':
        print('je suis dans le plan YZ**')
        print(OPA.SYSTEM['OpaCrystal'].phi)
        if OPA.Opa_PM[0] == 'e':
            print('je suis bien dans le cas 0 = e')
            # idler refractive index
            n[0] = crystal.n_angle(wavelength=OPA.SYSTEM['IdlerBeam'].lp,
                                   theta=OPA.SYSTEM['OpaCrystal'].theta,
                                   phi=OPA.SYSTEM['OpaCrystal'].phi,
                                   plan=OPA.SYSTEM['OpaCrystal'].plan)
        else:
            n[0] = crystal.nx(OPA.SYSTEM['IdlerBeam'].lp)

        if OPA.Opa_PM[1] == 'e':
            # signal refractive index
            n[1] = crystal.n_angle(OPA.SYSTEM['SignalBeam'].lp,
                                   theta=OPA.SYSTEM['OpaCrystal'].theta,
                                   phi=OPA.SYSTEM['OpaCrystal'].phi,
                                   plan=OPA.SYSTEM['OpaCrystal'].plan)
        else:
            n[1] = crystal.nx(OPA.SYSTEM['SignalBeam'].lp)

        if OPA.Opa_PM[2] == 'e':
            # pump refractive index
            n[2] = crystal.n_angle(OPA.SYSTEM['PumpBeam'].lp,
                                   theta=OPA.SYSTEM['OpaCrystal'].theta,
                                   phi=OPA.SYSTEM['OpaCrystal'].phi,
                                   plan=OPA.SYSTEM['OpaCrystal'].plan)
        else:
            n[2] = crystal.nx(OPA.SYSTEM['PumpBeam'].lp)
        return n

    elif OPA.SYSTEM['OpaCrystal'].plan == 'XZ':
        print('je suis dans le plan XZ**')
        print(OPA.SYSTEM['OpaCrystal'].phi)
        print(crystal.VZ2)
        if OPA.Opa_PM[0] == 'e':
            print('je suis bien dans le cas 0 = e')
            # idler refractive index
            n[0] = crystal.n_angle(wavelength=OPA.SYSTEM['IdlerBeam'].lp,
                                   theta=OPA.SYSTEM['OpaCrystal'].theta,
                                   phi=OPA.SYSTEM['OpaCrystal'].phi,
                                   plan=OPA.SYSTEM['OpaCrystal'].plan)
        else:
            n[0] = crystal.ny(OPA.SYSTEM['IdlerBeam'].lp)

        if OPA.Opa_PM[1] == 'e':
            # signal refractive index
            n[1] = crystal.n_angle(OPA.SYSTEM['SignalBeam'].lp,
                                   theta=OPA.SYSTEM['OpaCrystal'].theta,
                                   phi=OPA.SYSTEM['OpaCrystal'].phi,
                                   plan=OPA.SYSTEM['OpaCrystal'].plan)
        else:
            n[1] = crystal.ny(OPA.SYSTEM['SignalBeam'].lp)

        if OPA.Opa_PM[2] == 'e':
            # pump refractive index
            n[2] = crystal.n_angle(OPA.SYSTEM['PumpBeam'].lp,
                                   theta=OPA.SYSTEM['OpaCrystal'].theta,
                                   phi=OPA.SYSTEM['OpaCrystal'].phi,
                                   plan=OPA.SYSTEM['OpaCrystal'].plan)
        else:
            n[2] = crystal.ny(OPA.SYSTEM['PumpBeam'].lp)
        return n


def get_indices_uniaxe(OPA, crystal):
    print('Voici les indices de refraction pour uni')
    n_i, n_s, n_p = get_indice_name_uni(crystal, OPA)
    print(n_s)  # [0:15])
    print(n_p)  # [0:15])
    print(np.nan_to_num(n_i, nan=0))  # [0:15])
    beta_i, beta_s, beta_p = get_beta(crystal, OPA, np.nan_to_num(
        n_i, nan=0), np.nan_to_num(n_s, nan=0), np.nan_to_num(n_p, nan=0))
    beta1_i, beta1_s, beta1_p = get_beta1(OPA, beta_i, beta_s, beta_p)
    vg0i, vg0s, vg0p = get_velocity(beta1_i, beta1_s, beta1_p, OPA)
    print(np.round(vg0p * 100 / c, 2), np.round(vg0s * 100 / c, 2), np.round(vg0i * 100 / c, 2))
    length = 1.25e-3
    delay_t = delay_in_crystal(vg0p, vg0s,length)
    print(f'Delay of {np.round(delay_t,3)}fs between pump and signal after {length} m of propagation')

def delay_in_crystal(vg0p, vg0s,l):
    t = (l/vg0s - l/vg0p)*1e15
    return t

def get_velocity(beta1_i, beta1_s, beta1_p, OPA):
    vg0i = 1/beta1_i[int(OPA.SYSTEM['OpaFramework'].nt/2)]
    vg0s = 1/beta1_s[int(OPA.SYSTEM['OpaFramework'].nt/2)]
    vg0p = 1/beta1_p[int(OPA.SYSTEM['OpaFramework'].nt/2)]
    return vg0i, vg0s, vg0p


def get_beta1(OPA, beta_i, beta_s, beta_p):
    beta1_i = np.gradient(beta_i, OPA.SYSTEM['OpaFramework'].w)
    beta1_s = np.gradient(beta_s, OPA.SYSTEM['OpaFramework'].w)
    beta1_p = np.gradient(beta_p, OPA.SYSTEM['OpaFramework'].w)
    return beta1_i, beta1_s, beta1_p


def get_beta(crystal, OPA, n_i, n_s, n_p):
    beta_i = crystal.beta(
        n=n_i, omega_center=OPA.SYSTEM['IdlerBeam'].wc, omega_range=OPA.SYSTEM['OpaFramework'].w)
    beta_s = crystal.beta(
        n=n_s, omega_center=OPA.SYSTEM['SignalBeam'].wc, omega_range=OPA.SYSTEM['OpaFramework'].w)
    beta_p = crystal.beta(
        n=n_p, omega_center=OPA.SYSTEM['PumpBeam'].wc, omega_range=OPA.SYSTEM['OpaFramework'].w)
    return beta_i, beta_s, beta_p


def get_indice_name_uni(crystal, OPA):
    n = [0, 0, 0]
    if OPA.Opa_PM[0] == 'e':
        # idler refractive index
        n[0] = crystal.n_theta(OPA.SYSTEM['IdlerBeam'].lp,
                               OPA.SYSTEM['OpaCrystal'].theta)
    else:
        n[0] = crystal.nx(OPA.SYSTEM['IdlerBeam'].lp)

    if OPA.Opa_PM[1] == 'e':
        # idler refractive index
        n[1] = crystal.n_theta(OPA.SYSTEM['SignalBeam'].lp,
                               OPA.SYSTEM['OpaCrystal'].theta)
    else:
        n[1] = crystal.nx(OPA.SYSTEM['SignalBeam'].lp)

    if OPA.Opa_PM[2] == 'e':
        # idler refractive index
        n[2] = crystal.n_theta(OPA.SYSTEM['PumpBeam'].lp,
                               OPA.SYSTEM['OpaCrystal'].theta)
    else:
        n[2] = crystal.nx(OPA.SYSTEM['PumpBeam'].lp)
    return n


def get_indices_qpm(OPA, crystal):
    print('Voici les indices de refraction pour qpm')
    #crystal = get_material(OPA.SYSTEM['OpaCrystal'].name)
    n_p = crystal.n(OPA.SYSTEM['PumpBeam'].lp, T=140)
    betap = crystal.beta(
        n=n_p, omega_center=OPA.SYSTEM['PumpBeam'].wc, omega_range=OPA.SYSTEM['OpaFramework'].w)
    print(betap[0:15])


def get_group_velocity(OPA):
    crystal = get_material(OPA.SYSTEM['OpaCrystal'].name)