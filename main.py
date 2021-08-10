from lib_opa import *
from lib_opa.calculation import allComponents  # , get_indices
import matplotlib.pyplot as plt
from lib_opa.new_refrac import get_material
from lib_opa.calculation import generate_indices
# https://whimsical.com/XSrhVZ1abs7SHSRh2UMVdg


def main():
    print('New ------------------------------------------\n')
    # Convention idler - signal --> pump
    OPA1 = Opa('ooe')
    OPA1.set_polar
    F1 = (Framework(dt=0.1e-15,
                    T=2e-12,
                    dr=25e-6,
                    R=1125e-6))
    Pump = (Pulse(lambda_c=0.515e-6,
                  duration=100e-15,
                  beta2=00_000e-30,
                  delay=0,
                  energy=125e-6,
                  radius=600e-6,
                  Framework=F1))
    Signal = (Pulse(lambda_c=1.03e-6,
                    duration=125e-15,
                    beta2=0e-30,
                    delay=0,
                    energy=1e-9,
                    radius=600e-6,
                    Framework=F1))

    # Signal.pulse_duration()
    Pump.pulse_duration()
    # LBO = Crystal('LBO', theta=15, phi=0, plan = 'XZ')
    BBO = Crystal('BBO', theta=23.4)
    OPA1.add_component(F1, "OpaFramework")
    OPA1.add_component(Pump, "PumpBeam")
    OPA1.add_component(Signal, "SignalBeam")
    OPA1.add_component(BBO, "OpaCrystal")
    allComponents(OPA1)
    generate_indices(OPA1)
    Pump.pulse_duration()
    plt.plot(np.real(Pump.pulse_t))
    plt.show()
    # print(Pump.pulse_t)


if __name__ == "__main__":
    main()
