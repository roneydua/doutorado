# %%
#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@File    :   pure_reflection_and_pure_transmission_analysis.ipynb
@Time    :   2023/05/08 15:45:23
@Author  :   Roney D. Silva
@Contact :   roneyddasilva@gmail.com
'''
%reload_ext autoreload
%autoreload 2
%matplotlib inline

import locale

import matplotlib.pyplot as plt
import numpy as np
import sympy as sp
from IPython.core.interactiveshell import InteractiveShell
from ipywidgets import fixed, interact, interactive
from modeling.math_model_accel import AccelModel
from matplotlib import ticker
from IPython.display import display, Markdown

from bragg import Bragg
from common_functions import *
from pure_reflection_and_pure_transmission_analysis_Functions import *

plt.style.use("./common_functions/roney3.mplstyle")
locale.setlocale(locale.LC_ALL, "pt_BR.UTF-8")
my_colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
FIG_L = 6.29
FIG_A = (90.0) / 25.4
plt.rcParams["figure.dpi"] = 288
plt.rcParams["figure.figsize"] = (FIG_L, FIG_A)
# Create a FBG on lambda with 1550nm

accel = AccelModel()
max_deformation = accel.seismic_mass * 9.89 * 0.25 / accel.k / accel.fiber_length
deformation_vector = np.linspace(-max_deformation,
                                 max_deformation,
                                 5,
                                 endpoint=True)
gravity_span_vector = 4.0 * accel.k * accel.fiber_length * deformation_vector / accel.seismic_mass
update_pot_vs_gravity = True

# %% [markdown]
# # FBG com laser
# - FBG 1 com pico em 1549nm
# - FBG 2 com pico em 1549.25nm
# - Laser com potência de 1mw com pico em 1549.3nm
#   - _**O espectro do laser foi normalizado para facilitar a visualização.**_

# %%
bragg = Bragg(fbg_size= 2e-3,
              delta_n=4e-4,
              delta_span_wavelength=500,
              diff_of_peak=5,
              wavelength_peak=1549)
laser = 1e-3*calc_laser(w=bragg.wavelength_span,
                    center_w=1549.3 * 1e-9,
                    std=0.05 * 1e-9)
def simple_plot():
    fig, ax = plt.subplots(1, 2, num=1, figsize=(FIG_L, FIG_A))
    ax[0].plot(bragg.wavelength_span_nm, bragg.r0)
    total_power = np.trapz(x=bragg.wavelength_span, y=laser)
    ax[1].xaxis.set_major_locator(ticker.MultipleLocator(.5))

    ax[1].plot(bragg.wavelength_span_nm,
               laser,
               label="Potência total de " + '{:2.1f}'.format(total_power * 1e3) +
               r"\si{\milli\watt}")
    ax[1].legend()
    # plt.plot(bragg.wavelength_span_nm, laser )
    ax[0].set_ylabel("Refletividade")
    ax[1].set_ylabel(r"\si{\watt\per\meter}")
    ax[0].set_xlim([1548, 1550])
    ax[1].set_xlim([1548,1550])
    fig.supxlabel("$\\lambda$, \\si{\\nm}")


simple_plot()

# %% [markdown]
# ## Transmissão (Laser)

# %%
Trans = Transmition(bragg, laser)
interact(Trans.transmission,
         deformation=(deformation_vector.min(), deformation_vector.max(),
                      deformation_vector[1] - deformation_vector[0]))

if update_pot_vs_gravity:
    pot_transmission = calc_pot_transmission(bragg, deformation_vector, laser)
    # pot_vs_gravity(gravity_span_vector, pot_transmission, 'Transmissão',
    #             "../../images/transmission_vs_strain.pdf")


# %% [markdown]
# ## Reflexão (Laser)

# %%
Ref = Reflection(bragg, laser)
interact(Ref.reflection,
         deformation=(deformation_vector.min(), deformation_vector.max(),
                      deformation_vector[1] - deformation_vector[0]))

if update_pot_vs_gravity:
    pot_reflection = calc_pot_reflection(bragg, deformation_vector, laser)
    # pot_vs_gravity(gravity_span_vector, pot_reflection, 'Reflexão',
    #             "../../images/reflection_vs_strain.pdf")


# %% [markdown]
# ## Reflexão da reflexão (Laser)

# %%
bragg2 = Bragg(2e-3,
               delta_n=4e-4,
               delta_span_wavelength=500,
               diff_of_peak=5,
               wavelength_peak=1549.75)
RefRef = ReflectionReflection(bragg1=bragg, bragg2=bragg2, laser=laser)

interact(RefRef.reflection_reflection,
         deformation=(deformation_vector.min(), deformation_vector.max(),
                      deformation_vector[1] - deformation_vector[0]))

if update_pot_vs_gravity:
    pot_reflection_reflection = calc_pot_reflection_reflection(
        bragg, bragg2, deformation_vector, laser)

    pot_vs_gravity(gravity_span_vector, pot_reflection_reflection, "Reflexão da Reflexão", "../../images/reflection_reflection_vs_strain.pdf")

# %% [markdown]
# ## Transmissão da transmissão (Laser)

# %%
TransTrans = TransmissionTransmission(bragg1=bragg, bragg2=bragg2, laser=laser)


interact(TransTrans.transmission_transmission,
         deformation=(deformation_vector.min(), deformation_vector.max(),
                      deformation_vector[1] - deformation_vector[0]))
if update_pot_vs_gravity:
    pot_transmission_transmission = calc_pot_transmission_transmission(
        bragg, bragg2, deformation_vector, laser)
    # pot_vs_gravity(gravity_span_vector, pot_transmission_transmission,
    #             "Transmissão da transmissão",
    #             "../../images/transmission_transmission_vs_strain.pdf")


# %% [markdown]
# ## Reflexão da transmissão (Laser)

# %%
RefTrans = ReflectionTransmission(bragg1=bragg, bragg2=bragg2, laser=laser)

interact(RefTrans.reflection_transmission,
         deformation=(deformation_vector.min(), deformation_vector.max(),
                      deformation_vector[1] - deformation_vector[0]))

if update_pot_vs_gravity:
    pot_reflection_transmission = calc_pot_reflection_transmission(
        bragg, bragg2, deformation_vector, laser)

    pot_vs_gravity(gravity_span_vector, pot_reflection_transmission,
                "Reflexão da Transmissão",
                "../../images/reflection_transmission_vs_strain.pdf")

# %% [markdown]
# ## Transmissão da reflexão (Laser)

# %%
TransRef = TransmissionReflection(bragg1=bragg, bragg2=bragg2, laser=laser)

interact(TransRef.transmission_reflection,
         deformation=(deformation_vector.min(), deformation_vector.max(),
                      deformation_vector[1] - deformation_vector[0]))

if update_pot_vs_gravity:
    pot_transmission_reflection = calc_pot_transmission_reflection(
        bragg, bragg2, deformation_vector, laser)
    # pot_vs_gravity(gravity_span_vector, pot_transmission_reflection,
    #                 "Transmissão da reflexão",
    #                 "../../images/transmission_reflection_vs_strain.pdf")


# %% [markdown]
# # FBG com banda larga
# - A fonte de banda larga possui uma potência de 10mw com largura de 15nm.
# - As FBGs são as mesmas utilizada na simulação do Laser como fonte de luz.
# - 

# %%
banda_larga = 1e-3 * 10. / 15. / 1e-9 * np.ones(laser.size)
# banda_larga[250:] *= 0

def simple_plot():
    fig, ax = plt.subplots(1, 2, num=1, figsize=(FIG_L, FIG_A))
    ax[0].plot(bragg.wavelength_span_nm, bragg.r0)
    total_power = np.trapz(x=bragg.wavelength_span, y=banda_larga)
    # ax[1].xaxis.set_major_locator(ticker.MultipleLocator(.5))

    ax[1].plot(bragg.wavelength_span_nm,
               banda_larga,
               label="Potência total de " +
               '{:2.1f}'.format(total_power * 1e3) + r"\si{\milli\watt}")
    ax[1].legend()
    # plt.plot(bragg.wavelength_span_nm, laser )
    ax[0].set_ylabel("Refletividade")
    ax[1].set_ylabel(r"\si{\milli\watt\per\meter}")
    ax[0].set_xlim([1548, 1550])
    # ax[1].set_xlim([1548, 1550])
    fig.supxlabel("$\\lambda$, \\si{\\nm}")


simple_plot()

# %% [markdown]
# ## Reflexão da reflexão (Banda Larga)

# %%
RefRefBandaLarga = ReflectionReflection(bragg1=bragg, bragg2=bragg2, laser=banda_larga)


interact(RefRefBandaLarga.reflection_reflection,
         deformation=(deformation_vector.min(), deformation_vector.max(),
                      deformation_vector[1] - deformation_vector[0]))
if update_pot_vs_gravity:
    pot_reflection_reflection_banda_larga = calc_pot_reflection_reflection(
        bragg, bragg2, deformation_vector, banda_larga)

    # pot_vs_gravity(
    #     gravity_span_vector, pot_reflection_reflection_banda_larga,
    #     "Reflexão da Reflexão",
    #     "../../images/reflection_reflection_vs_strain_banda_larga.pdf")


# %% [markdown]
# ## Transmissão da transmissão (Banda larga)

# %%
TransTransBandaLarga = TransmissionTransmission(bragg1=bragg, bragg2=bragg2, laser=banda_larga)

interact(TransTransBandaLarga.transmission_transmission,
         deformation=(deformation_vector.min(), deformation_vector.max(),
                      deformation_vector[1] - deformation_vector[0]))

if update_pot_vs_gravity:
    pot_transmission_transmission_banda_larga = calc_pot_transmission_transmission(
        bragg, bragg2, deformation_vector, banda_larga)

    # pot_vs_gravity(
    #     gravity_span_vector, pot_transmission_transmission_banda_larga,
    #     "Transmissão da Transmissão",
    #     "../../images/transmission_transmission_vs_strain_banda_larga.pdf")


# %% [markdown]
# ## Reflexão da transmissão (Banda larga)

# %%
RefTransBandaLarga = ReflectionTransmission(bragg1=bragg, bragg2=bragg2, laser=banda_larga)

interact(RefTransBandaLarga.reflection_transmission,
         deformation=(deformation_vector.min(), deformation_vector.max(),
                      deformation_vector[1] - deformation_vector[0]))

if update_pot_vs_gravity:
    pot_reflection_transmission_banda_larga = calc_pot_reflection_transmission(
        bragg, bragg2, deformation_vector, banda_larga)

    # pot_vs_gravity(
    #     gravity_span_vector, pot_reflection_transmission_banda_larga,
    #     "Reflexão da Transmissão",
    #     "../../images/reflection_transmission_vs_strain_banda_larga.pdf")


# %% [markdown]
# ## Transmissão da Reflexão (Banda larga)

# %%
TransRefBandalarga = TransmissionReflection(bragg1=bragg, bragg2=bragg2, laser=banda_larga)

interact(TransRefBandalarga.transmission_reflection,
         deformation=(deformation_vector.min(), deformation_vector.max(),
                      deformation_vector[1] - deformation_vector[0]))

if update_pot_vs_gravity:
    pot_transmission_reflection_banda_larga = calc_pot_transmission_reflection(
        bragg, bragg2, deformation_vector, banda_larga)
    # pot_vs_gravity(
    #     gravity_span_vector, pot_transmission_reflection_banda_larga,
    #     "Transmissão da reflexão",
    #     "../../images/transmission_reflection_vs_strain_banda_larga.pdf")


# %%
def all_pot_vs_gravity_Laser():
    fig, ax = plt.subplots(3, 2, num=1,sharex=True, figsize=(FIG_L, 0.8*FIG_A))
    fig.supxlabel(r"Aceleração [g]")
    fig.supylabel("\\si{\\micro\\watt}")
    def plot(_ax, y,label):
        poly_coef = np.polyfit(x=gravity_span_vector/9.89, y=y * 1e6, deg=1)
        _ax.set_title(r'$p(g)[\si{\micro\watt}]=$' + '{:2.4f}'.format(poly_coef[0]) + r"$\cdot g+$"
                    '{:2.4f}'.format(poly_coef[1]))
        _ax.plot(gravity_span_vector/9.89,
                y*1e6,
                label=label)
        _ax.legend()
    # Plot reflection vs strain
    plot(ax[0,0],pot_transmission, "Transmissão")
    plot(ax[0,1],pot_reflection,'Reflexão')
    plot(ax[1,0],pot_reflection_reflection, 'Reflexão da reflexão')
    plot(ax[1,1],pot_transmission_transmission, 'Transmissão da transmissão')
    plot(ax[2,0],pot_reflection_transmission, 'Reflexão da transmissão')
    plot(ax[2,1],pot_transmission_reflection, 'Transmissão da reflexão')

    plt.savefig("../../images/all_interrogation_vs_strain_Laser.pdf", format="pdf")
    plt.close(fig=1)





def all_pot_vs_gravity_banda_larga():
    fig, ax = plt.subplots(2,
                           2,
                           num=1,
                           sharex=True,
                           figsize=(FIG_L, 0.8*FIG_A))
    fig.supxlabel(r"Aceleração [g]")
    fig.supylabel("\\si{\\micro\\watt}")

    def plot(_ax, y, label):
        poly_coef = np.polyfit(x=gravity_span_vector / 9.89, y=y * 1e6, deg=1)
        _ax.set_title(r'$p(g)[\si{\micro\watt}]=$' +
                      '{:2.4f}'.format(poly_coef[0]) + r"$\cdot g+$"
                      '{:2.4f}'.format(poly_coef[1]))
        _ax.plot(gravity_span_vector / 9.89, y * 1e6, label=label)
        _ax.legend()

    # Plot reflection vs strain
    plot(ax[0, 0], pot_reflection_reflection_banda_larga, 'Reflexão da reflexão')
    plot(ax[0, 1], pot_transmission_transmission_banda_larga, 'Transmissão da transmissão')
    plot(ax[1, 0], pot_reflection_transmission_banda_larga, 'Reflexão da transmissão')
    plot(ax[1, 1], pot_transmission_reflection_banda_larga, 'Transmissão da reflexão')

    plt.savefig("../../images/all_interrogation_vs_strain_banda_larga.pdf",
                format="pdf")
    plt.close(fig=1)


all_pot_vs_gravity_Laser()
all_pot_vs_gravity_banda_larga()

# %%



