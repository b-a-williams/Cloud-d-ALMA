from astropy import units as u
from astropy import constants as const
import math as m
import numpy as np

### Setting up constants and parameters using values from Walker et al. (2018) ###

mu = 2.33
k_B_SI = const.k_B.to(u.kg*u.m**2/u.s**2/u.K)
Pk = 10**8
T = 20 * u.K
T2 = 40 * u.K
mh = 1.67262192369e-27 * u.kg

pcrit = (4/3)*m.pi*((Pk*mu*mh)/T)

sound_speed_m_s = ((k_B_SI * T) / (mu*mh))**0.5
sound_speed_km_s = sound_speed_m_s.to(u.km/u.s)

d6_vd = (5 * u.km/u.s) / ((8*np.log(2))**0.5)

mass_upper_msol = 239 * u.solMass
mass_upper_kg = mass_upper_msol.to(u.kg)
mass_lower_msol = 69 * u.solMass
mass_lower_kg = mass_lower_msol.to(u.kg)
rad_pc = 0.16 * u.pc
rad_m = rad_pc.to(u.m)
vol_m3 = (4/3) * m.pi * (rad_m)**3

density_upper_kgm3 = mass_upper_kg/vol_m3
density_lower_kgm3 = mass_lower_kg/vol_m3

density_lower_m3 = density_lower_kgm3/mh
density_lower_cm3 = density_lower_m3.to(u.cm**-3)
density_upper_m3 = density_upper_kgm3/mh
density_upper_cm3 = density_upper_m3.to(u.cm**-3)

grav_cons = 6.674e-11 * (u.m**3/u.kg/u.s**2)

### Calculating thermal and turbulent lengths and masses ###

jeans_upper_pc = (0.4*(sound_speed_km_s/(0.2 * u.km/u.s))*(density_lower_cm3/(10**3 / u.cm**3))**(-0.5))* u.pc
jeans_lower_pc = (0.4*(sound_speed_km_s/(0.2 * u.km/u.s))*(density_upper_cm3/(10**3 / u.cm**3))**(-0.5))* u.pc

jeans_lower_AU = jeans_lower_pc.to(u.AU)
jeans_upper_AU = jeans_upper_pc.to(u.AU)

turb_lower_pc = (0.4*(d6_vd/(0.2 * u.km/u.s))*(density_upper_cm3/(10**3 / u.cm**3))**(-0.5))* u.pc
turb_upper_pc = (0.4*(d6_vd/(0.2 * u.km/u.s))*(density_lower_cm3/(10**3 / u.cm**3))**(-0.5))* u.pc

turb_lower_AU = turb_lower_pc.to(u.AU)
turb_upper_AU = turb_upper_pc.to(u.AU)

mj_upper = (2*((sound_speed_km_s/(0.2 * u.km/u.s))**3)*(density_lower_cm3/(10**3 / u.cm**3))**(-0.5))* u.solMass
mj_lower = (2*((sound_speed_km_s/(0.2 * u.km/u.s))**3)*(density_upper_cm3/(10**3 / u.cm**3))**(-0.5))* u.solMass

mt_upper = (2*((d6_vd/(0.2 * u.km/u.s))**3)*(density_lower_cm3/(10**3 / u.cm**3))**(-0.5))* u.solMass
mt_lower = (2*((d6_vd/(0.2 * u.km/u.s))**3)*(density_upper_cm3/(10**3 / u.cm**3))**(-0.5))* u.solMass