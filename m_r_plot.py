import numpy as np
from matplotlib import pyplot as plt

### Reading in our ALMA data ###

mass        = np.loadtxt('masses.txt')
radius      = np.loadtxt('radiipc.txt')
mass_5sig   = np.loadtxt('masses_5sig.txt')
radius_5sig = np.loadtxt('radiipc_5sig.txt')

### Reading in other data ###

msma     = np.array([51,78,251,2143,236,135,285,119,455,71,109,239,225,445,325])
rsma     = np.array([0.12,0.12,0.19,0.26,0.12,0.15,0.22,0.13,0.26,0.09,0.12,0.16,0.17,0.16,0.15])
louvet   = np.loadtxt('Louvet_cores.txt')
louvetr  = louvet[:,0]
louvetm  = louvet[:,1]
peretto  = np.loadtxt('Per_cores.txt')
perettor = peretto[:,0]
perettom = peretto[:,1]
lu = np.loadtxt('walker.txt',usecols=(8,10))
lur = lu[:,0]
lum = lu[:,1]
lurpc = lur*(4.848e-6)
walkerr = np.array([4847,1300,1300,1475,1275,1543,1711,1682,1179,945,1687,1236,1161,945,1139,2613,2206,1568,1302])
walkerm = np.array([64.2,16.2,18.0,2.9,0.9,1.3,3.5,2.0,0.8,0.6,4.4,2.2,2.2,1.4,2.3,9.6,3.1,1.8,1.0])
walkerrpc = walkerr*(4.848e-6)
massa    = 72
rada     = 0.04

### Setting up constant lines for plot ###

R        = np.linspace(1e-3, 3.5e1, int(1e5))
alpha    = 1
Pturb_k  = 1e9
mu       = 2.3
mh       = 1.67e-27
kb       = 1.38e-23
T        = 75
sm       = 1.989e30
pc       = 3.086e16
Pturb    = (Pturb_k)*(1e6)*kb
M        = (((4/3)*np.pi*alpha*Pturb*mu*mh)/(kb*T)) * (R*pc)**3

Pturb_k_disc = 1e5
T_disc       = 20
Pturb_disc   = (Pturb_k_disc)*(1e6)*kb
M_disc       = (((4/3)*np.pi*alpha*Pturb_disc*mu*mh)/(kb*T_disc)) * (R*pc)**3

Pturb_k_int = 1e7
T_int = 47
Pturb_int = (Pturb_k_int)*(1e6)*kb
M_int       = (((4/3)*np.pi*alpha*Pturb_int*mu*mh)/(kb*T)) * (R*pc)**3

Nh      = 2.8*(1.67e-27)
M_CD       = np.linspace(1, 1e5, int(1e5))
M_CD_si    = M_CD*sm
R_si    = R*pc

N_1e25  = 1e25
N_1e23  = 1e23

M_N_1e25 = Nh*((1e4)*np.pi*R_si*R_si*N_1e25)/sm
M_N_1e23 = Nh*((1e4)*np.pi*R_si*R_si*N_1e23)/sm


SD = 870
Mthresh = SD*(R**1.33)

### Creating plot ###

plt.figure(figsize=(24,16))
plt.xscale('log')
plt.yscale('log')
plt.xlim(0.002,0.6)
plt.ylim(0.1,4000)
plt.xlabel('Radius (pc)', fontsize=48)
plt.ylabel('Mass (M$_\odot$)', fontsize=48)
plt.xticks(fontsize=40)
plt.yticks(fontsize=40)
plt.plot(rsma, msma, marker='s', color='r', markersize=10, linewidth=0)
plt.plot(rada, massa, marker='s', mfc='w', mec='r', markersize=10, linewidth=0)
plt.plot(louvetr, louvetm, marker='^', color='b', markersize=10, linewidth=0)
plt.plot(perettor, perettom, marker='.', color='k', markersize=10, linewidth=0)
plt.plot(lurpc, lum, marker='+', color='g', markersize=10, linewidth=0)
plt.plot(walkerrpc, walkerm, marker='H', color='orange', markersize=5, linewidth=0)
plt.plot(radius, mass, marker='D', color='c', markersize=5, linewidth=0)
plt.plot(radius_5sig, mass_5sig, marker='D', color='deeppink', markersize=8, linewidth=0)
plt.plot(0.16, 239, marker='*', color='k', markersize=10, linewidth=0)
plt.plot([0.004848,0.004848], [0.01,10e4], color='r', ls='--', alpha=0.5, lw=2)
plt.plot(R, M/sm, color='k', ls='--', alpha=0.5, lw=2)
plt.plot(R, M_disc/sm, color='k', ls='--', alpha=0.5, lw=2)
plt.plot(R, M_int/sm, color='k', ls='--', alpha=0.5, lw=2)
plt.plot(R, M_N_1e25, color='k', ls='-.', alpha=0.5, lw=2)
plt.plot(R, M_N_1e23, color='k', ls='-.', alpha=0.5, lw=2)
plt.plot(R, Mthresh, color='k', ls='None', alpha=0.5, lw=2)
plt.fill_between(R, Mthresh, color='gainsboro', alpha=.5)
plt.legend(('Walker+ 18', 'Rathborne+ 14', 'Louvet+ 14', 'Peretto+ 13', 'Lu+ 20', 'Walker+ 21', 'This work (3$\sigma$)', 'This work (5$\sigma$)'), loc='lower right', fontsize='22', edgecolor='k')
plt.annotate(r'$\rho_{crit}$ (Disc)' '\nP$_e/k$ ~ 10$^5$ K cm$^{-3}$', (0.11, 0.2), fontsize=30, ha='center')
plt.annotate('P$_e/k$ ~ 10$^7$ K cm$^{-3}$', (0.00975, 0.15), fontsize=30, ha='center')
plt.annotate(r'$\rho_{crit}$ (CMZ)' '\nP$_e/k$ ~ 10$^9$ K cm$^{-3}$', (0.125, 800), fontsize=30, ha='center')
plt.annotate('10$^{25}$ cm$^{-2}$', (0.04, 2000), fontsize=30, ha='center')
plt.annotate('10$^{23}$ cm$^{-2}$', (0.4, 450), fontsize=30, ha='center')
plt.show()
plt.savefig('m_r_plot.png')