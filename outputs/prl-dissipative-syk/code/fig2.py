import json, numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 8, 'font.family': 'serif', 'mathtext.fontset': 'cm',
                     'axes.linewidth': 0.6, 'lines.linewidth': 1.1})
D = '../data/'; mu = 0.05
fig, ax = plt.subplots(1, 3, figsize=(7.0, 2.2))
cs = {0.5: '#c6dbef', 1.0: '#6baed6', 2.0: '#2171b5', 4.0: '#08306b'}
for b0, c in cs.items():
    r = np.array(json.load(open(D + f'therm_N14_mu0.05_b{b0}.json')))
    t = r[:, 0] * 4 * mu
    ax[0].semilogy(t, r[:, 1], 'o', ms=2, color=c)
    ax[0].semilogy(t, r[:, 3] / b0, '-', color=c, lw=0.9, label=rf'$\beta_0={b0:g}$')
    ax[1].plot(t, r[:, 4], '-', color=c, label=rf'$\beta_0={b0:g}$')
    ax[1].plot(t, r[:, 5], ':', color=c, lw=0.8)
tt = np.linspace(0, 3, 50)
ax[0].semilogy(tt, np.exp(-tt), 'k--', lw=0.7)
ax[0].set_xlabel(r'$q\mu t$'); ax[0].set_ylabel(r'$E(t)/E_0$ (dots), $\beta(t)/\beta_0$ (lines)')
ax[0].legend(frameon=False, fontsize=6); ax[0].text(0.85, 0.9, '(a)', transform=ax[0].transAxes)
ax[1].set_xlabel(r'$q\mu t$'); ax[1].set_ylabel(r'trace distance')
ax[1].text(0.85, 0.9, '(b)', transform=ax[1].transAxes)
ax[1].text(1.2, 0.42, r'$\|\rho(t)-\rho_{\beta_0}\|$', fontsize=6.5)
ax[1].text(1.2, 0.13, r'$\|\rho(t)-\rho_{\beta(t)}\|$', fontsize=6.5)
s = json.load(open(D + 'secular.json'))
Ns = np.array([int(k) for k in s if int(k) >= 14]); L = np.array([s[str(n)] for n in Ns])
for k, c in zip([1, 2, 3, 4], ['C3', '#2171b5', '#6baed6', '#9ecae1']):
    ax[2].plot(1 / Ns, L[:, k], 'o-', ms=3, color=c, lw=0.7, label=rf'$\lambda_{k}$')
for k in [1, 2, 3]:
    ax[2].axhline(k, color='gray', lw=0.5, ls=':')
ax[2].set_xlim(0, 0.075); ax[2].set_ylim(0.8, 3.2)
ax[2].set_xlabel(r'$1/N$'); ax[2].set_ylabel(r'secular rates $/q\mu$')
ax[2].legend(frameon=False, fontsize=6, loc='upper right'); ax[2].text(0.05, 0.9, '(c)', transform=ax[2].transAxes)
for a in ax: a.tick_params(direction='in', width=0.6, length=2.5)
fig.tight_layout(pad=0.3, w_pad=0.8)
fig.savefig('../figs/fig2.pdf'); fig.savefig('../figs/fig2.png', dpi=200)
