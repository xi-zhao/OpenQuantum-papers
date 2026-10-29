import json, glob, numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 8, 'font.family': 'serif', 'mathtext.fontset': 'cm',
                     'axes.linewidth': 0.6, 'lines.linewidth': 1.1})
D = '../data/'
g = np.array(json.load(open(D + 'gamma1_q4.json')))
mu_sd, G1 = g[:, 0], -g[:, 1]
fig, ax = plt.subplots(1, 3, figsize=(7.0, 2.25))
# (a) gap / mu
a = ax[0]
cols = {12: '#9ecae1', 14: '#6baed6', 16: '#3182bd', 18: '#08519c'}
for N, c in cols.items():
    d = np.array(json.load(open(D + f'grid_N{N}.json')))
    mus = np.unique(d[:, 1]); m = [d[d[:, 1] == x, 2].mean() / x for x in mus]
    a.plot(mus, m, 'o', ms=2.6, color=c, label=f'$N={N}$')
sel = mu_sd > 0.02
a.plot(mu_sd[sel], 2 * G1[sel] / mu_sd[sel], 'k-', label=r'$2\Gamma_1/\mu$ (SD, $N\!=\!\infty$)')
a.axhline(4, color='C3', lw=1.1, label=r'energy mode $q\mu$')
for k in [2, 3]:
    a.axhline(4 * k, color='C3', lw=0.6, ls=':')
a.axvline(0.3228, color='gray', lw=0.5, ls='--')
a.set_xlim(0, 0.82); a.set_ylim(1.5, 13)
a.set_xlabel(r'$\mu/J$'); a.set_ylabel(r'decay rate$/\mu$')
a.legend(frameon=False, fontsize=6, loc='upper right'); a.text(0.02, 0.92, '(a)', transform=a.transAxes)
# (b) mu_c(N)
b = ax[1]
Ns, m, e = [], [], []
for N in [10, 12, 14, 16, 18]:
    v = [r['mu_c'] for r in json.load(open(D + f'crossing_N{N}.json'))]
    Ns.append(N); m.append(np.mean(v)); e.append(np.std(v) / np.sqrt(len(v)) if len(v) > 1 else 0.005)
Ns = np.array(Ns); x = 1 / Ns
b.errorbar(x, m, e, fmt='o', ms=3, color='#08519c', capsize=1.5)
p = np.polyfit(x[1:], np.array(m)[1:], 1); xx = np.linspace(0, 0.1, 10)
b.plot(xx, np.polyval(p, xx), '--', color='#08519c', lw=0.7)
b.plot([0], [0.3228], 's', color='k', ms=4)
b.text(0.006, 0.322, r'SD ($N\!=\!\infty$): 0.323', fontsize=6.5, va='center')
b.set_xlim(0, 0.105); b.set_ylim(0, 0.37)
b.set_xlabel(r'$1/N$'); b.set_ylabel(r'$\mu_c$')
b.text(0.88, 0.92, '(b)', transform=b.transAxes)
# (c) large q
c = ax[2]
lq = json.load(open(D + 'largeq_muc.json'))
qs = np.array([int(k) for k in lq]); r = np.array([v[0] / v[1] for v in lq.values()])
c.plot(1 / qs, r, 'o', color='k', ms=3.5)
c.plot([0, 0.27], [1, 1 + 1.8 * 0.27], ':', color='gray', lw=0.7)
c.axhline(1, color='C3', lw=0.8)
c.set_xlabel(r'$1/q$'); c.set_ylabel(r'$\mu_c^{\rm SD}/\mu_c^{(q\to\infty)}$')
c.set_xlim(0, 0.27); c.set_ylim(0.9, 1.65)
c.text(0.03, 0.92, '(c)', transform=c.transAxes)
for axx in ax:
    axx.tick_params(direction='in', width=0.6, length=2.5)
fig.tight_layout(pad=0.3, w_pad=0.8)
fig.savefig('../figs/fig1.pdf'); fig.savefig('../figs/fig1.png', dpi=200)
print('mu_c fit intercept', p[-1])
