import json, numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 8, 'font.family': 'serif', 'mathtext.fontset': 'cm',
                     'axes.linewidth': 0.6, 'lines.linewidth': 1.1})
D = '../data/'
fig, ax = plt.subplots(1, 3, figsize=(7.0, 2.3))
# (a) ladder eigenvalue curves
cur = json.load(open(D + 'ladder_curves3.json'))
cols = {'0.1': '#08519c', '0.25': '#3182bd', '0.4': '#9ecae1'}
for mu, c in cols.items():
    d = cur[mu]; E = np.array(d['E']) / float(mu); lam = np.array(d['lam'])
    ax[0].plot(E, lam[:, 0], '-', color=c, label=rf'$\mu={float(mu):g}J$')
    ax[0].plot(E, lam[:, 1], ':', color=c, lw=0.8)
    ax[0].axvline(d['thr'] / float(mu), color=c, lw=0.5, ls='--')
ax[0].axhline(1, color='k', lw=0.5); ax[0].axvline(4, color='C3', lw=0.8)
ax[0].set_xlim(0, 16); ax[0].set_ylim(0, 3)
ax[0].set_xlabel(r'$E/\mu$'); ax[0].set_ylabel(r'ladder eigenvalue $k(E)$')
ax[0].legend(frameon=False, fontsize=6, loc='upper left'); ax[0].text(0.87, 0.92, '(a)', transform=ax[0].transAxes)
# (b) physical bound states vs threshold
b = ax[1]
g = np.array(json.load(open(D + 'gamma1_q4.json'))); mu_sd, G1 = g[:, 0], -g[:, 1]
sel = mu_sd > 0.015
b.plot(mu_sd[sel], 2 * G1[sel] / mu_sd[sel], 'k-', label=r'threshold $2\Gamma_1/\mu$')
b.axhline(4, color='C3', lw=0.8, label=r'energy mode $q\mu$')
bs = []
for f in ['ladder_branch_a.json', 'ladder_branch_b.json', 'ladder_branch_c.json']:
    for r in json.load(open(D + f)):
        for j, e in r['roots']:
            if j == 0: bs.append([r['mu'], e])
bs = np.array(bs)
b.plot(bs[:, 0], bs[:, 1], 'o', color='C3', ms=3, mfc='none', label='ladder bound state')
b.axvline(0.3228, color='gray', lw=0.5, ls='--')
b.set_xlim(0, 0.5); b.set_ylim(2, 16)
b.set_xlabel(r'$\mu/J$'); b.set_ylabel(r'rate$/\mu$')
b.legend(frameon=False, fontsize=6, loc='upper right'); b.text(0.04, 0.06, '(b)', transform=b.transAxes)
# (c) mu_c(N) fast-mode identification
c = ax[2]
Ns, m, e = [], [], []
for N in [10, 12, 14, 16, 18]:
    try: v = [r['mu_c'] for r in json.load(open(D + f'crossing2_N{N}.json'))]
    except FileNotFoundError: continue
    Ns.append(N); m.append(np.mean(v)); e.append(np.std(v) / np.sqrt(len(v)) if len(v) > 1 else 0.006)
Ns = np.array(Ns, float); m = np.array(m)
c.errorbar(1 / Ns, m, e, fmt='o', ms=3, color='#08519c', capsize=1.5, label='ED (fast-mode crossing)')
xx = np.linspace(0, 0.1, 100)
for pw, ls in [(1, '--'), (0.5, ':')]:
    pf = np.polyfit(Ns[1:] ** (-pw), m[1:], 1)
    c.plot(xx, np.polyval(pf, xx ** pw), ls, color='#08519c', lw=0.7, label=rf'linear in $N^{{-{pw:g}}}$: {pf[-1]:.2f}')
c.legend(frameon=False, fontsize=6, loc='lower left')
c.plot([0], [0.3228], 's', color='k', ms=4)
c.text(0.005, 0.334, r'ladder: $q\mu_c=2\Gamma_1$', fontsize=6.5)
c.set_xlim(0, 0.105); c.set_ylim(0, 0.36)
c.set_xlabel(r'$1/N$'); c.set_ylabel(r'$\mu_c(N)/J$')
c.text(0.88, 0.92, '(c)', transform=c.transAxes)
for a in ax: a.tick_params(direction='in', width=0.6, length=2.5)
fig.tight_layout(pad=0.3, w_pad=0.8)
fig.savefig('../figs/fig1_v2.pdf'); fig.savefig('../figs/fig1_v2.png', dpi=200)
print(list(zip(Ns, np.round(m, 4))))
