import json, numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 8, 'font.family': 'serif', 'mathtext.fontset': 'cm',
                     'axes.linewidth': 0.6, 'lines.linewidth': 1.1})
D = '../data/'
fig, ax = plt.subplots(1, 3, figsize=(7.0, 2.25))
# (a) mu->0 rates vs N
cr = json.load(open(D + 'criterion.json'))
Ns = np.array(sorted(int(k) for k in cr))
sl = np.array([np.mean(cr[str(n)]['slow']) for n in Ns]); fa = np.array([np.mean(cr[str(n)]['fast']) for n in Ns])
a = ax[0]
a.plot(Ns, 4 + 0 * Ns, 's-', color='C3', ms=3, label=r'$H$: $q$')
a.plot(Ns, Ns - 4, '^-', color='C2', ms=3, label=r'$PH$: $N-q$')
a.plot(Ns, sl, 'o-', color='#08519c', ms=3, label='next slow (secular)')
a.plot(Ns, fa, 'v-', color='#6baed6', ms=3, label='slowest fast')
a.plot(Ns, Ns / 2, ':', color='gray', lw=0.7, label=r'$N/2$')
a.set_xlabel(r'$N$'); a.set_ylabel(r'$\lim_{\mu\to0}$ rate$/\mu$')
a.set_ylim(1.5, 10.5); a.set_xticks([6, 8, 10, 12, 14]); a.legend(frameon=False, fontsize=5.6, loc='upper left')
a.text(0.86, 0.06, '(a)', transform=a.transAxes)
# (b) secular tower
s = json.load(open(D + 'secular.json'))
Nk = np.array([int(k) for k in s if int(k) >= 14]); L = np.array([s[str(n)] for n in Nk])
for k, c in zip([1, 2, 3, 4], ['C3', '#08519c', '#3182bd', '#9ecae1']):
    ax[1].plot(1 / Nk, L[:, k], 'o-', ms=3, color=c, lw=0.7, label=rf'$\lambda_{k}$')
for k in [1, 2, 3]: ax[1].axhline(k, color='gray', lw=0.5, ls=':')
ax[1].set_xlim(0, 0.075); ax[1].set_ylim(0.8, 3.3)
ax[1].set_xlabel(r'$1/N$'); ax[1].set_ylabel(r'secular rates $/q\mu$')
ax[1].legend(frameon=False, fontsize=6, loc='upper right'); ax[1].text(0.05, 0.9, '(b)', transform=ax[1].transAxes)
# (c) B concentration
NB = np.array([10, 12, 14, 16, 18, 20.]); fl = np.array([0.4693, 0.4339, 0.3735, 0.3333, 0.3018, 0.2721])
ax[2].loglog(NB, fl, 'o', color='#08519c', ms=3.5)
p = np.polyfit(np.log(NB), np.log(fl), 1); xx = np.linspace(9, 40, 20)
ax[2].loglog(xx, np.exp(np.polyval(p, np.log(xx))), '--', color='#08519c', lw=0.7, label=rf'$\propto N^{{{p[0]:.2f}}}$')
from matplotlib.ticker import NullFormatter, FixedLocator
ax[2].xaxis.set_minor_formatter(NullFormatter()); ax[2].yaxis.set_minor_formatter(NullFormatter())
ax[2].xaxis.set_major_locator(FixedLocator([10, 20, 40])); ax[2].set_xticklabels(['10', '20', '40'])
ax[2].yaxis.set_major_locator(FixedLocator([0.2, 0.3, 0.5])); ax[2].set_yticklabels(['0.2', '0.3', '0.5'])
ax[2].set_xlabel(r'$N$'); ax[2].set_ylabel(r'$\|B-\bar B\|/\bar B$')
ax[2].legend(frameon=False, fontsize=6.5); ax[2].text(0.86, 0.9, '(c)', transform=ax[2].transAxes)
for x in ax: x.tick_params(direction='in', width=0.6, length=2.5, which='both')
fig.tight_layout(pad=0.3, w_pad=0.8)
fig.savefig('../figs/fig2_v2.pdf'); fig.savefig('../figs/fig2_v2.png', dpi=200)
