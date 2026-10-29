import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
C=dict(blue='#2a78d6',orange='#eb6834',aqua='#1baf7a',yellow='#eda100',ink='#0b0b0b',ink2='#52514e',grid='#d9d8d4',fill='#cde2fb')
def setup():
    plt.rcParams.update({'font.size':8,'font.family':'serif','mathtext.fontset':'cm','axes.linewidth':0.6,
        'axes.edgecolor':C['ink2'],'xtick.color':C['ink2'],'ytick.color':C['ink2'],'axes.labelcolor':C['ink'],
        'xtick.major.width':0.6,'ytick.major.width':0.6,'xtick.major.size':2.5,'ytick.major.size':2.5,
        'legend.frameon':False,'legend.fontsize':6.8,'lines.linewidth':1.4,'lines.markersize':4,'savefig.dpi':300})
def tag(ax,s,x=0.03,y=0.93): ax.text(x,y,s,transform=ax.transAxes,fontsize=8.5,color=C['ink'],fontweight='bold')
