"""Deterministic figures from the saved data; no curvature recomputation."""
import csv,json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import sympy as s

ORDER=('T3','Nil','Sol','S3')
LABEL={'T3':r'$T^3$','Nil':'Nil','Sol':'Sol','S3':r'$S^3$'}
COLOR={'T3':'#333333','Nil':'#1679ad','Sol':'#bf6a24','S3':'#9b2f56'}

def number(x):return float(s.sympify(x,locals={'q':s.S.One,'log':s.log,'sqrt':s.sqrt,'pi':s.pi}))
def render(data_dir,output_dir):
    data=Path(data_dir);out=Path(output_dir);out.mkdir(parents=True,exist_ok=False)
    rows=json.loads((data/'classification.json').read_text(encoding='utf-8'))['geometries']
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':200})
    def save(fig,name):
        fig.savefig(out/(name+'.png'),bbox_inches='tight',facecolor='white')
        fig.savefig(out/(name+'.pdf'),bbox_inches='tight',metadata={'Creator':'Paper06 saved-data renderer','CreationDate':None,'ModDate':None})
        plt.close(fig)
    fig,ax=plt.subplots(figsize=(9.4,4.8));ax.axis('off')
    boxes=[(.5,.9,'Fixed class-A geometry and quotient\nStructure matrix $n$, reference volume'),
           (.15,.56,'Common scale + torsion\nAction, endpoints,\nreduced orbits'),
           (.5,.56,'Torsionless shape extension\nTorus cone; Nil/Sol evolution\nfrom tracefree Ricci'),
           (.85,.56,'Static frame comparison\nSame flat torus,\nrotated coframe'),
           (.15,.14,'Restricted $\\chi$ dependence\nNY, Pontryagin, internal dual,\nrelative primitive'),
           (.5,.14,'Additional structure data\n$J_2,J_3$: quadratic/cubic invariants;\nabsolute CS, quotient volume'),
           (.85,.14,'Representation dependence\nrank $n$ and closed forms\nversus $b_1$')]
    for x,y,txt in boxes:ax.text(x,y,txt,ha='center',va='center',fontsize=9,bbox=dict(boxstyle='round,pad=0.6',facecolor='#f0f4f7',edgecolor='#5b6c78'),transform=ax.transAxes)
    for x0,y0,x1,y1 in [(.4,.81,.17,.69),(.5,.81,.5,.69),(.6,.81,.83,.69),(.15,.43,.15,.27),(.2,.43,.45,.27),(.85,.43,.85,.27)]:
        ax.annotate('',xy=(x1,y1),xytext=(x0,y0),xycoords='axes fraction',arrowprops=dict(arrowstyle='->',color='#5b6c78'))
    save(fig,'fig01_structure')
    fig,ax=plt.subplots(figsize=(6.8,4.2))
    xx=np.linspace(-.4,1.12,100);ax.plot(xx,-9*xx,color='#777777',lw=1,label=r'$C=-9\chi$ (static family)')
    for name,offset in zip(ORDER,[(15,18),(-26,-18),(-8,20),(-26,-20)]):
        x=number(rows[name]['chi']);y=number(rows[name]['C_static']);ax.scatter(x,y,s=60,color=COLOR[name],zorder=3)
        ax.annotate(LABEL[name],(x,y),xytext=offset,textcoords='offset points',color=COLOR[name])
    ax.set(xlabel=r'$\chi$',ylabel=r'$C$',xlim=(-.45,1.15),ylim=(-11,5));ax.grid(alpha=.2);ax.legend(frameon=False);fig.tight_layout();save(fig,'fig04_static_bridge')
    fig,axs=plt.subplots(1,2,figsize=(9.2,3.7))
    loc=np.arange(4)
    for j in range(3):axs[0].bar(loc+(j-1)*.23,[number(rows[n]['n'][j]) for n in ORDER],width=.23,label=fr'$n_{j+1}$')
    axs[0].axhline(0,color='#777777',lw=.7);axs[0].set_xticks(loc,[LABEL[n] for n in ORDER]);axs[0].set_ylabel('Structure eigenvalues');axs[0].legend(ncol=3,frameon=False,fontsize=9)
    norm=[np.sqrt(sum(number(v)**2 for row in rows[n]['tracefree'] for v in row)) for n in ORDER]
    axs[1].bar(loc,norm,color=[COLOR[n] for n in ORDER],width=.55);axs[1].set_xticks(loc,[LABEL[n] for n in ORDER]);axs[1].set_ylabel(r'$q^2\|{}^{(3)}\mathrm{Ric}_{\rm TF}\|$')
    axs[0].set_title('(a) Class A in every row',fontsize=11);axs[1].set_title('(b) Omitted coframe equation',fontsize=11)
    fig.tight_layout(w_pad=3);save(fig,'fig03_structure_residual')
    with (data/'orbits.csv').open(encoding='utf-8') as f:orbits=list(csv.DictReader(f))
    fig,axs=plt.subplots(1,3,figsize=(10,3.5),sharey=True)
    axs[0].fill_between([0,4],0,3.6,color='#eeeeee',hatch='///',edgecolor='#cccccc');axs[0].text(2,1.6,'No real\nconstraint surface',ha='center',va='center')
    for n in ('T3','Nil','Sol'):
        for sheet in ([0] if n=='T3' else [-1,1]):
            dat=[r for r in orbits if r['geometry']==n and int(r['sheet'])==sheet];ax=axs[1 if n=='T3' else 2]
            ax.plot([float(r['tau']) for r in dat],[float(r['q']) for r in dat],color=COLOR[n],ls='-' if sheet>=0 else '--',label=LABEL[n] if sheet>=0 else None)
    for ax,title in zip(axs,[r'$\chi>0$: $S^3$',r'$\chi=0$: $T^3$',r'$\chi<0$: Nil, Sol']):
        ax.set(xlabel=r'Proper time $\tau$',title=title,xlim=(0,4),ylim=(0,3.6));ax.grid(alpha=.15)
    axs[0].set_ylabel(r'Scale $q$');axs[2].legend(frameon=False);fig.tight_layout();save(fig,'fig02_reduced_orbits')
    with (data/'pulse.csv').open(encoding='utf-8') as f:pulse=list(csv.DictReader(f))
    fig,axs=plt.subplots(2,1,figsize=(7.8,5.4),sharex=True)
    for n in ORDER:
        rr=[r for r in pulse if r['geometry']==n];tt=[float(r['t']) for r in rr]
        axs[0].plot(tt,[float(r['Pontryagin']) for r in rr],label=LABEL[n],color=COLOR[n])
        axs[1].plot(tt,[float(r['relative_slice']) for r in rr],color=COLOR[n])
    axs[0].set_ylabel(r'Local $\mathcal{P}$');axs[1].set_ylabel(r'Relative primitive $\phi$');axs[1].set_xlabel(r'$t=\tau$')
    axs[0].legend(ncol=4,frameon=False);axs[0].set_title(r'Off-shell pulse: $q=N=1$, $V=0$, $\eta=(2/5)e^{-t^2}$',fontsize=11)
    for ax in axs:ax.axhline(0,color='#999999',lw=.6);ax.grid(alpha=.15)
    fig.tight_layout();save(fig,'fig05_dynamic_endpoints')
    cone=json.loads((data/'torus_deformations.json').read_text(encoding='utf-8'))['mixed_shape']['figure_slice']
    h=np.linspace(*cone['h_range'],cone['h_samples']);angle=np.linspace(0,2*np.pi,cone['angle_samples'])
    hh,aa=np.meshgrid(h,angle);rr=np.sqrt(cone['normalization'])*np.abs(hh)
    fig=plt.figure(figsize=(7.5,5.4));ax=fig.add_subplot(111,projection='3d')
    ax.plot_surface(rr*np.cos(aa),rr*np.sin(aa),hh,color='#4c95b4',alpha=.55,linewidth=0,rstride=3,cstride=2)
    ax.plot([0,0],[0,0],cone['h_range'],color='#b24b3e',ls='--',lw=2,label='Pure scale (only origin allowed)')
    ax.scatter([0],[0],[0],color='#222222',s=25)
    ax.set(xlabel=r'$\Sigma_1$',ylabel=r'$\Sigma_2$',zlabel=r'$h$',title=r'Two-shear slice: $\Sigma_1^2+\Sigma_2^2=6h^2$')
    ax.view_init(elev=22,azim=-56);ax.legend(loc='upper left',fontsize=9)
    fig.text(.5,.025,r'$\Sigma=\Sigma_1E_1+\Sigma_2E_2$, $E_1=\mathrm{diag}(1,-1,0)/\sqrt{2}$, $E_2=\mathrm{diag}(1,1,-2)/\sqrt{6}$',ha='center',fontsize=10)
    fig.subplots_adjust(bottom=.12,top=.88,left=.0,right=.9);save(fig,'fig06_torus_cone')
