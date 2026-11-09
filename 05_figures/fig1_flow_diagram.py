import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
fig,ax=plt.subplots(figsize=(9,13)); ax.set_xlim(0,18); ax.set_ylim(1.5,30); ax.axis('off')
C='#1f4e79'; MX,MW,EX,EW,H=1,7,10,7.5,1.2
def box(x,y,w,h,t,bold=False,fc='#eef3f8'):
    ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle="round,pad=0.02,rounding_size=0.15",fc=fc,ec=C,lw=1.2))
    ax.text(x+w/2,y+h/2,t,ha='center',va='center',fontsize=9,fontweight='bold' if bold else 'normal')
def arrow(p,q): ax.annotate('',xy=q,xytext=p,arrowprops=dict(arrowstyle='->',color=C,lw=1.2))
def flow(title,mains,excl,ytop):
    ax.text(9,ytop+H+0.7,title,ha='center',fontsize=11.5,fontweight='bold')
    y=ytop; box(MX,y,MW,H,mains[0],True)
    for ex,nxt in zip(excl,mains[1:]):
        yn=y-2.3; cx=MX+MW/2
        arrow((cx,y),(cx,yn+H))
        box(EX,yn+H+0.3,EW,1.0,ex,fc='#fbf1e6'); arrow((cx,yn+H+0.8),(EX,yn+H+0.8))
        box(MX,yn,MW,H,nxt,True); y=yn
    return y
y=flow("MIMIC-IV (primary cohort)",
 ["1,045 ICU-stay records, belonging to\n795 admissions with a repair code","791 admissions","780 admissions","764 admissions","740 admissions analysed"],
 ["Excluded: 4 coded with both open and\nendovascular repair","Excluded: 11 coded only with the\nlocation-unspecified code 39.78","Excluded: 16 with end-stage renal disease\nor chronic dialysis","Excluded: 24 with no creatinine in the\nfirst 24 hours"],27)
cx=MX+MW/2
box(MX-0.3,y-2.3,3.6,1.2,"Open repair\n261 admissions"); box(MX+3.7,y-2.3,3.6,1.2,"Endovascular repair\n479 admissions")
arrow((cx,y),(MX+1.5,y-1.1)); arrow((cx,y),(MX+5.5,y-1.1))
y2=flow("eICU-CRD (exploratory; repair type not recorded)",
 ["446 eligible stays","395 stays","384 stays","382 stays analysed"],
 ["Excluded: 51 dissecting aneurysms\nwithout a location","Excluded: 11 with no creatinine","Excluded: 2 with end-stage renal disease"],11.8)
box(MX-0.3,y2-2.3,3.6,1.2,"Abdominal\n232 stays"); box(MX+3.7,y2-2.3,3.6,1.2,"Thoracic\n150 stays")
arrow((cx,y2),(MX+1.5,y2-1.1)); arrow((cx,y2),(MX+5.5,y2-1.1))
plt.savefig('Figure1_flow.png',dpi=170,bbox_inches='tight'); print(y,y2)
