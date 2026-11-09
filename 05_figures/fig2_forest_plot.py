# Figure 2 (5 Nov 2026): forest plot rebuilt from the Table 2 rows (hierarchy baseline primary). Values are from the S7 table (S7_rows.csv).
import csv, re, matplotlib, matplotlib.ticker, matplotlib.transforms as mt
matplotlib.use("Agg")
import matplotlib.pyplot as plt
V={}
for r in csv.DictReader(open('S7_rows.csv')):
    m=re.fullmatch(r'(\d+\.\d+) \((\d+\.\d+)–(\d+\.\d+)\)',r['estimate'])
    if m: V[r['ID']]=tuple(map(float,m.groups()))
spec=[("Primary (hierarchy baseline)",None),("Adjusted","H1"),("Unadjusted","H1u"),("Overlap-weighted","H8"),
("Baseline definition",None),("First-24-hour baseline alone, all admissions","R1"),("Pre-admission subset: first-24-hour baseline","P0"),("Pre-admission subset: pre-admission baseline","P1"),
("Population",None),("First admission per patient","H4"),("Excluding ruptured repairs","H5"),
("Model",None),("Plus creatinine draws per day","H7"),("Plus any vasopressor","H6"),
("Outcome definition",None),("Stage 2 or 3","H2"),("Stage 3","H3"),("Absolute rise 0.3 mg/dL (no baseline)","O3"),("Dialysis (no baseline)","D1"),("In-hospital death","M1"),
("Aneurysm location",None),("Abdominal","H9a"),("Thoracic","H9b"),("Thoracoabdominal","H9c")]
fig=plt.figure(figsize=(9,8)); ax=fig.add_axes([0.42,0.08,0.34,0.90]); y=len(spec); bt=mt.blended_transform_factory(ax.transAxes,ax.transData)
for name,k in spec:
    y-=1
    if k is None: ax.text(-1.2,y,name,transform=bt,fontweight='bold',fontsize=9,va='center'); continue
    e,lo,hi=V[k]; ax.plot([lo,hi],[y,y],color='#1f4e79',lw=1.6); ax.plot(e,y,'s',color='#1f4e79',ms=5)
    ax.text(-1.2,y,"   "+name,transform=bt,fontsize=8.5,va='center'); ax.text(1.04,y,f"{e:.2f} ({lo:.2f} to {hi:.2f})",transform=bt,fontsize=8,va='center')
ax.axvline(1,color='grey',lw=0.8,ls='--'); ax.set_xscale('log'); ax.set_xlim(0.4,6); ax.set_xticks([0.5,1,2,4]); ax.set_xticklabels(['0.5','1','2','4'])
ax.xaxis.set_minor_locator(matplotlib.ticker.NullLocator()); ax.set_yticks([]); ax.set_ylim(-1,len(spec))
for s in ('left','right','top'): ax.spines[s].set_visible(False)
ax.set_xlabel("Risk ratio, open vs endovascular repair (log scale)"); ax.text(1.04,len(spec)-0.2,"RR (95% CI)",transform=bt,fontsize=8.5,fontweight='bold',va='center')
plt.savefig('Figure2_forest_v2.png',dpi=200); print(len(spec))
