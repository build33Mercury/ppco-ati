from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SD = ROOT / "source_data"
OUT = ROOT / "submission_figures"
OUT.mkdir(exist_ok=True)

def save(name):
    plt.tight_layout(); plt.savefig(OUT/f"{name}.png", dpi=300, bbox_inches="tight"); plt.savefig(OUT/f"{name}.pdf", bbox_inches="tight"); plt.close()

fig=plt.figure(figsize=(8.5,5.2)); ax=fig.add_axes([0,0,1,1]); ax.axis("off")
boxes=[(0.05,0.68,0.23,0.16,"Archived source\nrepresentation"),(0.38,0.68,0.23,0.16,"ATI support test\nexact / bounded /\nnon-identifiable / UNKNOWN"),(0.72,0.74,0.22,0.12,"Exact recovery\nATI0"),(0.72,0.53,0.22,0.12,"Model/prior-assisted\nprediction ATI3/ATI4"),(0.72,0.32,0.22,0.12,"Non-identifiable / abstain\nATI5/ATI6"),(0.38,0.24,0.23,0.16,"Prospective archive\ndesign for declared\ntarget family")]
for x,y,w,h,t in boxes:
    ax.add_patch(plt.Rectangle((x,y),w,h,fill=False,linewidth=1.4)); ax.text(x+w/2,y+h/2,t,ha="center",va="center",fontsize=10)
for a,b in [((0.28,0.76),(0.38,0.76)),((0.61,0.76),(0.72,0.80)),((0.61,0.74),(0.72,0.59)),((0.61,0.71),(0.72,0.38)),((0.495,0.68),(0.495,0.40))]:
    ax.annotate("",xy=b,xytext=a,arrowprops=dict(arrowstyle="->",lw=1.4))
ax.text(0.50,0.08,"Certificates apply to declared operator bytes and models;\nprediction does not upgrade information support.",ha="center",fontsize=10); save("Figure1_ATI_claim_calibration_framework")

df=pd.read_csv(SD/"Figure2_covariance_RFE.csv"); x=np.arange(len(df)); w=0.36; plt.figure(figsize=(7.5,4.5)); plt.bar(x-w/2,df["ABIDE"],w,label="ABIDE TRAIN"); plt.bar(x+w/2,df["HCP"],w,label="HCP external"); plt.xticks(x,df["Method"],rotation=20,ha="right"); plt.ylabel("Mean covariance relative Frobenius error"); plt.ylim(0,1.08); plt.legend(); save("Figure2_covariance_prediction_error")

df=pd.read_csv(SD/"Figure3_overlap_minus_PPCO.csv"); y=np.arange(len(df)); est=df["Estimate"].to_numpy(); lo=df["CI_low"].to_numpy(); hi=df["CI_high"].to_numpy(); plt.figure(figsize=(7.2,3.8)); plt.errorbar(est,y,xerr=np.vstack([est-lo,hi-est]),fmt="o",capsize=5); plt.axvline(0,linewidth=1); plt.yticks(y,df["Cohort"]); plt.xlabel("Overlap minus PPCO covariance RFE\n(negative values favor overlap)"); plt.gca().invert_yaxis(); save("Figure3_PPCO_benchmark_contrasts")

df=pd.read_csv(SD/"Figure4_operator_perturbation_self_controls.csv"); plt.figure(figsize=(8.2,4.6)); plt.barh(df["Perturbation family"],df["Fraction_exact"]); plt.xlim(0,1.05); plt.xlabel("Fraction of perturbed self controls remaining exact"); plt.gca().invert_yaxis(); save("Figure4_operator_perturbation_sensitivity")

df=pd.read_csv(SD/"Figure5_scaling_measurements.csv"); plt.figure(figsize=(7.2,4.4));
for case,g in df[df["Latent_n"]==902629].groupby("Case"):
    plt.plot(g["Source_rows"],g["Median_seconds"],marker="o",label=case.replace("_"," "))
plt.xlabel("Source operator rows"); plt.ylabel("Median analysis time (s)"); plt.legend(); save("Figure5_sparse_Gram_scaling")
