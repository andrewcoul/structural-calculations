# Source for notebooks/14_seismic_loads.ipynb — rebuild with: python build.py

# %% [markdown]
# # Seismic Loads (ELF)
# **ASCE 7-22 Chapters 11 and 12, Equivalent Lateral Force procedure (Method 1) · IBC 2024**
#
# Computes Seismic Design Category, approximate fundamental period, seismic response coefficient Cs
# with all code limits (Eq. 12.8-2 through 12.8-6), base shear, vertical force distribution,
# story shears, and overturning moments. Enter ASCE 7-22 design values SDS, SD1, S1, and TL
# directly from the ASCE 7 Hazard Tool or USGS Seismic Design Geodatabase (no Fa/Fv factors).
#
# **How to use:** Run the **Engine** cell once, edit the form fields below, and run **Runtime → Run all** (`Ctrl+F9` / `⌘F9`).
#
# > Design aid only. A licensed engineer must verify all inputs and results.

# %%
#@title 1. Site design values (from the ASCE 7 Hazard Tool)
SDS = 1.0  #@param {type:"number"}
SD1 = 0.6  #@param {type:"number"}
S1 = 0.6  #@param {type:"number"}
TL = 8.0  #@param {type:"number"}
site_class = "D"  #@param ["A", "B", "BC", "C", "CD", "D", "DE", "E", "Default"]
risk_category = "II"  #@param ["I", "II", "III", "IV"]

# %%
#@title 2. Seismic force-resisting system (Table 12.2-1)
R = 6.5  #@param {type:"number"}
Omega0 = 3.0  #@param {type:"number"}
Cd = 4.0  #@param {type:"number"}
rho = 1.0  #@param {type:"number"}
period_system = "All other structural systems"  #@param ["Steel moment-resisting frame", "Concrete moment-resisting frame", "Steel EBF / BRBF", "All other structural systems"]
T_computed = 0.0  #@param {type:"number"}
use_12813 = False  #@param {type:"boolean"}
#@markdown Example: R = 6.5, Omega0 = 3, Cd = 4 for light-frame wood shear walls (wood structural panels).
#@markdown *T_computed*: period from a structural model in seconds (0 = use Ta). *use_12813*: cap SDS at max(1.0, 0.7 SDS) for Cs; use only when all 12.8.1.3 conditions are met.

# %%
#@title 3. Levels (bottom to top, separated by ;)
story_heights = "12; 10; 10"  #@param {type:"string"}
level_weights = "220; 210; 150"  #@param {type:"string"}
#@markdown Story heights in ft (first story first). Level weights in kip: effective seismic weight assigned to each floor/roof level (12.7.2).

# %%
#@title 4. Calculate
from structcalc.calcs import seismic_loads

out = seismic_loads.run(
    SDS=SDS, SD1=SD1, S1=S1, TL=TL, risk_category=risk_category, R=R, Omega0=Omega0, Cd=Cd,
    rho=rho, system=period_system, T_computed=T_computed, story_heights=story_heights,
    story_weights=level_weights, use_12813=use_12813, site_class=site_class,
)
out["report"].show()

# %%
#@title 5. Force distribution plot
fig = seismic_loads.plot(out)

# %%
#@title 6. Save report as HTML (optional)
save_report = False  #@param {type:"boolean"}
if save_report:
    path = out["report"].save_html("seismic_loads.html")
    print("saved", path)
    try:
        from google.colab import files
        files.download(path)
    except ImportError:
        pass
