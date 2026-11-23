# Source for notebooks/13_wind_loads.ipynb — rebuild with: python build.py

# %% [markdown]
# # Wind Loads
# **ASCE 7-22 Chapters 26, 27 (MWFRS directional procedure), and 30 (C&C walls) · IBC 2024**
#
# Rectangular building with a flat, gable, or hip roof. Computes velocity pressures using the
# ASCE 7-22 Kz constants, external pressure coefficients from Fig. 27.3-1, net design pressures
# for both internal pressure signs, base shear and overturning in both principal directions
# (including the 27.1.5 minimum), and C&C wall pressures for Zones 4 and 5.
# Pressures are strength-level (load case W).
#
# **How to use:** Run the **Engine** cell once, edit the form fields below, and run **Runtime → Run all** (`Ctrl+F9` / `⌘F9`). Get V for your Risk Category from the ASCE 7 Hazard Tool.
#
# > Design aid only. A licensed engineer must verify all inputs and results.

# %%
#@title 1. Site and wind parameters
V = 115.0  #@param {type:"number"}
risk_category = "II"  #@param ["I", "II", "III", "IV"]
exposure = "C"  #@param ["B", "C", "D"]
Kzt = 1.0  #@param {type:"number"}
ground_elevation_ft = 0.0  #@param {type:"number"}
use_Ke = True  #@param {type:"boolean"}
Kd = 0.85  #@param {type:"number"}
G = 0.85  #@param {type:"number"}
enclosure = "Enclosed"  #@param ["Enclosed", "Partially enclosed", "Partially open", "Open"]
#@markdown V in mph (3-second gust, Risk Category specific map). Kd = 0.85 for buildings. G = 0.85 for rigid buildings.

# %%
#@title 2. Building geometry (ft, degrees)
roof_type = "Gable"  #@param ["Flat", "Gable", "Hip"]
roof_angle = 20.0  #@param {type:"number"}
eave_height = 20.0  #@param {type:"number"}
B_along_ridge = 100.0  #@param {type:"number"}
L_across_ridge = 60.0  #@param {type:"number"}
heights = ""  #@param {type:"string"}
cc_effective_area = 20.0  #@param {type:"number"}
#@markdown *heights*: optional list of heights for the qz table, separated by `;` (blank = automatic).
#@markdown *cc_effective_area*: effective wind area of the wall component in ft².

# %%
#@title 3. Calculate
from structcalc.calcs import wind_loads

out = wind_loads.run(
    V=V, exposure=exposure, Kzt=Kzt, ground_elev=ground_elevation_ft, use_Ke=use_Ke,
    Kd=Kd, G=G, enclosure=enclosure, roof_type=roof_type, theta=roof_angle,
    h_eave=eave_height, B_ridge=B_along_ridge, L_across=L_across_ridge,
    heights=heights, cc_area=cc_effective_area, risk_category=risk_category,
)
out["report"].show()

# %%
#@title 4. Wind pressures, C&C zones, and velocity pressure profile
fig = wind_loads.plot(out)

# %%
#@title 5. Save report as HTML (optional)
save_report = False  #@param {type:"boolean"}
if save_report:
    path = out["report"].save_html("wind_loads.html")
    print("saved", path)
    try:
        from google.colab import files
        files.download(path)
    except ImportError:
        pass
