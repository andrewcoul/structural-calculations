# Source for notebooks/15_snow_loads.ipynb — rebuild with: python build.py

# %% [markdown]
# # Snow Loads
# **ASCE 7-22 Chapter 7 · IBC 2024**
#
# Flat and sloped roof snow loads using ASCE 7-22 strength-level ground snow loads (no snow
# importance factor), thermal factor from Table 7.3-3, slope factor from Fig. 7.4-1, minimum
# roof snow load, rain-on-snow surcharge, unbalanced gable/hip loads, leeward and windward roof
# step drifts with winter wind parameter W2, and sliding snow.
#
# **How to use:** Run the **Engine** cell once, edit the form fields below, and run **Runtime → Run all** (`Ctrl+F9` / `⌘F9`). Get pg and W2 for your Risk Category from the ASCE 7 Hazard Tool.
#
# > Design aid only. A licensed engineer must verify all inputs and results.

# %%
#@title 1. Site
pg = 30.0  #@param {type:"number"}
W2 = 0.55  #@param {type:"number"}
risk_category = "II"  #@param ["I", "II", "III", "IV"]
terrain = "C"  #@param ["B", "C", "D", "Above tree line (windswept mountains)", "Alaska, no trees within 2 mi"]
exposure = "Partially exposed"  #@param ["Fully exposed", "Partially exposed", "Sheltered"]
#@markdown pg in psf (strength level, Risk Category specific). W2 = winter wind parameter (about 0.25 to 0.65).

# %%
#@title 2. Roof
thermal = "Heated, unventilated roof (Table 7.3-3)"  #@param ["Heated, unventilated roof (Table 7.3-3)", "Kept just above freezing / cold ventilated roof", "Unheated or open air", "Freezer building", "Continuously heated greenhouse"]
roof_R_value = 30.0  #@param {type:"number"}
roof_type = "Flat"  #@param ["Flat", "Monoslope", "Gable", "Hip"]
roof_slope_deg = 0.0  #@param {type:"number"}
slippery_surface = False  #@param {type:"boolean"}
eave_to_ridge_ft = 20.0  #@param {type:"number"}
short_prismatic_rafters = False  #@param {type:"boolean"}
#@markdown *short_prismatic_rafters*: simply supported prismatic members from eave to ridge with W <= 20 ft (7.6.1 simplified unbalanced case).

# %%
#@title 3. Roof step drift (lower roof)
include_drift = True  #@param {type:"boolean"}
upper_roof_length_ft = 100.0  #@param {type:"number"}
lower_roof_length_ft = 40.0  #@param {type:"number"}
step_height_ft = 10.0  #@param {type:"number"}
#@markdown Upper roof length = leeward fetch; lower roof length = windward fetch, both measured in the drift direction.

# %%
#@title 4. Sliding snow from an upper sloped roof
include_sliding = False  #@param {type:"boolean"}
upper_roof_slope_deg = 26.57  #@param {type:"number"}
upper_eave_to_ridge_ft = 20.0  #@param {type:"number"}
upper_slippery = False  #@param {type:"boolean"}

# %%
#@title 5. Calculate
from structcalc.calcs import snow_loads

out = snow_loads.run(
    pg=pg, W2=W2, risk_category=risk_category, terrain=terrain, exposure=exposure,
    thermal=thermal, R_value=roof_R_value, roof_type=roof_type, theta=roof_slope_deg,
    slippery=slippery_surface, eave_to_ridge=eave_to_ridge_ft,
    prismatic_short=short_prismatic_rafters, drift=include_drift,
    upper_len=upper_roof_length_ft, lower_len=lower_roof_length_ft, step_h=step_height_ft,
    sliding=include_sliding, upper_theta=upper_roof_slope_deg,
    upper_eave_to_ridge=upper_eave_to_ridge_ft, upper_slippery=upper_slippery,
)
out["report"].show()

# %%
#@title 6. Drift surcharge plot
fig = snow_loads.plot(out, lower_roof_length_ft)

# %%
#@title 7. Save report as HTML (optional)
save_report = False  #@param {type:"boolean"}
if save_report:
    path = out["report"].save_html("snow_loads.html")
    print("saved", path)
    try:
        from google.colab import files
        files.download(path)
    except ImportError:
        pass
