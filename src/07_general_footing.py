# Source for notebooks/07_general_footing.ipynb — rebuild with: python build.py

# %% [markdown]
# # General Footing
# **ACI 318-19 · ASCE 7-22 (ASD for soil, LRFD for concrete) · IBC 2024**
#
# Isolated rectangular spread footing under a concentric rectangular column or pedestal.
# Checks soil bearing with biaxial eccentricity (exact partial-contact solution, no tension),
# overturning, sliding, and uplift, plus concrete flexure at both column faces, minimum steel,
# bar spacing, short-direction band steel, one-way shear, two-way (punching) shear with
# unbalanced moment, bearing, and development length.
#
# **How to use:** Run the **Engine** cell once, edit the form fields below, and run **Runtime → Run all** (`Ctrl+F9` / `⌘F9`).
#
# > Design aid only. A licensed engineer must verify all inputs and results.

# %%
#@title 1. Footing geometry and materials
Lx = 8.0  #@param {type:"number"}
Ly = 8.0  #@param {type:"number"}
thickness = 24.0  #@param {type:"number"}
col_x = 18.0  #@param {type:"number"}
col_y = 18.0  #@param {type:"number"}
pedestal_height = 0.0  #@param {type:"number"}
soil_cover = 1.0  #@param {type:"number"}
fc = 3000.0  #@param {type:"number"}
fy = 60000.0  #@param {type:"number"}
gamma_concrete = 150.0  #@param {type:"number"}
gamma_soil = 110.0  #@param {type:"number"}
#@markdown Lx, Ly, soil cover in **ft**. Thickness, column sizes, pedestal height in **in**. f'c, fy in **psi**; unit weights in **pcf**.

# %%
#@title 2. Reinforcement (bottom mat)
bar_x = 6  #@param [4, 5, 6, 7, 8, 9, 10, 11] {type:"raw"}
spacing_x = 9.0  #@param {type:"number"}
bar_y = 6  #@param [4, 5, 6, 7, 8, 9, 10, 11] {type:"raw"}
spacing_y = 9.0  #@param {type:"number"}
cover = 3.0  #@param {type:"number"}
#@markdown Bars along x are in the bottom layer. Spacing and cover in **in**.

# %%
#@title 3. Soil
allowable_pressure = 3.0  #@param {type:"number"}
allowable_is_net = False  #@param {type:"boolean"}
increase_for_wind_seismic = 1.0  #@param {type:"number"}
friction_coefficient = 0.35  #@param {type:"number"}
passive_pressure = 0.0  #@param {type:"number"}
fs_overturning_required = 1.0  #@param {type:"number"}
fs_sliding_required = 1.0  #@param {type:"number"}
#@markdown Allowable pressure in **ksf**. Passive pressure as an equivalent fluid in **pcf** (0 = ignore).
#@markdown FS required applies to ASD combinations (0.6D+0.6W, 0.6D+0.7E): 1.0 is usual; enter 1.5 for service-load practice.

# %%
#@title 4. Service loads at the top of the pedestal
#@markdown Each case: `P, Mx, My, Vx, Vy` in **kip** and **kip-ft** (P + down). +My and +Vx move the resultant toward +x; +Mx and +Vy toward +y. Leave blank for none.
load_D = "80, 0, 0, 0, 0"  #@param {type:"string"}
load_L = "60, 0, 0, 0, 0"  #@param {type:"string"}
load_Lr = ""  #@param {type:"string"}
load_S = ""  #@param {type:"string"}
load_W = "0, 0, 40, 5, 0"  #@param {type:"string"}
load_E = ""  #@param {type:"string"}

# %%
#@title 5. Calculate
from structcalc.calcs import general_footing

col_loads = {c: general_footing.parse_load_string(s) for c, s in
             {"D": load_D, "L": load_L, "Lr": load_Lr, "S": load_S, "W": load_W, "E": load_E}.items()
             if s.strip()}
out = general_footing.run(
    Lx=Lx, Ly=Ly, h=thickness, col_x=col_x, col_y=col_y, ped_height=pedestal_height,
    soil_cover=soil_cover, gamma_c=gamma_concrete, gamma_s=gamma_soil,
    qa=allowable_pressure, qa_increase_WE=increase_for_wind_seismic, qa_is_net=allowable_is_net,
    mu=friction_coefficient, passive_pcf=passive_pressure,
    fs_overturn=fs_overturning_required, fs_sliding=fs_sliding_required,
    fc=fc, fy=fy, bar_x=bar_x, s_x=spacing_x, bar_y=bar_y, s_y=spacing_y, cover=cover,
    col_loads=col_loads,
)
out["report"].show()

# %%
#@title 6. Soil pressure plan, critical shear planes, and section elevation
fig = general_footing.plot(out)

# %%
#@title 7. Save report as HTML (optional)
save_report = False  #@param {type:"boolean"}
if save_report:
    path = out["report"].save_html("general_footing.html")
    print("saved", path)
    try:
        from google.colab import files
        files.download(path)
    except ImportError:
        pass
