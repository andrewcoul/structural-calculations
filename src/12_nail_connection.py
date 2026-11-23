# Source for notebooks/12_nail_connection.ipynb — rebuild with: python build.py

# %% [markdown]
# # Nailed / Wood Screw Connection
# **NDS 2024 Chapter 12 (ASD or LRFD) · lateral, withdrawal, and combined loads**
#
# Single-shear nailed or wood-screwed connection with lumber or wood structural panel side members.
# Evaluates all yield modes (NDS Table 12.3.1A), withdrawal (NDS 12.2), combined lateral and
# withdrawal loads (NDS 12.4.1), toe-nail, end-grain, and diaphragm factors, penetration checks,
# and the required number of fasteners.
#
# **How to use:** Run the **Engine** cell once, edit the form fields below, and run **Runtime → Run all** (`Ctrl+F9` / `⌘F9`).
#
# > Design aid only. A licensed engineer must verify all inputs and results.

# %%
#@title 1. Fastener
fastener = "Common nail"  #@param ["Common nail", "Box nail", "Sinker nail", "Wood screw"]
size = "16d"  #@param ["6d", "8d", "10d", "12d", "16d", "20d", "30d", "40d", "50d", "60d", "#6", "#8", "#10", "#12", "#14"]
screw_length = 2.5  #@param {type:"number"}
D_override = 0.0  #@param {type:"number"}
Fyb_override = 0.0  #@param {type:"number"}
#@markdown Nail sizes 6d-60d (NDS Table L4); wood screw gauges #6-#14 (NDS Table L3) with *screw_length* in inches.
#@markdown Overrides: 0 = use the table value (D in inches, Fyb in psi).

# %%
#@title 2. Members (inches)
side_member = "Douglas Fir-Larch"  #@param ["Douglas Fir-Larch", "Douglas Fir-Larch (North)", "Douglas Fir-South", "Hem-Fir", "Hem-Fir (North)", "Southern Pine", "Spruce-Pine-Fir", "Spruce-Pine-Fir (South)", "Western Woods", "Redwood (close grain)", "Eastern Softwoods", "Plywood - Structural I / Marine", "Plywood - other grades", "OSB - all grades", "Mixed Maple", "Red Oak", "Mixed Oak", "Yellow Poplar"]
t_side = 1.5  #@param {type:"number"}
main_member = "Douglas Fir-Larch"  #@param ["Douglas Fir-Larch", "Douglas Fir-Larch (North)", "Douglas Fir-South", "Hem-Fir", "Hem-Fir (North)", "Southern Pine", "Spruce-Pine-Fir", "Spruce-Pine-Fir (South)", "Western Woods", "Redwood (close grain)", "Eastern Softwoods", "Mixed Maple", "Red Oak", "Mixed Oak", "Yellow Poplar"]
t_main = 3.5  #@param {type:"number"}
install = "Side grain"  #@param ["Side grain", "End grain", "Toe-nail"]
diaphragm = False  #@param {type:"boolean"}
#@markdown Side member = the member that the head bears on. For toe-nails, side thickness is not used.

# %%
#@title 3. Factors and loads
method = "ASD"  #@param ["ASD", "LRFD"]
CD = 1.0  #@param {type:"number"}
CM = 1.0  #@param {type:"number"}
Ct = 1.0  #@param {type:"number"}
lam = 1.0  #@param {type:"number"}
n_fasteners = 4  #@param {type:"integer"}
lateral_load = 400.0  #@param {type:"number"}
withdrawal_load = 0.0  #@param {type:"number"}
#@markdown Loads in lb on the whole group (ASD-level for ASD, factored for LRFD).

# %%
#@title 4. Calculate
from structcalc.calcs import nail_connection

out = nail_connection.run(
    fastener=fastener, size=size, screw_length=screw_length,
    D_override=D_override, Fyb_override=Fyb_override,
    side_member=side_member, main_member=main_member, t_side=t_side, t_main=t_main,
    install=install, diaphragm=diaphragm, CD=CD, CM=CM, Ct=Ct, method=method, lam=lam,
    n_fasteners=n_fasteners, lateral_load=lateral_load, withdrawal_load=withdrawal_load,
)
out["report"].show()

# %%
#@title 5. Fastener cross-section and lateral capacity vs. penetration
fig = nail_connection.plot(out)

# %%
#@title 6. Save report as HTML (optional)
save_report = False  #@param {type:"boolean"}
if save_report:
    path = out["report"].save_html("nail_connection.html")
    print("saved", path)
    try:
        from google.colab import files
        files.download(path)
    except ImportError:
        pass
