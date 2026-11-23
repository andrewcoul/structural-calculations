# Source for notebooks/05_concrete_beam.ipynb — rebuild with: python build.py

# %% [markdown]
# # Concrete Beam
# **ACI 318-19 · ASCE 7-22 LRFD load combinations · IBC 2024**
#
# Rectangular or T-beam with single-span, cantilever, overhang, fixed-end, or continuous supports.
# Checks positive and negative flexure by strain compatibility (doubly reinforced, φ from net
# tensile strain), minimum steel, strain limit εt ≥ 0.004, shear with stirrups (including size
# effect when Av < Av,min), stirrup spacing, bar fit, crack-control spacing, and immediate plus
# long-term deflection using ACI 318-19 effective moment of inertia Ie.
#
# **How to use:** Run the **Engine** cell once, edit the form fields below, and run **Runtime → Run all** (`Ctrl+F9` / `⌘F9`).
#
# > Design aid only. A licensed engineer must verify all inputs and results.

# %%
#@title 1. Section and materials
section_type = "Rectangular"  #@param ["Rectangular", "T-beam"]
b = 12.0  #@param {type:"number"}
h = 24.0  #@param {type:"number"}
bf = 48.0  #@param {type:"number"}
hf = 4.0  #@param {type:"number"}
#@markdown Dimensions in **in**. `bf` and `hf` (effective flange width and thickness, ACI 6.3.2) are used only for a T-beam.
fc = 4000.0  #@param {type:"number"}
fy = 60000.0  #@param {type:"number"}
fyt = 60000.0  #@param {type:"number"}
concrete_weight = 150.0  #@param {type:"number"}
lam = 1.0  #@param {type:"number"}
#@markdown Strengths in **psi**, unit weight in **pcf**. `lam` = lightweight factor (1.0 normalweight).

# %%
#@title 2. Reinforcement
n_bot = 3  #@param {type:"integer"}
bar_bot = 8  #@param [3, 4, 5, 6, 7, 8, 9, 10, 11] {type:"raw"}
n_top = 2  #@param {type:"integer"}
bar_top = 5  #@param [3, 4, 5, 6, 7, 8, 9, 10, 11] {type:"raw"}
cover = 1.5  #@param {type:"number"}
stirrup = 3  #@param [3, 4, 5] {type:"raw"}
legs = 2  #@param {type:"integer"}
s_stirrup = 10.0  #@param {type:"number"}
agg = 0.75  #@param {type:"number"}
#@markdown Clear cover to stirrup, stirrup spacing (0 = none), and max aggregate size (`agg`) in **in**. One bar layer per face.

# %%
#@title 3. Geometry
span_type = "Simple span"  #@param ["Simple span", "Cantilever", "Overhang one end (right)", "Overhang both ends", "Fixed-fixed", "Fixed-pinned (propped)", "Custom supports"]
span = 20.0  #@param {type:"number"}
overhang_left = 0.0  #@param {type:"number"}
overhang_right = 0.0  #@param {type:"number"}
custom_supports = "0:pin, 20:pin, 40:pin"  #@param {type:"string"}
#@markdown Lengths in **ft**. *custom_supports* applies only when *span_type* is **Custom supports** (`x:pin` or `x:fixed`).

# %%
#@title 4. Loads (service level, + down)
#@markdown Full-length uniform loads in **plf** (self weight is added to D when checked):
w_D = 1000.0  #@param {type:"number"}
w_L = 1200.0  #@param {type:"number"}
w_Lr = 0.0  #@param {type:"number"}
w_S = 0.0  #@param {type:"number"}
self_weight = True  #@param {type:"boolean"}
#@markdown Extra loads (separate with `;`): `P, x, kip, case` · `w, x1, x2, kip/ft, case` · `t, x1, x2, w1, w2, case` · `M, x, kip-ft, case`
extra_loads = ""  #@param {type:"string"}

# %%
#@title 5. Deflection
deflection_method = "Effective Ie (ACI 24.2.3.5)"  #@param ["Effective Ie (ACI 24.2.3.5)", "Gross Ig"]
sustained_live = 0.25  #@param {type:"number"}
xi = 2.0  #@param {type:"number"}
live_limit = 360  #@param {type:"integer"}
long_term_limit = 240  #@param [240, 480] {type:"raw"}
#@markdown `sustained_live` = fraction of live load that is sustained. `xi` = time-dependent factor (2.0 for 5 years or more, ACI 24.2.4.1.3).
#@markdown Long-term limit: L/240 when nonstructural elements are not likely damaged, L/480 when they are (ACI Table 24.2.2).

# %%
#@title 6. Calculate
from structcalc.calcs import concrete_beam

out = concrete_beam.run(
    section_type=section_type, b=b, h=h, bf=bf, hf=hf, fc=fc, fy=fy, fyt=fyt, lam=lam,
    wc=concrete_weight, n_bot=n_bot, bar_bot=bar_bot, n_top=n_top, bar_top=bar_top,
    cover=cover, stirrup=stirrup, legs=legs, s_stirrup=s_stirrup, agg=agg,
    span_type=span_type, span=span, overhang_left=overhang_left,
    overhang_right=overhang_right, custom_supports=custom_supports,
    w={"D": w_D / 1000, "L": w_L / 1000, "Lr": w_Lr / 1000, "S": w_S / 1000},
    extra_loads=extra_loads, self_weight=self_weight,
    deflection_method=deflection_method, sustained_live=sustained_live, xi=xi,
    live_limit=live_limit, lt_limit=long_term_limit,
)
out["report"].show()

# %%
#@title 7. Loads and diagrams (governing positive-moment combination)
from structcalc.beam import plot

name = out["governing"]["pos"]
res = next(r for n, _, r in out["combos"] if n == name)
fig = plot(res, title=f"{name} (deflection with gross Ig, factored loads)", beam=out["beam"])

# %%
#@title 8. Save report as HTML (optional)
save_report = False  #@param {type:"boolean"}
if save_report:
    path = out["report"].save_html("concrete_beam.html")
    print("saved", path)
    try:
        from google.colab import files
        files.download(path)
    except ImportError:
        pass
