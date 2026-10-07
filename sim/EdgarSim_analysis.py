#!/usr/bin/env python3
"""
Simplified background analysis (extended: + binned scan vs Z)
"""

import awkward as ak
import numpy as np
import matplotlib.pyplot as plt
import os

plt.rcParams["font.size"] = 14

# ==============================================================================
# Configuration
# ==============================================================================

Q_bb             = 2039  # keV
window_half_width = 55   # keV
ROI_low          = Q_bb - window_half_width  # 1984 keV
ROI_high         = Q_bb + window_half_width  # 2094 keV

npe_threshold         = 4      # photoelectrons from SiPMs
aoe_class_threshold   = -1.80  # from aoe_class_paras.yaml
m1_energy_threshold   = 5      # keV

print("="*80)
print("Simplified Background Analysis in ROI around Q_ββ")
print("="*80)
print(f"Q_ββ = {Q_bb} keV")
print(f"ROI  = [{ROI_low}, {ROI_high}] keV")
print(f"AC UAr threshold : NPE ≤ {npe_threshold}")
print(f"PSD cut          : AoE_class > {aoe_class_threshold}")
print(f"M1               : exactly 1 Ge detector with energy > {m1_energy_threshold} keV")
print("="*80)

# ==============================================================================
# Load data
# ==============================================================================

base_path  = "/lfs/l9/legend/users/esanchez/Simulations/L1000_simulations"
input_file = f"{base_path}/output//Reentrant_tube_noLight/merged_lh5/Bi214_EFCu_RT_nolayer_optical_map_vtx_update.parquet"

print("\nLoading data...")
sim = ak.from_parquet(input_file)
print(f"Total events simulated: {len(sim):,}")

output_dir = f"{base_path}/plots/Reentrant_tube_noLight/update/Bi214"
os.makedirs(output_dir, exist_ok=True)

# ==============================================================================
# Build cut masks
# ==============================================================================

print("\nBuilding cut masks...")

# Vertex z cut: event-level mask
z_cut_pass = sim.vtx_z < 2

# PSD: per-hit mask
psd_pass = sim.AoE_class > aoe_class_threshold

# AC UAr: event-level mask
total_npe_per_event = ak.sum(sim.npe, axis=-1)
uar_veto_pass       = total_npe_per_event <= npe_threshold

# M1: event-level mask
n_ge_fired = ak.sum(sim.energy > m1_energy_threshold, axis=-1)
m1_pass    = n_ge_fired == 1

# Combined masks (all include z cut)
m1_pass     = m1_pass     & z_cut_pass
uar_veto_pass = uar_veto_pass & z_cut_pass
m1_uar_pass = m1_pass & uar_veto_pass

# ==============================================================================
# Helper
# ==============================================================================

def count_hits_in_roi(energies, roi_low, roi_high, mask=None):
    if mask is not None:
        energies = energies[mask]
    flat = ak.flatten(energies)
    return int(ak.sum((flat >= roi_low) & (flat <= roi_high)))

def count_events_in_roi(energies, roi_low, roi_high, mask=None):
    if mask is not None:
        energies = energies[mask]
    return int(ak.sum(ak.any((energies >= roi_low) & (energies <= roi_high), axis=-1)))

def survival_fraction_error(n_pass, n_total):
    if n_total == 0:
        return 0.0
    f = n_pass / n_total
    return np.sqrt(f * (1 - f) / n_total)

# ==============================================================================
# PLOT 1: Ge spectra — M1 / M1+UAr / M1+UAr+PSD
# ==============================================================================

print("\n1. Plotting Ge spectra with progressive cuts...")

ge_m1          = ak.flatten(sim.energy[m1_pass])
ge_m1          = ge_m1[ge_m1 > 0]

ge_m1_uar      = ak.flatten(sim.energy[m1_uar_pass])
ge_m1_uar      = ge_m1_uar[ge_m1_uar > 0]

psd_m1_uar     = psd_pass[m1_uar_pass]
ge_m1_uar_psd  = ak.flatten(sim.energy[m1_uar_pass][psd_m1_uar])
ge_m1_uar_psd  = ge_m1_uar_psd[ge_m1_uar_psd > 0]

fig, ax = plt.subplots(figsize=(10, 7))
for data, label, color in [
    (ge_m1,         'M1',                'C0'),
    (ge_m1_uar,     f'M1 + AC UAr (NPE ≤ {npe_threshold})',              'C1'),
    (ge_m1_uar_psd, f'M1 + AC UAr + PSD (AoE_class > {aoe_class_threshold})', 'C3'),
]:
    if len(data) > 0:
        ax.hist(ak.to_numpy(data), bins=300, range=(1100, 3000),
                histtype='step', linewidth=2, label=label, color=color, alpha=0.85)

ax.axvline(ROI_low,  color='red', linestyle='--', linewidth=1.5, alpha=0.7)
ax.axvline(ROI_high, color='red', linestyle='--', linewidth=1.5, alpha=0.7)
ax.axvspan(ROI_low, ROI_high, alpha=0.08, color='red',
           label=f'ROI [{ROI_low}–{ROI_high} keV]')
ax.set_xlabel('Energy [keV]')
ax.set_ylabel('Counts')
ax.set_yscale('log')
ax.set_xlim(1100, 3000)
ax.legend(fontsize=12)
ax.set_title('Ge spectrum — Bi-214 in EFCu (Reentrant tube)')
plt.tight_layout()
plt.savefig(f"{output_dir}/ge_spectra_cuts.png", dpi=150)
print("  Saved: ge_spectra_cuts.png")
plt.close()

# ==============================================================================
# PLOT 2: Vertex position — x vs y
# ==============================================================================

print("\n2. Plotting vertex x vs y...")

vtx_x = ak.to_numpy(sim.vtx_x)
vtx_y = ak.to_numpy(sim.vtx_y)
vtx_z = ak.to_numpy(sim.vtx_z)

valid = np.isfinite(vtx_x) & np.isfinite(vtx_y) & np.isfinite(vtx_z) & (vtx_z < 2)
vtx_x_v = vtx_x[valid]
vtx_y_v = vtx_y[valid]
vtx_z_v = vtx_z[valid]
vtx_r_v = np.sqrt(vtx_x_v**2 + vtx_y_v**2)

fig, ax = plt.subplots(figsize=(8, 7))
h = ax.hist2d(vtx_x_v, vtx_y_v,
              bins=200,
              norm=plt.matplotlib.colors.LogNorm(),
              cmap='viridis')
plt.colorbar(h[3], ax=ax, label='Counts')
ax.set_xlabel('Vertex x [m]')
ax.set_ylabel('Vertex y [m]')
ax.set_title('Primary vertex position — x vs y')
ax.set_aspect('equal')
plt.tight_layout()
plt.savefig(f"{output_dir}/vertex_x_vs_y.png", dpi=150)
print("  Saved: vertex_x_vs_y.png")
plt.close()

# ==============================================================================
# PLOT 3: Vertex position — z vs r  (z on x-axis, r on y-axis)
# ==============================================================================

print("\n3. Plotting vertex z vs r...")

fig, ax = plt.subplots(figsize=(8, 7))
h = ax.hist2d(vtx_z_v, vtx_r_v,       # FIXED: was (vtx_r_v, vtx_z_v)
              bins=200,
              norm=plt.matplotlib.colors.LogNorm(),
              cmap='viridis')
plt.colorbar(h[3], ax=ax, label='Counts')
ax.set_xlabel('Vertex z [m]')          # FIXED: was 'Vertex r [m]'
ax.set_ylabel('Vertex r [m]')          # FIXED: was 'Vertex z [m]'
ax.set_title('Primary vertex position — z vs r')
plt.tight_layout()
plt.savefig(f"{output_dir}/vertex_r_vs_z.png", dpi=150)
print("  Saved: vertex_r_vs_z.png")
plt.close()

# ==============================================================================
# PLOT 3b: Vertex position — x vs z, z on vertical axis
# ==============================================================================

print("\n3b. Plotting vertex x vs z...")

fig, ax = plt.subplots(figsize=(9, 10))
h = ax.hist2d(vtx_x_v, vtx_z_v,
              bins=200,
              norm=plt.matplotlib.colors.LogNorm(),
              cmap='viridis')
plt.colorbar(h[3], ax=ax, label='Counts')
ax.set_xlabel('Vertex x [m]')
ax.set_ylabel('Vertex z [m]')
ax.set_title('Primary vertex position — x vs z')
plt.tight_layout()
plt.savefig(f"{output_dir}/vertex_x_vs_z.png", dpi=150)
print("  Saved: vertex_x_vs_z.png")
plt.close()

# ==============================================================================
# PLOTS 3c-3f: Underground LAr energy deposition vs vertex position
#
# SCOPE NOTE: this merged parquet only has the primary DECAY vertex position
# and the per-event TOTAL underground LAr energy (sim.lar_energy, summed
# over all LAr hits in that event) -- not individual LAr energy-deposit STEP
# positions (those live only in the raw stp files, used separately for a
# different dataset). So this shows "this event's total LAr energy vs. where
# its primary decay happened", not "where in the LAr the energy landed".
# Same z < 2 m fiducial as the rest of this script (via vtx_x_v/y_v/z_v/r_v,
# already filtered above).
# ==============================================================================

print("\n3c. Plotting underground LAr energy vs vertex position...")

lar_energy_per_event_all = ak.to_numpy(ak.sum(sim.lar_energy, axis=-1))
lar_energy_v = lar_energy_per_event_all[valid]  # same fiducial mask as vtx_x_v/y_v/z_v/r_v

# --- Mean/total LAr energy vs vertex z ---
n_bins_lar_z = 60
z_bins_lar = np.linspace(vtx_z_v.min(), vtx_z_v.max(), n_bins_lar_z)
z_bin_centers_lar = 0.5 * (z_bins_lar[:-1] + z_bins_lar[1:])
total_lar_per_zbin, _ = np.histogram(vtx_z_v, bins=z_bins_lar, weights=lar_energy_v)
n_per_zbin_lar, _ = np.histogram(vtx_z_v, bins=z_bins_lar)
with np.errstate(invalid='ignore', divide='ignore'):
    mean_lar_per_zbin = np.where(n_per_zbin_lar > 0, total_lar_per_zbin / n_per_zbin_lar, np.nan)

fig, axes = plt.subplots(1, 2, figsize=(18, 6))
axes[0].plot(z_bin_centers_lar, total_lar_per_zbin, 'o-', color='C0')
axes[0].set_xlabel('Vertex z [m]')
axes[0].set_ylabel('Total LAr energy [keV]')
axes[0].set_title('Total underground LAr energy vs. vertex z')
axes[1].plot(z_bin_centers_lar, mean_lar_per_zbin, 'o-', color='C1')
axes[1].set_xlabel('Vertex z [m]')
axes[1].set_ylabel('Mean LAr energy per event [keV]')
axes[1].set_title('Mean underground LAr energy vs. vertex z')
plt.tight_layout()
plt.savefig(f"{output_dir}/lar_energy_vs_vertex_z.png", dpi=150)
print("    Saved: lar_energy_vs_vertex_z.png")
plt.close()

# --- Mean/total LAr energy vs vertex r ---
print("3d. Plotting underground LAr energy vs vertex r...")
n_bins_lar_r = 40
r_bins_lar = np.linspace(0, vtx_r_v.max(), n_bins_lar_r)
r_bin_centers_lar = 0.5 * (r_bins_lar[:-1] + r_bins_lar[1:])
total_lar_per_rbin, _ = np.histogram(vtx_r_v, bins=r_bins_lar, weights=lar_energy_v)
n_per_rbin_lar, _ = np.histogram(vtx_r_v, bins=r_bins_lar)
with np.errstate(invalid='ignore', divide='ignore'):
    mean_lar_per_rbin = np.where(n_per_rbin_lar > 0, total_lar_per_rbin / n_per_rbin_lar, np.nan)

fig, axes = plt.subplots(1, 2, figsize=(18, 6))
axes[0].plot(r_bin_centers_lar, total_lar_per_rbin, 'o-', color='C0')
axes[0].set_xlabel('Vertex r [m]')
axes[0].set_ylabel('Total LAr energy [keV]')
axes[0].set_title('Total underground LAr energy vs. vertex r')
axes[1].plot(r_bin_centers_lar, mean_lar_per_rbin, 'o-', color='C1')
axes[1].set_xlabel('Vertex r [m]')
axes[1].set_ylabel('Mean LAr energy per event [keV]')
axes[1].set_title('Mean underground LAr energy vs. vertex r')
plt.tight_layout()
plt.savefig(f"{output_dir}/lar_energy_vs_vertex_r.png", dpi=150)
print("    Saved: lar_energy_vs_vertex_r.png")
plt.close()

# --- 2D maps: total LAr energy per vertex-position bin ---
print("3e. Plotting LAr energy map vs vertex position (x-y)...")
fig, ax = plt.subplots(figsize=(9, 8))
h = ax.hist2d(vtx_x_v, vtx_y_v, bins=150, weights=lar_energy_v,
              norm=plt.matplotlib.colors.LogNorm(), cmap='viridis')
plt.colorbar(h[3], ax=ax, label='Total LAr energy [keV]')
ax.set_xlabel('Vertex x [m]')
ax.set_ylabel('Vertex y [m]')
ax.set_aspect('equal')
ax.set_title('LAr energy map vs. vertex position (x-y)')
plt.tight_layout()
plt.savefig(f"{output_dir}/lar_energy_xy_map.png", dpi=150)
print("    Saved: lar_energy_xy_map.png")
plt.close()

print("3f. Plotting LAr energy map vs vertex position (r-z)...")
fig, ax = plt.subplots(figsize=(8, 10))
h = ax.hist2d(vtx_r_v, vtx_z_v, bins=150, weights=lar_energy_v,
              norm=plt.matplotlib.colors.LogNorm(), cmap='viridis')
plt.colorbar(h[3], ax=ax, label='Total LAr energy [keV]')
ax.set_xlabel('Vertex r [m]')
ax.set_ylabel('Vertex z [m]')
ax.set_title('LAr energy map vs. vertex position (r-z)')
plt.tight_layout()
plt.savefig(f"{output_dir}/lar_energy_rz_map.png", dpi=150)
print("    Saved: lar_energy_rz_map.png")
plt.close()

# ==============================================================================
# PLOTS 3g-3l (NEW): Number of ROI events per position bin, x-y / x-z / r-z,
# BEFORE cuts and AFTER cuts (M1 + AC UAr + PSD).
#
# Same 2D-map construction as Plots 3e/3f (hist2d with `weights=`), but here
# the weight is a per-event 0/1 flag ("did this event have a ROI hit") rather
# than an energy value -- so the colorbar directly reads as "number of ROI
# events in this spatial bin", not a summed energy.
#
# "Before cuts": has_roi_hit = any Ge hit landing in [ROI_low, ROI_high],
# with only the same z<2 fiducial cut used throughout this script -- no M1/
# AC UAr/PSD applied. This matches n_no_cut/ev_no_cut used in the summary
# table below (same event selection, just resolved spatially here).
#
# "After cuts": has_roi_hit_after_cuts = a ROI hit that ALSO passes PSD
# (per-hit) on an event that ALSO passes M1 + AC UAr (event-level) -- the
# same full selection as n_m1_uar_psd/ev_m1_uar_psd in the summary table.
#
# Both flags are computed once here (indexed like the full sim, i.e. before
# the `valid` fiducial mask), then sliced with the same `valid` mask used
# for vtx_x_v/y_v/z_v/r_v so every array here has consistent length/order.
# ==============================================================================

print("\n3g. Computing ROI-hit flags (before and after cuts)...")

has_roi_hit_full = ak.to_numpy(
    ak.any((sim.energy >= ROI_low) & (sim.energy <= ROI_high), axis=-1)
).astype(float)

roi_psd_hit_mask = (sim.energy >= ROI_low) & (sim.energy <= ROI_high) & psd_pass
has_roi_hit_after_cuts_full = (
    ak.to_numpy(ak.any(roi_psd_hit_mask, axis=-1)) & ak.to_numpy(m1_uar_pass)
).astype(float)

# Slice to the same fiducial ("valid") mask as vtx_x_v/y_v/z_v/r_v above.
roi_before_v = has_roi_hit_full[valid]
roi_after_v  = has_roi_hit_after_cuts_full[valid]

print(f"  ROI events before cuts (fiducial z<2): {int(roi_before_v.sum()):,}")
print(f"  ROI events after  cuts (fiducial z<2): {int(roi_after_v.sum()):,}")


def plot_roi_count_map(x, y, weights, xlabel, ylabel, title, filename, equal_aspect=False,
                        figsize=(9, 8)):
    """hist2d of ROI-event counts (weights is a 0/1 flag array) vs. two position
    coordinates. LogNorm colorbar, matching the LAr-energy-map plots' style."""
    fig, ax = plt.subplots(figsize=figsize)
    h = ax.hist2d(x, y, bins=150, weights=weights,
                  norm=plt.matplotlib.colors.LogNorm(), cmap='viridis')
    plt.colorbar(h[3], ax=ax, label='Events with ROI hit')
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if equal_aspect:
        ax.set_aspect('equal')
    ax.set_title(title)
    plt.tight_layout()
    plt.savefig(f"{output_dir}/{filename}", dpi=150)
    print(f"    Saved: {filename}")
    plt.close()


print("3h. Plotting ROI-event map, x-y, before cuts...")
plot_roi_count_map(vtx_x_v, vtx_y_v, roi_before_v,
                    'Vertex x [m]', 'Vertex y [m]',
                    'ROI events (before cuts) vs. vertex position (x-y)',
                    'roi_events_xy_before_cuts.png', equal_aspect=True)

print("3i. Plotting ROI-event map, x-z, before cuts...")
plot_roi_count_map(vtx_x_v, vtx_z_v, roi_before_v,
                    'Vertex x [m]', 'Vertex z [m]',
                    'ROI events (before cuts) vs. vertex position (x-z)',
                    'roi_events_xz_before_cuts.png', figsize=(9, 10))

print("3j. Plotting ROI-event map, r-z, before cuts...")
plot_roi_count_map(vtx_r_v, vtx_z_v, roi_before_v,
                    'Vertex r [m]', 'Vertex z [m]',
                    'ROI events (before cuts) vs. vertex position (r-z)',
                    'roi_events_rz_before_cuts.png', figsize=(8, 10))

print("3k. Plotting ROI-event map, x-y, after cuts (M1 + AC UAr + PSD)...")
plot_roi_count_map(vtx_x_v, vtx_y_v, roi_after_v,
                    'Vertex x [m]', 'Vertex y [m]',
                    'ROI events (M1 + AC UAr + PSD) vs. vertex position (x-y)',
                    'roi_events_xy_after_cuts.png', equal_aspect=True)

print("3l. Plotting ROI-event map, x-z, after cuts (M1 + AC UAr + PSD)...")
plot_roi_count_map(vtx_x_v, vtx_z_v, roi_after_v,
                    'Vertex x [m]', 'Vertex z [m]',
                    'ROI events (M1 + AC UAr + PSD) vs. vertex position (x-z)',
                    'roi_events_xz_after_cuts.png', figsize=(9, 10))

print("3m. Plotting ROI-event map, r-z, after cuts (M1 + AC UAr + PSD)...")
plot_roi_count_map(vtx_r_v, vtx_z_v, roi_after_v,
                    'Vertex r [m]', 'Vertex z [m]',
                    'ROI events (M1 + AC UAr + PSD) vs. vertex position (r-z)',
                    'roi_events_rz_after_cuts.png', figsize=(8, 10))

# ==============================================================================
# PLOTS 4-5: Probability of Ge energy deposition vs vertex position (z and r)
# Three curves per plot:
#   - P(any Ge hit)        : at least one Ge detector with energy > 0
#   - P(ROI hit)           : at least one Ge hit in [ROI_low, ROI_high] keV
#   - P(ROI | all cuts)    : ROI hit surviving M1 + AC UAr + PSD
# ==============================================================================

print("\n4. Computing detection probability profiles...")

vtx_z_all = ak.to_numpy(sim.vtx_z)
vtx_x_all = ak.to_numpy(sim.vtx_x)
vtx_y_all = ak.to_numpy(sim.vtx_y)
vtx_r_all = np.sqrt(vtx_x_all**2 + vtx_y_all**2)

# Event-level boolean flags
has_ge_hit   = ak.to_numpy(ak.sum(sim.energy > 0,       axis=-1) > 0)
has_roi_hit  = ak.to_numpy(ak.sum((sim.energy >= ROI_low) & (sim.energy <= ROI_high), axis=-1) > 0)

# Apply z < 2 fiducial to all
z_mask = vtx_z_all < 2

def efficiency_profile(coord, condition, z_mask, bins):
    """Compute efficiency = N(condition & bin) / N(bin) with Poisson uncertainty."""
    coord_m   = coord[z_mask]
    cond_m    = condition[z_mask]
    n_total,  edges = np.histogram(coord_m,          bins=bins)
    n_pass,   _     = np.histogram(coord_m[cond_m],  bins=edges)
    centers   = 0.5 * (edges[:-1] + edges[1:])
    with np.errstate(invalid='ignore', divide='ignore'):
        eff = np.where(n_total > 0, n_pass / n_total, np.nan)
        err = np.where(n_total > 0, np.sqrt(eff * (1 - eff) / n_total), np.nan)
    return centers, eff, err

n_bins_z = 60
n_bins_r = 40
z_bins   = np.linspace(vtx_z_all[z_mask].min(), 2, n_bins_z)
r_bins   = np.linspace(0, vtx_r_all[z_mask].max(), n_bins_r)


# --- Plot 4a/4b: P(any Ge hit | M1 + AC UAr + PSD) vs z and r ---
print("  4a. P(any Ge hit after cuts) vs vertex z...")
fig, ax = plt.subplots(figsize=(9, 5))
centers, eff, err = efficiency_profile(vtx_z_all, has_ge_hit & ak.to_numpy(m1_uar_pass) & ak.to_numpy(ak.any(psd_pass, axis=-1)), z_mask, z_bins)
ax.plot(centers, eff, color='C0', linewidth=2, label='Any Ge hit | M1 + AC UAr + PSD')
ax.fill_between(centers, eff - err, eff + err, alpha=0.2, color='C0')
ax.set_xlabel('Vertex z [m]')
ax.set_ylabel('Probability')
ax.set_title('P(any Ge hit | M1 + AC UAr + PSD) vs vertex z — Bi-214 EFCu')
ax.legend(fontsize=11)
ax.set_ylim(bottom=0)
plt.tight_layout()
plt.savefig(f"{output_dir}/prob_ge_cuts_vs_z.png", dpi=150)
print("    Saved: prob_ge_cuts_vs_z.png")
plt.close()

print("  4b. P(any Ge hit after cuts) vs vertex r...")
fig, ax = plt.subplots(figsize=(9, 5))
centers, eff, err = efficiency_profile(vtx_r_all, has_ge_hit & ak.to_numpy(m1_uar_pass) & ak.to_numpy(ak.any(psd_pass, axis=-1)), z_mask, r_bins)
ax.plot(centers, eff, color='C0', linewidth=2, label='Any Ge hit | M1 + AC UAr + PSD')
ax.fill_between(centers, eff - err, eff + err, alpha=0.2, color='C0')
ax.set_xlabel('Vertex r [m]')
ax.set_ylabel('Probability')
ax.set_title('P(any Ge hit | M1 + AC UAr + PSD) vs vertex r — Bi-214 EFCu')
ax.legend(fontsize=11)
ax.set_ylim(bottom=0)
plt.tight_layout()
plt.savefig(f"{output_dir}/prob_ge_cuts_vs_r.png", dpi=150)
print("    Saved: prob_ge_cuts_vs_r.png")
plt.close()

# --- Plot 5a/5b: P(ROI hit, no cuts) vs z and r ---
print("  5a. P(ROI hit, no cuts) vs vertex z...")
fig, ax = plt.subplots(figsize=(9, 5))
centers, eff, err = efficiency_profile(vtx_z_all, has_roi_hit, z_mask, z_bins)
ax.plot(centers, eff, color='C1', linewidth=2, label=f'ROI hit [{ROI_low}–{ROI_high} keV], no cuts')
ax.fill_between(centers, eff - err, eff + err, alpha=0.2, color='C1')
ax.set_xlabel('Vertex z [m]')
ax.set_ylabel('Probability')
ax.set_title(f'P(ROI hit, no cuts) vs vertex z — Bi-214 EFCu')
ax.legend(fontsize=11)
ax.set_ylim(bottom=0)
plt.tight_layout()
plt.savefig(f"{output_dir}/prob_roi_nocuts_vs_z.png", dpi=150)
print("    Saved: prob_roi_nocuts_vs_z.png")
plt.close()

print("  5b. P(ROI hit, no cuts) vs vertex r...")
fig, ax = plt.subplots(figsize=(9, 5))
centers, eff, err = efficiency_profile(vtx_r_all, has_roi_hit, z_mask, r_bins)
ax.plot(centers, eff, color='C1', linewidth=2, label=f'ROI hit [{ROI_low}–{ROI_high} keV], no cuts')
ax.fill_between(centers, eff - err, eff + err, alpha=0.2, color='C1')
ax.set_xlabel('Vertex r [m]')
ax.set_ylabel('Probability')
ax.set_title(f'P(ROI hit, no cuts) vs vertex r — Bi-214 EFCu')
ax.legend(fontsize=11)
ax.set_ylim(bottom=0)
plt.tight_layout()
plt.savefig(f"{output_dir}/prob_roi_nocuts_vs_r.png", dpi=150)
print("    Saved: prob_roi_nocuts_vs_r.png")
plt.close()

# ==============================================================================
# Summary table
# ==============================================================================

ge_no_cut    = sim.energy[z_cut_pass]
n_no_cut     = count_hits_in_roi(ge_no_cut,       ROI_low, ROI_high)
ev_no_cut    = count_events_in_roi(ge_no_cut,      ROI_low, ROI_high)

n_m1         = count_hits_in_roi(sim.energy[m1_pass],     ROI_low, ROI_high)
ev_m1        = count_events_in_roi(sim.energy[m1_pass],   ROI_low, ROI_high)

n_m1_uar     = count_hits_in_roi(sim.energy[m1_uar_pass],   ROI_low, ROI_high)
ev_m1_uar    = count_events_in_roi(sim.energy[m1_uar_pass], ROI_low, ROI_high)

n_m1_uar_psd = count_hits_in_roi(sim.energy[m1_uar_pass],   ROI_low, ROI_high, mask=psd_m1_uar)
ev_m1_uar_psd= count_events_in_roi(sim.energy[m1_uar_pass], ROI_low, ROI_high, mask=psd_m1_uar)

# Additional scenarios for full table (all include z cut)
psd_pass_only       = psd_pass & z_cut_pass
uar_pass_sim        = sim.energy[uar_veto_pass]   # uar_veto_pass already includes z_cut
psd_pass_uar        = psd_pass[uar_veto_pass]

n_psd_only          = count_hits_in_roi(sim.energy,    ROI_low, ROI_high, mask=psd_pass_only)
ev_psd_only         = count_events_in_roi(sim.energy,  ROI_low, ROI_high, mask=psd_pass_only)

n_uar_only          = count_hits_in_roi(uar_pass_sim,  ROI_low, ROI_high)
ev_uar_only         = count_events_in_roi(uar_pass_sim, ROI_low, ROI_high)

n_uar_psd           = count_hits_in_roi(uar_pass_sim,  ROI_low, ROI_high, mask=psd_pass_uar)
ev_uar_psd          = count_events_in_roi(uar_pass_sim, ROI_low, ROI_high, mask=psd_pass_uar)

def surv(n):
    return n / n_no_cut * 100 if n_no_cut > 0 else 0.0

def surv_err(n):
    return survival_fraction_error(n, n_no_cut) * 100

# ==============================================================================
# Detailed per-scenario printout (matching the other Reentrant_tube script's
# style): total ROI hits/events before cuts, then after each individual and
# combined cut, with survival % (+/- binomial error) and rejection %.
# ==============================================================================

print("\n" + "="*80)
print("SCENARIO 1: NO CUTS")
print("="*80)
print(f"Ge hits in ROI [{ROI_low}-{ROI_high} keV]: {n_no_cut:,}")
print(f"Events with Ge hits in ROI: {ev_no_cut:,}")

print("\n" + "="*80)
print(f"SCENARIO 2: PSD ONLY (AoE_class > {aoe_class_threshold})")
print("="*80)
print(f"Ge hits in ROI: {n_psd_only:,}  |  Events: {ev_psd_only:,}")
print(f"Survival: {surv(n_psd_only):.2f} ± {surv_err(n_psd_only):.2f}%  |  "
      f"Rejection: {100-surv(n_psd_only):.2f}%")

print("\n" + "="*80)
print("SCENARIO 3: AC UAr ONLY")
print("="*80)
print(f"Ge hits in ROI: {n_uar_only:,}  |  Events: {ev_uar_only:,}")
print(f"Survival: {surv(n_uar_only):.2f} ± {surv_err(n_uar_only):.2f}%  |  "
      f"Rejection: {100-surv(n_uar_only):.2f}%")

print("\n" + "="*80)
print("SCENARIO 4: M1 ONLY (single Ge hit)")
print("="*80)
print(f"Ge hits in ROI: {n_m1:,}  |  Events: {ev_m1:,}")
print(f"Survival: {surv(n_m1):.2f} ± {surv_err(n_m1):.2f}%  |  "
      f"Rejection: {100-surv(n_m1):.2f}%")

print("\n" + "="*80)
print("SCENARIO 5: M1 + AC UAr")
print("="*80)
print(f"Ge hits in ROI: {n_m1_uar:,}  |  Events: {ev_m1_uar:,}")
print(f"Survival: {surv(n_m1_uar):.2f} ± {surv_err(n_m1_uar):.2f}%  |  "
      f"Rejection: {100-surv(n_m1_uar):.2f}%")

print("\n" + "="*80)
print("SCENARIO 6: AC UAr + PSD")
print("="*80)
print(f"Ge hits in ROI: {n_uar_psd:,}  |  Events: {ev_uar_psd:,}")
print(f"Survival: {surv(n_uar_psd):.2f} ± {surv_err(n_uar_psd):.2f}%  |  "
      f"Rejection: {100-surv(n_uar_psd):.2f}%")

print("\n" + "="*80)
print("SCENARIO 7: M1 + AC UAr + PSD")
print("="*80)
print(f"Ge hits in ROI: {n_m1_uar_psd:,}  |  Events: {ev_m1_uar_psd:,}")
print(f"Survival: {surv(n_m1_uar_psd):.2f} ± {surv_err(n_m1_uar_psd):.2f}%  |  "
      f"Rejection: {100-surv(n_m1_uar_psd):.2f}%")

rows = [
    ("No cuts",           n_no_cut,     ev_no_cut,     100.0,               0.0,                 0.0),
    ("PSD only",          n_psd_only,   ev_psd_only,   surv(n_psd_only),    surv_err(n_psd_only),    100-surv(n_psd_only)),
    ("AC UAr only",       n_uar_only,   ev_uar_only,   surv(n_uar_only),    surv_err(n_uar_only),    100-surv(n_uar_only)),
    ("M1 only",           n_m1,         ev_m1,         surv(n_m1),          surv_err(n_m1),          100-surv(n_m1)),
    ("M1 + AC UAr",       n_m1_uar,     ev_m1_uar,     surv(n_m1_uar),      surv_err(n_m1_uar),      100-surv(n_m1_uar)),
    ("AC UAr + PSD",      n_uar_psd,    ev_uar_psd,    surv(n_uar_psd),     surv_err(n_uar_psd),     100-surv(n_uar_psd)),
    ("M1 + AC UAr + PSD", n_m1_uar_psd, ev_m1_uar_psd, surv(n_m1_uar_psd), surv_err(n_m1_uar_psd), 100-surv(n_m1_uar_psd)),
]

print("\n" + "="*95)
print("SUMMARY TABLE")
print("="*95)
print(f"{'Cut Scenario':<35} {'ROI Hits':<12} {'Events in Q_ββ':<18} {'Survival %':<22} {'Rejection %'}")
print("-"*95)
for label, hits, events, s, s_err, r in rows:
    print(f"{label:<35} {hits:<12,} {events:<18,} {s:.2f} ± {s_err:.2f}{'%':<10} {r:.2f}%")
print("="*95)

print("\nOriginal analysis complete. Continuing with binned scan vs Z...")
print("="*80)

# ==============================================================================
# NEW: Binned scan vs Z (50 cm bins)
#
# Unlike the probability PROFILES above (Plots 4-5, which use fine, roughly-
# equal-count bins purely for a smooth-looking curve), this is an explicit
# fixed-width (50 cm) bin-by-bin table+plot, giving N_gen, ROI hits before/
# after the FULL cut chain (M1 + AC UAr + PSD, reusing m1_uar_pass and
# psd_pass already computed above), survival fraction, and an absolute BI
# rate per bin -- directly comparable in structure to the non-optical-map
# EFCu/Reentrant_tube "BI vs Z" analysis.
#
# ASSUMPTIONS:
#   - Same z < 2 m fiducial cut used throughout this script is applied here
#     too, for consistency with everything above.
#   - Activity/mass: "Detector Mount (EFCu)" from the BOM table (matches
#     this dataset's source): 238U = (0.19 +/- 0.10) uBq/kg, mass = 16.2 kg.
#     Only 238U is used (214Bi is a 238U-chain daughter in secular
#     equilibrium). Update EFCU_ACTIVITY_UBQ_PER_KG / EFCU_MASS_KG below if
#     this isn't the right component/mass for this optical-map dataset.
#   - BI vs Z is an ABSOLUTE rate (Bq / counts-per-year), NOT normalized to
#     cts/(keV*kg*yr) -- that needs the Ge active mass and a livetime.
# ==============================================================================

print("\n" + "="*80)

BIN_WIDTH_Z_M = 0.10  # 50 cm -- change this to rebin

print(f"Binned scan vs Z ({BIN_WIDTH_Z_M*100:.0f} cm bins, M1 + AC UAr + PSD)")
print("="*80)

EFCU_ACTIVITY_UBQ_PER_KG = 0.19
EFCU_ACTIVITY_ERR_UBQ_PER_KG = 0.10
EFCU_MASS_KG = 16.2

SECONDS_PER_YEAR = 365.25 * 24 * 3600
N_total_generated = len(sim)

print(f"Z bin width: {BIN_WIDTH_Z_M*100:.0f} cm")
print(f"Fiducial cut applied: vtx_z < 2 m (same as rest of script)")
print(f"Cut used for 'after cuts': M1 + AC UAr + PSD")

# --- Overall (no Z binning) survival, for reference -- same as n_m1_uar_psd above ---
n_roi_no_cut_all_z = n_no_cut  # already computed above, with z_cut_pass applied
n_roi_full_all_z   = n_m1_uar_psd
surv_overall_z     = surv(n_m1_uar_psd)
surv_overall_z_err = survival_fraction_error(n_roi_full_all_z, n_roi_no_cut_all_z) * 100
print(f"\nOverall (no Z binning): ROI hits before cuts = {n_roi_no_cut_all_z:,}, "
      f"after M1+AC UAr+PSD = {n_roi_full_all_z:,}")
print(f"Overall survival: {surv_overall_z:.2f} +/- {surv_overall_z_err:.2f}%")

# --- Bin by vtx_z (respecting the same z < 2 fiducial used throughout) ---
vtx_z_full_z = ak.to_numpy(sim.vtx_z)  # meters, signed
valid_z_z = np.isfinite(vtx_z_full_z) & (vtx_z_full_z < 2)
z_min_z = float(np.min(vtx_z_full_z[valid_z_z]))
z_max_z = float(np.max(vtx_z_full_z[valid_z_z]))
print(f"\nvtx_z range (z<2 fiducial): [{z_min_z:.3f}, {z_max_z:.3f}] m")

n_bins_z_scan = int(np.ceil((z_max_z - z_min_z) / BIN_WIDTH_Z_M))
bin_edges_z_m = z_min_z + BIN_WIDTH_Z_M * np.arange(n_bins_z_scan + 1)
bin_centers_z_cm = (bin_edges_z_m[:-1] + bin_edges_z_m[1:]) / 2.0 * 100.0
bin_labels_z_cm = [f"{bin_edges_z_m[i]*100:+.0f} to {bin_edges_z_m[i+1]*100:+.0f}"
                    for i in range(n_bins_z_scan)]

print(f"Number of Z bins: {n_bins_z_scan}")

n_gen_per_zbin        = np.zeros(n_bins_z_scan, dtype=int)
n_roi_no_cut_per_zbin = np.zeros(n_bins_z_scan, dtype=int)
n_roi_full_per_zbin   = np.zeros(n_bins_z_scan, dtype=int)

for i in range(n_bins_z_scan):
    lo, hi = bin_edges_z_m[i], bin_edges_z_m[i + 1]
    in_bin = valid_z_z & (vtx_z_full_z >= lo) & \
             (vtx_z_full_z < hi if i < n_bins_z_scan - 1 else vtx_z_full_z <= hi)

    n_gen_per_zbin[i] = int(np.sum(in_bin))
    n_roi_no_cut_per_zbin[i] = int(count_events_in_roi(
        sim.energy[in_bin], ROI_low, ROI_high
    ))

    # "after cuts" = M1 + AC UAr + PSD, restricted to this bin
    m1_uar_bin = ak.to_numpy(m1_uar_pass)[in_bin]
    sim_bin = sim[in_bin]
    sel_bin = sim_bin[m1_uar_bin]
    # psd_pass is per-hit (jagged, variable length per event) -- slice it
    # directly with awkward indexing, never ak.to_numpy() it (that crashes:
    # "cannot convert to RegularArray" on a ragged boolean array).
    psd_sel_bin = psd_pass[in_bin][m1_uar_bin]
    n_roi_full_per_zbin[i] = int(count_events_in_roi(
        sel_bin.energy, ROI_low, ROI_high, mask=psd_sel_bin
    ))

surv_per_zbin = np.divide(n_roi_full_per_zbin, n_roi_no_cut_per_zbin,
                           out=np.zeros(n_bins_z_scan), where=n_roi_no_cut_per_zbin > 0) * 100
surv_err_per_zbin = np.array([
    survival_fraction_error(a, b) * 100
    for a, b in zip(n_roi_full_per_zbin, n_roi_no_cut_per_zbin)
])

print("\nbin [cm]              N_gen    ROI(no cuts)   ROI(M1+AC UAr+PSD)   survival %")
print("-" * 85)
for i in range(n_bins_z_scan):
    print(f"{bin_labels_z_cm[i]:<20} {n_gen_per_zbin[i]:<8,} {n_roi_no_cut_per_zbin[i]:<15,} "
          f"{n_roi_full_per_zbin[i]:<20,} {surv_per_zbin[i]:.2f} +/- {surv_err_per_zbin[i]:.2f}")

# --- Plots: ROI counts vs Z, survival fraction vs Z ---
fig, ax = plt.subplots(figsize=(11, 7))
ax.errorbar(bin_centers_z_cm, n_roi_no_cut_per_zbin, yerr=np.sqrt(n_roi_no_cut_per_zbin),
            fmt='o-', label='No cuts', color='C0', capsize=3)
ax.errorbar(bin_centers_z_cm, n_roi_full_per_zbin, yerr=np.sqrt(n_roi_full_per_zbin),
            fmt='s-', label='M1 + AC UAr + PSD', color='C3', capsize=3)
ax.set_xlabel('vertex z [cm]')
ax.set_ylabel('Events with a Ge hit in ROI')
ax.set_yscale('log')
ax.legend()
plt.tight_layout()
plt.savefig(f"{output_dir}/scan_roi_counts_vs_z.png", dpi=150)
print(f"\nSaved: {output_dir}/scan_roi_counts_vs_z.png")
plt.close()

# --- Standalone: background BEFORE cuts only, vs Z ---
fig, ax = plt.subplots(figsize=(11, 7))
ax.errorbar(bin_centers_z_cm, n_roi_no_cut_per_zbin, yerr=np.sqrt(n_roi_no_cut_per_zbin),
            fmt='o-', color='C0', capsize=4, markersize=7)
ax.set_xlabel('vertex z [cm]')
ax.set_ylabel('Events with a Ge hit in ROI, before cuts\n(raw simulated counts)')
plt.tight_layout()
plt.savefig(f"{output_dir}/scan_roi_counts_no_cuts_vs_z.png", dpi=150)
print(f"Saved: {output_dir}/scan_roi_counts_no_cuts_vs_z.png")
plt.close()

# --- Standalone: background BEFORE cuts, vs DISTANCE FROM BOTTOM, with the
#     EFCu/OFCH-Cu layer interface marked ---
#
# Reentrant tube wall structure (from the mechanical drawing), bottom to top:
#   494mm  bottom cap (same EFCu material)
#   1685mm EFCu, straight wall  <- this dataset's source material
#   2005mm OFCH Cu
#   2062mm 316L SS
#   = 6247mm total
#
# Since this simulation's source is EFCu, the relevant interface is where
# the EFCu layer ends and OFCH Cu begins: 1685mm above the bottom of the
# straight EFCu wall. Our own vtx_z sampling only spans the EFCu confinement
# region (a few meters, not the full 6.247m tube), so "bottom" here is taken
# as the lowest OBSERVED vertex in the data (z_min_z) -- i.e. this assumes
# the confinement's own lower edge coincides with the bottom of the straight
# EFCu wall. If instead it coincides with the very bottom of the dome
# (including the 494mm cap), the interface would sit at 494+1685=2179mm
# instead of 1685mm -- check against the geometry if this matters.
dist_from_bottom_cm = bin_centers_z_cm - z_min_z * 100.0
EFCU_OFCH_INTERFACE_CM = 168.5 + 49.4  # 1685 mm, per the drawing

fig, ax = plt.subplots(figsize=(11, 7))
ax.errorbar(dist_from_bottom_cm, n_roi_no_cut_per_zbin, yerr=np.sqrt(n_roi_no_cut_per_zbin),
            fmt='o-', color='C0', capsize=4, markersize=7)
ax.axvline(EFCU_OFCH_INTERFACE_CM, color='red', linestyle='--', linewidth=1.5,
           label=f'EFCu / OFCH Cu interface ({EFCU_OFCH_INTERFACE_CM:.0f} cm)')
ax.set_xlabel('Distance from bottom [cm]\n(bottom = lowest observed vertex in this dataset)')
ax.set_ylabel('Events with a Ge hit in ROI, before cuts\n(raw simulated counts)')
ax.legend()
plt.tight_layout()
plt.savefig(f"{output_dir}/scan_roi_counts_no_cuts_vs_bottom_distance.png", dpi=150)
print(f"Saved: {output_dir}/scan_roi_counts_no_cuts_vs_bottom_distance.png")
plt.close()

fig, ax = plt.subplots(figsize=(11, 7))
ax.errorbar(bin_centers_z_cm, surv_per_zbin, yerr=surv_err_per_zbin,
            fmt='o-', color='C0', capsize=4, markersize=7)
ax.axhline(surv_overall_z, color='red', linestyle='--', linewidth=1.5,
           label=f'Overall average ({surv_overall_z:.2f}%)')
ax.fill_between([bin_centers_z_cm.min(), bin_centers_z_cm.max()],
                 surv_overall_z - surv_overall_z_err, surv_overall_z + surv_overall_z_err,
                 color='red', alpha=0.15)
ax.set_xlabel('vertex z [cm]')
ax.set_ylabel('Survival fraction in ROI [%]\n(M1 + AC UAr + PSD vs. no cuts)')
ax.legend()
plt.tight_layout()
plt.savefig(f"{output_dir}/scan_survival_fraction_vs_z.png", dpi=150)
print(f"Saved: {output_dir}/scan_survival_fraction_vs_z.png")
plt.close()

# --- BI vs Z (absolute rate) ---
total_activity_bq_z = EFCU_ACTIVITY_UBQ_PER_KG * 1e-6 * EFCU_MASS_KG
total_activity_err_bq_z = EFCU_ACTIVITY_ERR_UBQ_PER_KG * 1e-6 * EFCU_MASS_KG
print(f"\nTotal EFCu 214Bi activity (238U-chain, secular equilibrium, "
      f"mass {EFCU_MASS_KG} kg): {total_activity_bq_z:.3e} +/- {total_activity_err_bq_z:.3e} Bq")

activity_per_generated_event_z = total_activity_bq_z / N_total_generated
activity_per_generated_event_z_err = total_activity_err_bq_z / N_total_generated
print(f"Activity per simulated decay: {activity_per_generated_event_z:.3e} "
      f"+/- {activity_per_generated_event_z_err:.3e} Bq/event "
      f"(N_total_generated = {N_total_generated:,})")

bi_bq_per_zbin = activity_per_generated_event_z * n_roi_full_per_zbin
bi_bq_err_per_zbin = np.array([
    activity_per_generated_event_z * np.sqrt(n) if n > 0 else 0.0
    for n in n_roi_full_per_zbin
])
bi_cpy_per_zbin = bi_bq_per_zbin * SECONDS_PER_YEAR
bi_cpy_err_per_zbin = bi_bq_err_per_zbin * SECONDS_PER_YEAR

print("\nbin [cm]              BI [Bq]              BI [counts/yr]")
print("-" * 70)
for i in range(n_bins_z_scan):
    print(f"{bin_labels_z_cm[i]:<20} {bi_bq_per_zbin[i]:.3e} +/- {bi_bq_err_per_zbin[i]:.3e}   "
          f"{bi_cpy_per_zbin[i]:.3e} +/- {bi_cpy_err_per_zbin[i]:.3e}")

bi_total_bq_z = activity_per_generated_event_z * n_roi_full_all_z
bi_total_cpy_z = bi_total_bq_z * SECONDS_PER_YEAR
print(f"\nTotal (sum over Z bins): {bi_total_bq_z:.3e} Bq = {bi_total_cpy_z:.3e} counts/yr")

fig, ax = plt.subplots(figsize=(11, 7))
ax.errorbar(bin_centers_z_cm, bi_cpy_per_zbin, yerr=bi_cpy_err_per_zbin,
            fmt='o-', color='C4', capsize=4, markersize=7)
ax.set_xlabel('vertex z [cm]')
ax.set_ylabel('Background rate after M1 + AC UAr + PSD, in ROI [counts/yr]')
ax.set_yscale('log')
plt.tight_layout()
plt.savefig(f"{output_dir}/scan_bi_vs_z.png", dpi=150)
print(f"\nSaved: {output_dir}/scan_bi_vs_z.png")
plt.close()

# --- Save scan results to file ---
scan_results_file = f"{output_dir}/scan_vs_z_results.txt"
with open(scan_results_file, 'w') as f:
    f.write("="*80 + "\n")
    f.write("Binned scan vs Z -- Reentrant_tube_noLight / EFCu optical map, 214Bi\n")
    f.write("Cuts: M1 + AC UAr + PSD\n")
    f.write("="*80 + "\n")
    f.write(f"Z bin width: {BIN_WIDTH_Z_M*100:.0f} cm\n")
    f.write(f"vtx_z range (z<2 fiducial): [{z_min_z:.3f}, {z_max_z:.3f}] m\n")
    f.write(f"Total events: {N_total_generated:,}\n")
    f.write(f"Overall survival (M1+AC UAr+PSD / no cuts): {surv_overall_z:.2f} +/- "
            f"{surv_overall_z_err:.2f}%\n")
    f.write(f"EFCu 214Bi activity: {total_activity_bq_z:.3e} +/- {total_activity_err_bq_z:.3e} Bq "
            f"(238U={EFCU_ACTIVITY_UBQ_PER_KG} uBq/kg, mass={EFCU_MASS_KG} kg)\n\n")
    f.write(f"{'bin [cm]':<20} {'N_gen':<10} {'ROI(no cuts)':<15} {'ROI(M1+AC+PSD)':<17} "
            f"{'survival %':<20} {'BI [counts/yr]'}\n")
    f.write("-"*115 + "\n")
    for i in range(n_bins_z_scan):
        f.write(f"{bin_labels_z_cm[i]:<20} {n_gen_per_zbin[i]:<10} {n_roi_no_cut_per_zbin[i]:<15} "
                f"{n_roi_full_per_zbin[i]:<17} {surv_per_zbin[i]:.2f} +/- {surv_err_per_zbin[i]:<12.2f} "
                f"{bi_cpy_per_zbin[i]:.3e} +/- {bi_cpy_err_per_zbin[i]:.3e}\n")
    f.write(f"\nTotal BI (all Z): {bi_total_bq_z:.3e} Bq = {bi_total_cpy_z:.3e} counts/yr\n")
    f.write("\nNOTE: absolute rate (counts/yr), NOT normalized to cts/(keV*kg*yr) -- \n")
    f.write("that would need the Ge active mass and a livetime/exposure.\n")

print(f"\nScan-vs-Z results saved to: {scan_results_file}")

print("\nDone.")
print("="*80)
