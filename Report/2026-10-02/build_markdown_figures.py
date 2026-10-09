"""Recreate the five notebook diagnostics for the S1, S2, and S3 talk scope."""

from pathlib import Path
import shutil
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from astropy.table import Table
from astropy.coordinates import SkyCoord
import astropy.units as u

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
ASSETS = HERE / "assets"
MATCH = ROOT / "output/COSMOS1/singlebandphot/matching"
TEST_FIG = MATCH / "test_figures"
BASE_FIG = MATCH / "figures"
ASSETS.mkdir(exist_ok=True)

for name in [
    "S1_star_offsets_master_1144.png", "S1_star_offsets_master_246.png", "S1_star_offsets_master_1570.png",
    "S2_all_offsets_master_802.png", "S2_all_offsets_master_176.png", "S2_all_offsets_master_763.png",
    "S3_galaxy_offsets_master_497.png", "S3_galaxy_offsets_master_567.png", "S3_galaxy_offsets_master_1317.png",
    "sdss_point_source_normalized_residuals.png",
]:
    shutil.copy2(TEST_FIG / name, ASSETS / name)
for name in [
    "singleband_input_summary.png", "7dt_matching_summary.png",
    "finalcat_radius_diagnostics.png", "master_population_summary.png",
]:
    shutil.copy2(BASE_FIG / name, ASSETS / name)

full = Table.read(MATCH / "master_full20_sdss.fits")
diag = Table.read(MATCH / "matching_diagnostics_7dt.fits")
bands = [f"m{n}" for n in range(400, 876, 25)]
wave = np.array([int(b[1:]) for b in bands])
pixel_scale = 0.505
selected = diag[np.asarray(diag["selected"], dtype=bool)]
by_master = {}
for row in selected:
    mid = int(row["master_id"])
    if mid not in by_master:
        by_master[mid] = {}
    by_master[mid][str(row["band"]).strip()] = row

sdss_type = np.ma.filled(full["sdss_type"], -1).astype(int)
sdss_clean = np.ma.filled(full["sdss_clean"], 0).astype(int)
sdss_mode = np.ma.filled(full["sdss_mode"], 0).astype(int)
sdss_matched = np.asarray(full["sdss_matched"], dtype=bool)
sdss_candidates = np.asarray(full["sdss_n_candidates"], dtype=int)
clean_unique = sdss_matched & (sdss_clean == 1) & (sdss_mode == 1) & (sdss_candidates == 1)
sample_masks = {
    "S1 star": clean_unique & (sdss_type == 6),
    "S2 all": np.ones(len(full), dtype=bool),
    "S3 galaxy": clean_unique & (sdss_type == 3),
}
colors = {"S1 star": "#1f77b4", "S2 all": "#171717", "S3 galaxy": "#d62728"}

rows_by_master = {}
for row in full:
    mid = int(row["master_id"])
    members = by_master.get(mid, {})
    if set(members) != set(bands):
        raise ValueError(f"master {mid} does not have 20 selected members")
    member_rows = [members[b] for b in bands]
    center = SkyCoord(float(row["master_ra"]) * u.deg, float(row["master_dec"]) * u.deg)
    band_coords = SkyCoord(
        np.array([float(m["source_ra"]) for m in member_rows]) * u.deg,
        np.array([float(m["source_dec"]) for m in member_rows]) * u.deg,
    )
    dra, ddec = center.spherical_offsets_to(band_coords)
    dra, ddec = dra.to_value(u.arcsec), ddec.to_value(u.arcsec)
    sig = np.array([float(m["sigma_band"]) for m in member_rows]) * pixel_scale
    weight = 1 / sig**2
    remaining = weight.sum() - weight
    loo_dra = dra - (np.sum(weight * dra) - weight * dra) / remaining
    loo_ddec = ddec - (np.sum(weight * ddec) - weight * ddec) / remaining
    radius = np.hypot(loo_dra, loo_ddec)
    sigma_pred = np.sqrt(sig**2 + 1 / remaining)
    rows_by_master[mid] = {
        "ra": dra, "dec": ddec, "r": radius, "d": radius / sigma_pred,
        "source_id": np.array([int(m["source_id"]) for m in member_rows]),
    }

snr_by_band = {}
for band in bands:
    raw = Table.read(MATCH.parent / f"cat_cosmos1_{band}.fits")
    snr_by_band[band] = dict(zip(np.asarray(raw["source_id"], dtype=int), np.asarray(raw["snrratio"], dtype=float)))

samples = {}
for name, mask in sample_masks.items():
    ids = np.asarray(full["master_id"], dtype=int)[mask]
    arrays = {k: np.concatenate([rows_by_master[int(mid)][k] for mid in ids]) for k in ["ra", "dec", "r", "d"]}
    arrays["wave"] = np.tile(wave, len(ids))
    arrays["band"] = np.tile(np.array(bands), len(ids))
    arrays["source_id"] = np.concatenate([rows_by_master[int(mid)]["source_id"] for mid in ids])
    arrays["snr"] = np.array([
        snr_by_band[b].get(int(sid), np.nan)
        for b, sid in zip(arrays["band"], arrays["source_id"])
    ])
    samples[name] = arrays
    print(name, len(ids), len(arrays["d"]),
          f"median_D={np.median(arrays['d']):.4f}",
          f"frac_D_lt_3={np.mean(arrays['d'] < 3):.4f}",
          f"median_R={np.median(arrays['r']):.4f}")

plt.rcParams.update({"font.size": 12, "axes.titlesize": 15, "figure.titlesize": 17})

def save(fig, stem):
    fig.tight_layout()
    fig.savefig(ASSETS / f"{stem}.png", dpi=180, bbox_inches="tight")
    plt.close(fig)

# 1. Pooled RA and Dec offsets.
all_offsets = np.concatenate([np.r_[s["ra"], s["dec"]] for s in samples.values()])
bins = np.linspace(*np.percentile(all_offsets, [0.5, 99.5]), 61)
fig, axes = plt.subplots(1, 2, figsize=(14, 4.8), sharey=True)
for name, s in samples.items():
    for ax, key in zip(axes, ["ra", "dec"]):
        ax.hist(s[key], bins=bins, density=True, histtype="step", linewidth=1.6, color=colors[name], label=name)
for ax in axes: ax.axvline(0, color="0.6", lw=0.8)
axes[0].set(title="RA offsets", xlabel=r"$\Delta\alpha\cos\delta$ [arcsec]", ylabel="Density")
axes[1].set(title="Dec offsets", xlabel=r"$\Delta\delta$ [arcsec]")
axes[0].legend(frameon=False)
save(fig, "s1_s2_s3_pooled_delta")

# 2. Band-wise NMAD.
def nmad(x): return 1.4826 * np.median(np.abs(x - np.median(x)))
fig, axes = plt.subplots(1, 2, figsize=(14, 4.8), sharex=True)
for name, s in samples.items():
    for ax, key in zip(axes, ["ra", "dec"]):
        vals = [nmad(s[key][s["band"] == b]) for b in bands]
        ax.plot(wave, vals, "o-", ms=4, color=colors[name], label=name)
axes[0].set(title="RA scatter by band", xlabel="Wavelength [nm]", ylabel="RA NMAD [arcsec]")
axes[1].set(title="Dec scatter by band", xlabel="Wavelength [nm]", ylabel="Dec NMAD [arcsec]")
axes[0].legend(frameon=False)
save(fig, "s1_s2_s3_nmad_vs_wavelength")

# 3. Pooled normalized residuals and the unit Rayleigh curve.
dmax = max(4.5, np.percentile(np.concatenate([s["d"] for s in samples.values()]), 99))
dbins = np.linspace(0, dmax, 61)
dgrid = np.linspace(0, dmax, 500)
fig, ax = plt.subplots(figsize=(10, 4.8))
for name, s in samples.items():
    ax.hist(s["d"], bins=dbins, density=True, histtype="step", linewidth=1.7, color=colors[name], label=f"{name} (n={len(s['d'])})")
ax.plot(dgrid, dgrid*np.exp(-0.5*dgrid**2), "--", color="#2ca02c", label="Unit Rayleigh")
ax.axvline(3, color="0.45", ls=":", label="D = 3")
ax.set(xlabel=r"LOO normalized radial residual $D_{\mathrm{LOO}}$", ylabel="Density", title="Pooled normalized residuals")
ax.legend(frameon=False)
save(fig, "s1_s2_s3_pooled_D")

# 4. Median normalized residual by band.
fig, ax = plt.subplots(figsize=(11, 4.8))
for name, s in samples.items():
    medians = [np.median(s["d"][s["band"] == b]) for b in bands]
    ax.plot(wave, medians, "o-", ms=4, color=colors[name], label=name)
ax.axhline(np.sqrt(2*np.log(2)), color="#2ca02c", ls="--", label="Unit Rayleigh median")
ax.set(xlabel="Wavelength [nm]", ylabel=r"Median $D_{\mathrm{LOO}}$", title="Median normalized residual by band")
ax.legend(frameon=False)
save(fig, "s1_s2_s3_median_D_vs_wavelength")

# 5. SourceXtractor++ S/N and raw LOO radial residual.
all_snr = np.concatenate([s["snr"] for s in samples.values()])
all_r = np.concatenate([s["r"] for s in samples.values()])
valid_all = np.isfinite(all_snr) & (all_snr > 0) & np.isfinite(all_r)
snr_min, snr_max = np.percentile(all_snr[valid_all], [0.2, 99.8])
radius_max = np.percentile(all_r[valid_all], 99.5)
snr_bins = np.geomspace(snr_min, snr_max, 16)
snr_centers = np.sqrt(snr_bins[:-1] * snr_bins[1:])
fig, axes = plt.subplots(1, 3, figsize=(17, 4.9), sharex=True, sharey=True)
for ax, (name, s) in zip(axes, samples.items()):
    valid = np.isfinite(s["snr"]) & (s["snr"] > 0) & np.isfinite(s["r"])
    sc = ax.scatter(s["snr"][valid], s["r"][valid], c=s["wave"][valid], cmap="viridis", vmin=400, vmax=875,
                    s=5, alpha=0.16, linewidths=0, rasterized=True)
    med, lo, hi = [], [], []
    for left, right in zip(snr_bins[:-1], snr_bins[1:]):
        sel = valid & (s["snr"] >= left) & (s["snr"] < right)
        if sel.sum() >= 5:
            med.append(np.median(s["r"][sel]))
            lo.append(np.percentile(s["r"][sel], 16))
            hi.append(np.percentile(s["r"][sel], 84))
        else:
            med.append(np.nan); lo.append(np.nan); hi.append(np.nan)
    med, lo, hi = np.array(med), np.array(lo), np.array(hi)
    good = np.isfinite(med)
    ax.fill_between(snr_centers[good], lo[good], hi[good], color="black", alpha=0.12)
    ax.plot(snr_centers[good], med[good], "o-", color="black", ms=3, lw=1.2)
    ax.set_xscale("log")
    ax.set(xlim=(snr_min, snr_max), ylim=(0, radius_max), title=f"{name}: {valid.sum():,} measurements", xlabel="SourceXtractor++ S/N")
    ax.grid(alpha=0.15)
axes[0].set_ylabel(r"$R_{\mathrm{LOO}}$ [arcsec]")
fig.colorbar(sc, ax=axes, label="Wavelength [nm]", pad=0.02, shrink=0.8)
fig.savefig(ASSETS / "s1_s2_s3_snr_vs_Rloo.png", dpi=180, bbox_inches="tight")
plt.close(fig)
