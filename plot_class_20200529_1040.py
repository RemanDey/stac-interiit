
import glob
import os
import re
from datetime import datetime, timezone

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import pandas as pd
from astropy.io import fits

GAIN_EV = 13.5
BASE = "/home/remandey/my-programs/isro-interiit/ch2_cla_l1_2020_05/cla/data/calibrated/2020/05/29"
START_BOUND = datetime(2020, 5, 29, 10, 40, 0, tzinfo=timezone.utc)
END_BOUND = datetime(2020, 5, 29, 10, 43, 0, tzinfo=timezone.utc)
OUTDIR = "/home/remandey/my-programs/isro-interiit/plots_20200529_1040"

FNAME_RE = re.compile(r"ch2_cla_l1_(\d{8}T\d{6})(\d{3})_(\d{8}T\d{6})(\d{3})\.fits")

LINES = {  # keV : label
    1.25: "Mg-Ka", 1.49: "Al-Ka\n(inst.+XRF)", 1.74: "Si-Ka",
    2.31: "S-Ka", 3.69: "Ca-Ka", 4.51: "Ti-Ka",
    6.40: "Fe-Ka", 8.04: "Cu-Ka (inst.)",
}


def fname_to_times(path):
    m = FNAME_RE.search(os.path.basename(path))
    s = datetime.strptime(m.group(1), "%Y%m%dT%H%M%S").replace(tzinfo=timezone.utc)
    e = datetime.strptime(m.group(3), "%Y%m%dT%H%M%S").replace(tzinfo=timezone.utc)
    return s, e


def list_interval_files(base=BASE, start=START_BOUND, end=END_BOUND):
    files = sorted(glob.glob(os.path.join(base, "ch2_cla_l1_*.fits")))
    sel = [f for f in files if (lambda se: se[0] < end and se[1] > start)(fname_to_times(f))]
    return sel


def read_class_l1(path):
    with fits.open(path) as hdul:
        data = hdul[1].data
        hdr = hdul[1].header
    ch = np.asarray(data["CHANNEL"], dtype=int)
    counts = np.asarray(data["COUNTS"], dtype=float)
    energy_kev = ch * GAIN_EV / 1000.0
    meta = {k: hdr.get(k) for k in
            ["STARTIME", "ENDTIME", "MID_UTC", "EXPOSURE", "GAIN", "SCD_USED",
             "SAT_LAT", "SAT_LON", "BORE_LAT", "BORE_LON", "SAT_ALT",
             "SOLARANG", "PHASEANG", "TEMP", "LST_HR", "LST_MIN", "LST_SEC"]}
    return ch, counts, energy_kev, meta


def band_counts(counts, emin, emax, energy):
    m = (energy >= emin) & (energy <= emax)
    return counts[m].sum()


def main():
    os.makedirs(OUTDIR, exist_ok=True)
    files = list_interval_files()
    print(f"Selected {len(files)} files overlapping 10:40-10:43 UTC")
    assert len(files) >= 20, f"expected ~22 files, got {len(files)}"

    rows, spectra = [], []
    for f in files:
        ch, counts, en, meta = read_class_l1(f)
        spectra.append(counts)
        s, e = fname_to_times(f)
        mid = datetime.strptime(meta["MID_UTC"], "%Y-%m-%dT%H:%M:%S.%f").replace(tzinfo=timezone.utc) \
            if meta["MID_UTC"] else s + (e - s) / 2
        rows.append({
            "file": os.path.basename(f), "t_start": s, "t_end": e, "t_mid": mid,
            "total_0p5_10": band_counts(counts, 0.5, 10.0, en),
            "soft_0p5_2": band_counts(counts, 0.5, 2.0, en),
            "hard_7_15": band_counts(counts, 7.0, 15.0, en),
            **{k: meta[k] for k in ["SAT_LAT", "SAT_LON", "BORE_LAT", "BORE_LON",
                                    "SAT_ALT", "SOLARANG", "PHASEANG", "TEMP", "EXPOSURE"]},
        })
    df = pd.DataFrame(rows).sort_values("t_mid")
    cube = np.vstack(spectra)
    summed = cube.sum(axis=0)
    ch = np.arange(2048)
    en = ch * GAIN_EV / 1000.0
    exposure = float(df["EXPOSURE"].sum())
    print(f"Summed exposure = {exposure:.0f} s over {len(df)} x 8 s spectra")
    print(f"Lat {df['BORE_LAT'].min():.2f}..{df['BORE_LAT'].max():.2f}, "
          f"Lon {df['BORE_LON'].min():.2f}..{df['BORE_LON'].max():.2f}")
    df.to_csv(os.path.join(OUTDIR, "geometry_table.csv"), index=False)

    # ---- P0: single 8 s example ----
    ch0, c0, en0, m0 = read_class_l1(files[len(files) // 2])
    fig, ax = plt.subplots(figsize=(9, 4.5))
    ax.step(en0, np.where(c0 <= 0, 0.5, c0), where="mid")
    ax.set(xlim=(0.5, 10), yscale="log", xlabel="Energy (keV)", ylabel="Counts / 8 s / ch",
           title=f"CLASS L1 single 8 s example\n{os.path.basename(files[len(files)//2])}")
    fig.tight_layout()
    fig.savefig(os.path.join(OUTDIR, "P0_example_8s.png"), dpi=150)
    plt.close(fig)

    # ---- P1: integrated 3-min spectrum ----
    sci = (ch >= 37) & (ch <= 800)
    err = np.sqrt(np.where(summed > 0, summed, 1.0))
    fig, ax = plt.subplots(figsize=(10, 5.2))
    ax.errorbar(en[sci], np.where(summed[sci] <= 0, 0.5, summed[sci]), yerr=err[sci],
                fmt=".", ms=3, ecolor="gray", elinewidth=0.6, capsize=0, label="summed counts")
    ax.set(xlim=(0.5, 10), yscale="log", xlabel="Energy (keV)",
           ylabel=f"Counts per 13.5 eV channel / {exposure:.0f} s",
           title=f"CLASS integrated spectrum 2020-05-29 10:40-10:43 UTC ({len(df)}x8 s, {exposure:.0f} s)")
    ymax = np.max(summed[sci]) * 2.2
    for kev, lab in LINES.items():
        ax.axvline(kev, color="r", ls="--", lw=0.8, alpha=0.7)
        ax.text(kev, ymax, f" {lab}\n {kev:.2f} keV", rotation=0, fontsize=7,
                va="top", ha="left", color="darkred")
    ax.set_ylim(0.8, ymax * 1.6)
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3, which="both")
    fig.tight_layout()
    fig.savefig(os.path.join(OUTDIR, "P1_integrated_spectrum.png"), dpi=150)
    plt.close(fig)

    # ---- P2: lightcurves ----
    fig, axs = plt.subplots(3, 1, figsize=(10, 7), sharex=True)
    t = pd.to_datetime(df["t_mid"])
    for ax, col, ttl in zip(axs,
                            ["total_0p5_10", "soft_0p5_2", "hard_7_15"],
                            ["Total 0.5-10 keV (counts/8 s)", "Soft 0.5-2 keV (counts/8 s)",
                             "Hard 7-15 keV particle monitor (counts/8 s)"]):
        ax.step(t, df[col], where="mid")
        ax.set_ylabel(ttl, fontsize=9)
        ax.grid(alpha=0.3)
    axs[-1].set_xlabel("UTC 2020-05-29")
    axs[-1].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M:%S"))
    fig.suptitle("CLASS L1 count-rate vs time (per-8 s files, 10:40-10:43 UTC)")
    fig.tight_layout()
    fig.savefig(os.path.join(OUTDIR, "P2_lightcurves.png"), dpi=150)
    plt.close(fig)
    fig, axs = plt.subplots(3, 1, figsize=(10, 7), sharex=True)
    axs[0].plot(t, df["BORE_LAT"], "o-", ms=3, label="Boresight")
    axs[0].plot(t, df["SAT_LAT"], "s-", ms=2, alpha=0.6, label="Sub-satellite")
    axs[0].set_ylabel("Latitude (deg)")
    axs[0].legend(fontsize=8)
    axs[1].plot(t, df["BORE_LON"], "o-", ms=3)
    axs[1].set_ylabel("Boresight lon (deg)")
    axs[2].plot(t, df["SOLARANG"], "o-", ms=3, label="Solar ang")
    axs[2].plot(t, df["PHASEANG"], "s-", ms=2, alpha=0.6, label="Phase ang")
    axs[2].set_ylabel("Angle (deg)")
    axs[2].legend(fontsize=8)
    axs[2].set_xlabel("UTC 2020-05-29")
    for ax in axs:
        ax.grid(alpha=0.3)
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M:%S"))
    altr = f"Alt {df['SAT_ALT'].mean():.1f} km | Temp {df['TEMP'].mean():.1f} C | SCDs: {read_class_l1(files[0])[3]['SCD_USED']}"
    fig.suptitle(f"CLASS ground track 10:40-10:43 UTC\n{altr}")
    fig.tight_layout()
    fig.savefig(os.path.join(OUTDIR, "P3_geometry.png"), dpi=150)
    plt.close(fig)

    print("Wrote:", sorted(os.listdir(OUTDIR)))


if __name__ == "__main__":
    main()
