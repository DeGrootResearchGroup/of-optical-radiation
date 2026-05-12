#!/usr/bin/env python3
"""Generate result figures for the Chiu 1999 pilot UV channel case.

Produces, in postProcessing/figures/:

  01_dose_distribution.png  histogram of per-particle dose vs the
                            bimodal shape in Chiu Fig. 8b
  02_G_contour.png          fluence-rate G from the DOM solve
                            (analogue of Chiu Fig. 4's PSS contour)
  03_U_streamlines.png      mean velocity magnitude + streamlines
                            (analogue of Chiu Fig. 5's LDV vectors)
  04_trajectories.png       sample of the radiationDose Lagrangian
                            tracks coloured by accumulated dose
                            (analogue of Chiu Fig. 7)

Dependencies:

    pip install matplotlib numpy scipy vtk

Inputs:

    postProcessing/radiationDose/<time>/doseDistribution.csv
    postProcessing/radiationDose/<time>/trajectories.vtk
    VTK/uvChannelChiu1999_<flow_time>.vtk
    VTK/uvChannelChiu1999_<dom_time>.vtk

The VTK/ files come from the OpenFOAM utility:

    foamToVTK -time "<flow_time>,<dom_time>" -fields "(U G k epsilon)"

Allrun-DOM does not call foamToVTK automatically; run it once after the
case finishes and then run this script.
"""

import csv
import glob
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import vtk
from matplotlib.collections import LineCollection
from matplotlib.colors import LogNorm
from scipy.interpolate import griddata
from vtk.util.numpy_support import vtk_to_numpy


# -- Geometry constants (kept in sync with make_mesh.py) ------------------

ROW_PITCH_X  = 0.125
LAMP_PITCH_Y = 0.075
LAMP_RADIUS  = 0.0125

N_ROWS        = 5
LAMPS_PER_ROW = 4


def lamp_centres():
    x_lamps = np.arange(N_ROWS) * ROW_PITCH_X
    y_centred = (np.arange(LAMPS_PER_ROW) - (LAMPS_PER_ROW - 1) / 2.0) * LAMP_PITCH_Y
    y_even = y_centred
    y_odd  = y_centred + LAMP_PITCH_Y / 2.0
    out = []
    for ix, xc in enumerate(x_lamps):
        for yc in (y_odd if ix % 2 == 1 else y_even):
            out.append((xc, yc))
    return out


def draw_lamps(ax, color='gray', edgecolor='black', lw=0.6):
    for xc, yc in lamp_centres():
        ax.add_patch(
            plt.Circle((xc * 100, yc * 100), LAMP_RADIUS * 100,
                       facecolor=color, edgecolor=edgecolor, lw=lw, zorder=10)
        )


# -- VTK helpers ----------------------------------------------------------

def read_unstructured_vtk(path):
    reader = vtk.vtkUnstructuredGridReader()
    reader.SetFileName(path)
    reader.ReadAllScalarsOn()
    reader.ReadAllVectorsOn()
    reader.Update()
    grid = reader.GetOutput()
    cell_data = grid.GetCellData()
    fields = {
        cell_data.GetArrayName(i): vtk_to_numpy(cell_data.GetArray(i))
        for i in range(cell_data.GetNumberOfArrays())
    }
    centres_filter = vtk.vtkCellCenters()
    centres_filter.SetInputData(grid)
    centres_filter.Update()
    centres = vtk_to_numpy(centres_filter.GetOutput().GetPoints().GetData())
    return centres, fields


def read_polydata_trajectories(path):
    reader = vtk.vtkPolyDataReader()
    reader.SetFileName(path)
    reader.ReadAllScalarsOn()
    reader.Update()
    polydata = reader.GetOutput()
    pts = vtk_to_numpy(polydata.GetPoints().GetData())
    dose = vtk_to_numpy(polydata.GetPointData().GetArray("dose_mJcm2"))
    lines = polydata.GetLines()
    lines.InitTraversal()
    trajs = []
    id_list = vtk.vtkIdList()
    n_traj = polydata.GetNumberOfLines()
    for i in range(n_traj):
        lines.GetNextCell(id_list)
        ids = [id_list.GetId(j) for j in range(id_list.GetNumberOfIds())]
        trajs.append((pts[ids], dose[ids]))
    return trajs


# -- Figure builders ------------------------------------------------------

def figure_dose_distribution(case_dir, out_path):
    csv_paths = sorted(glob.glob(
        os.path.join(case_dir, "postProcessing", "radiationDose",
                     "*", "doseDistribution.csv")))
    if not csv_paths:
        return None
    csv_path = csv_paths[-1]
    doses, escaped = [], 0
    with open(csv_path) as f:
        for row in csv.DictReader(f):
            doses.append(float(row["dose_mJ_cm2"]))
            if row["endReason"].strip() == "escaped":
                escaped += 1
    doses = np.array(doses)
    mean_dose = doses.mean()

    fig, ax = plt.subplots(figsize=(7, 4))
    ax.hist(doses, bins=60, range=(0, 60),
            color='steelblue', edgecolor='black', linewidth=0.5)
    ax.set_xlabel("Dose, $D$ (mJ/cm$^2$)")
    ax.set_ylabel("Number of particles per bin")
    ax.set_title(
        f"Chiu 1999 pilot UV channel, V = 24 cm/s\n"
        f"n = {len(doses)}, escape = {100 * escaped / len(doses):.1f} %, "
        f"$\\overline{{D}}$ = {mean_dose:.1f} mJ/cm$^2$"
    )
    ax.axvline(mean_dose, color='r', ls='--', lw=1,
               label=f"mean = {mean_dose:.1f}")
    ax.legend()
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    return out_path


def figure_G_contour(case_dir, out_path):
    vtk_paths = sorted(glob.glob(
        os.path.join(case_dir, "VTK", "uvChannelChiu1999_*.vtk")),
        key=lambda p: int(os.path.basename(p).split("_")[1].split(".")[0]),
    )
    if not vtk_paths:
        return None
    # Latest-time VTK has G; use it.
    centres, fields = read_unstructured_vtk(vtk_paths[-1])
    if "G" not in fields:
        return None
    G = fields["G"]
    xi = np.linspace(-0.1, 0.6, 700)
    yi = np.linspace(-0.165, 0.165, 200)
    XI, YI = np.meshgrid(xi, yi)
    GI = griddata((centres[:, 0], centres[:, 1]), G, (XI, YI),
                  method='linear', fill_value=np.nan)

    fig, ax = plt.subplots(figsize=(11, 4))
    levels = np.linspace(0, max(40, float(np.nanpercentile(GI, 99))), 21)
    cf = ax.contourf(XI * 100, YI * 100, GI, levels=levels,
                     cmap='inferno', extend='max')
    draw_lamps(ax, color='gray', edgecolor='white', lw=0.6)
    ax.set_xlabel("Streamwise x (cm)")
    ax.set_ylabel("Transverse y (cm)")
    ax.set_aspect('equal')
    ax.set_xlim(-10, 60)
    ax.set_ylim(-16.5, 16.5)
    ax.set_title("Fluence rate $G$ (W/m$^2$) from DOM")
    fig.colorbar(cf, ax=ax, label="$G$ (W/m$^2$)", shrink=0.85)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    return out_path


def figure_U_streamlines(case_dir, out_path):
    vtk_paths = sorted(glob.glob(
        os.path.join(case_dir, "VTK", "uvChannelChiu1999_*.vtk")),
        key=lambda p: int(os.path.basename(p).split("_")[1].split(".")[0]),
    )
    # Pick the earliest time that has U (the flow snapshot, before the
    # DOM step bumped the time forward).
    centres, fields = None, None
    for path in vtk_paths:
        c, f = read_unstructured_vtk(path)
        if "U" in f:
            centres, fields = c, f
            break
    if fields is None:
        return None
    U = fields["U"]

    xi = np.linspace(-0.1, 0.6, 700)
    yi = np.linspace(-0.165, 0.165, 200)
    XI, YI = np.meshgrid(xi, yi)
    UI = griddata((centres[:, 0], centres[:, 1]), U[:, 0], (XI, YI),
                  method='linear', fill_value=np.nan)
    VI = griddata((centres[:, 0], centres[:, 1]), U[:, 1], (XI, YI),
                  method='linear', fill_value=np.nan)
    Umag = np.sqrt(UI**2 + VI**2)

    fig, ax = plt.subplots(figsize=(11, 4))
    cf = ax.contourf(XI * 100, YI * 100, Umag, levels=21, cmap='viridis')
    ax.streamplot(XI * 100, YI * 100, UI, VI, density=2.5,
                  color='white', linewidth=0.5, arrowsize=0.6)
    draw_lamps(ax, color='gray', edgecolor='black', lw=0.6)
    ax.set_xlabel("Streamwise x (cm)")
    ax.set_ylabel("Transverse y (cm)")
    ax.set_aspect('equal')
    ax.set_xlim(-10, 60)
    ax.set_ylim(-16.5, 16.5)
    ax.set_title("Mean velocity $|U|$ (m/s) and streamlines")
    fig.colorbar(cf, ax=ax, label="$|U|$ (m/s)", shrink=0.85)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    return out_path


def figure_trajectories(case_dir, out_path, sample_size=150, seed=0):
    vtk_paths = sorted(glob.glob(
        os.path.join(case_dir, "postProcessing", "radiationDose",
                     "*", "trajectories.vtk")))
    if not vtk_paths:
        return None
    trajs = read_polydata_trajectories(vtk_paths[-1])
    rng = np.random.default_rng(seed)
    indices = rng.choice(len(trajs), size=min(sample_size, len(trajs)),
                         replace=False)
    sample = [trajs[i] for i in indices]
    all_doses = np.concatenate([d for _, d in sample])
    pos = all_doses[all_doses > 0]
    vmin = float(pos.min()) if pos.size else 0.1
    vmax = float(all_doses.max())

    fig, ax = plt.subplots(figsize=(11, 4))
    for pts, dose in sample:
        segs = np.stack([pts[:-1, :2] * 100, pts[1:, :2] * 100], axis=1)
        mid_dose = (dose[:-1] + dose[1:]) / 2 + 1e-3
        lc = LineCollection(segs, array=mid_dose, cmap='plasma',
                            norm=LogNorm(vmin=vmin, vmax=vmax),
                            linewidth=0.7, alpha=0.7)
        ax.add_collection(lc)
    draw_lamps(ax, color='gray', edgecolor='black', lw=0.5)
    ax.set_xlabel("Streamwise x (cm)")
    ax.set_ylabel("Transverse y (cm)")
    ax.set_aspect('equal')
    ax.set_xlim(-50, 100)
    ax.set_ylim(-16.5, 16.5)
    ax.set_title(f"Sample particle trajectories ({len(sample)} of "
                 f"{len(trajs)}), coloured by accumulated dose")
    fig.colorbar(lc, ax=ax, label="Dose (mJ/cm$^2$)", shrink=0.85)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    return out_path


def figure_convergence_dose(case_dir, out_path):
    """Overlay dose distributions from a mesh-convergence sweep.

    Reads convergence/cell_*mm/postProcessing/radiationDose/<time>/
    doseDistribution.csv (each subcase is a self-contained OpenFOAM
    case spawned by Allrun-convergence) and plots a smoothed histogram
    per mesh.
    """
    sweep_root = os.path.join(case_dir, "convergence")
    if not os.path.isdir(sweep_root):
        return None

    runs = []
    for d in sorted(os.listdir(sweep_root)):
        sub_dir = os.path.join(sweep_root, d)
        if not os.path.isdir(sub_dir):
            continue
        match = d.replace("cell_", "").replace("mm", "")
        try:
            size_mm = float(match)
        except ValueError:
            continue
        # Find the latest dose CSV in this subcase's postProcessing dir.
        csv_paths = sorted(glob.glob(os.path.join(
            sub_dir, "postProcessing", "radiationDose", "*",
            "doseDistribution.csv",
        )))
        if not csv_paths:
            continue
        doses = []
        escaped = 0
        with open(csv_paths[-1]) as f:
            for row in csv.DictReader(f):
                doses.append(float(row["dose_mJ_cm2"]))
                if row["endReason"].strip() == "escaped":
                    escaped += 1
        runs.append((size_mm, np.array(doses), escaped))

    if len(runs) < 2:
        return None

    # Sort coarse -> fine
    runs.sort(key=lambda r: -r[0])

    fig, (ax, ax_tab) = plt.subplots(
        1, 2, figsize=(12, 4.2),
        gridspec_kw={"width_ratios": [3, 1.4]},
    )

    bins = np.linspace(0, 60, 61)
    centres = 0.5 * (bins[:-1] + bins[1:])
    colours = plt.cm.viridis(np.linspace(0.15, 0.85, len(runs)))
    for (size_mm, doses, esc), c in zip(runs, colours):
        h, _ = np.histogram(doses, bins=bins)
        ax.plot(centres, h, color=c, lw=1.5,
                label=f"{size_mm:g} mm  (n={len(doses)}, "
                      f"$\\overline{{D}}$={doses.mean():.1f})")
    ax.set_xlabel("Dose, $D$ (mJ/cm$^2$)")
    ax.set_ylabel("Particles per 1 mJ/cm$^2$ bin")
    ax.set_title("Mesh-convergence sweep at V = 24 cm/s")
    ax.grid(alpha=0.3)
    ax.legend(title="cell size")

    # Summary table on the right axes
    ax_tab.axis("off")
    rows = [("cell (mm)", "mean D", "log(N=1)", "log(N=3)", "esc %")]
    for size_mm, doses, esc in runs:
        # Local series-event log reduction
        def log_red(k, n):
            x = k * doses
            x = np.clip(x, 0, 700)
            terms = np.zeros_like(x)
            fact = 1.0
            xi = np.ones_like(x)
            for i in range(n):
                terms += xi / fact
                xi *= x
                fact *= (i + 1)
            survival = (np.exp(-x) * terms).mean()
            return -np.log10(max(survival, 1e-30))
        rows.append((
            f"{size_mm:g}",
            f"{doses.mean():.2f}",
            f"{log_red(0.5887, 1):.2f}",
            f"{log_red(1.107, 3):.2f}",
            f"{100*esc/len(doses):.1f}",
        ))
    tbl = ax_tab.table(cellText=rows[1:], colLabels=rows[0],
                       loc='center', cellLoc='center')
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(9)
    tbl.scale(1.0, 1.4)
    for k in range(len(rows[0])):
        tbl[(0, k)].set_facecolor('#e0e0e0')

    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    return out_path


def main():
    case_dir = os.path.dirname(os.path.abspath(__file__))
    out_dir = os.path.join(case_dir, "postProcessing", "figures")
    os.makedirs(out_dir, exist_ok=True)

    builders = [
        ("01_dose_distribution.png", figure_dose_distribution),
        ("02_G_contour.png",         figure_G_contour),
        ("03_U_streamlines.png",     figure_U_streamlines),
        ("04_trajectories.png",      figure_trajectories),
        ("05_convergence_dose.png",  figure_convergence_dose),
    ]

    wrote = []
    for fname, builder in builders:
        out_path = os.path.join(out_dir, fname)
        try:
            result = builder(case_dir, out_path)
        except Exception as e:
            sys.stderr.write(f"plot_results: {fname}: {type(e).__name__}: {e}\n")
            continue
        if result is None:
            sys.stderr.write(
                f"plot_results: {fname}: required inputs not found; "
                f"have you run Allrun-DOM and foamToVTK?\n"
            )
            continue
        wrote.append(result)

    if not wrote:
        sys.stderr.write(
            "plot_results: no figures produced. Run, in order:\n"
            "  ./Allrun-DOM\n"
            "  foamToVTK -time '<flow_time>,<dom_time>' "
            "-fields '(U G k epsilon)'\n"
            "  python3 plot_results.py\n"
        )
        return 1

    for p in wrote:
        print(p)
    return 0


if __name__ == "__main__":
    sys.exit(main())
