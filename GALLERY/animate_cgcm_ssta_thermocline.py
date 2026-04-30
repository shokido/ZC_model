#!/usr/bin/env python3
"""Animate CGCM SSTA and thermocline depth in a 2x1 Cartopy map.

Example:
    python animate_cgcm_ssta_thermocline.py
    python animate_cgcm_ssta_thermocline.py --stride 3 --output cgcm_ssta_h.gif
    python animate_cgcm_ssta_thermocline.py --frames 120 --output cgcm_ssta_h.mp4
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
import numpy as np
import xarray as xr

try:
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    from cartopy.mpl.ticker import LatitudeFormatter, LongitudeFormatter
except ImportError:
    ccrs = None
    cfeature = None
    LatitudeFormatter = None
    LongitudeFormatter = None


REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT = (
    REPO_ROOT
    / "OUTPUTS"
    / "CGCM"
    / "avg_cgcm_full_eqpac_30_ocn_clm_H120_cd1.4_dt3600_c10day.nc"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Animate SSTA and thermocline depth from a CGCM NetCDF output file."
    )
    parser.add_argument(
        "--input",
        type=Path,
        default=DEFAULT_INPUT,
        help=f"Input NetCDF file. Default: {DEFAULT_INPUT}",
    )
    parser.add_argument(
        "--ssta-var",
        default="ssta",
        help="SST anomaly variable. Default: ssta.",
    )
    parser.add_argument(
        "--thermocline-var",
        default="h_ocn_sw",
        help="Thermocline depth variable. Default: h_ocn_sw.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Optional animation output path, e.g. cgcm_ssta_thermocline.gif or .mp4.",
    )
    parser.add_argument(
        "--stride",
        type=int,
        default=1,
        help="Use every Nth time record. Default: 1.",
    )
    parser.add_argument(
        "--frames",
        type=int,
        default=None,
        help="Maximum number of frames to animate after applying stride.",
    )
    parser.add_argument(
        "--interval",
        type=int,
        default=120,
        help="Delay between frames in milliseconds. Default: 120.",
    )
    parser.add_argument(
        "--fps",
        type=int,
        default=8,
        help="Frames per second when saving. Default: 8.",
    )
    parser.add_argument(
        "--ssta-vmin",
        type=float,
        default=-2.0,
        help="Lower SSTA color limit. Default: -3.",
    )
    parser.add_argument(
        "--ssta-vmax",
        type=float,
        default=2.0,
        help="Upper SSTA color limit. Default: 3.",
    )
    parser.add_argument(
        "--thermocline-vmin",
        type=float,
        default=-20,
        help="Lower thermocline depth color limit. Default: robust symmetric value.",
    )
    parser.add_argument(
        "--thermocline-vmax",
        type=float,
        default=20,
        help="Upper thermocline depth color limit. Default: robust symmetric value.",
    )
    parser.add_argument(
        "--levels",
        type=int,
        default=25,
        help="Number of contour levels. Default: 25.",
    )
    parser.add_argument(
        "--lat-min",
        type=float,
        default=-25.0,
        help="Southern latitude limit. Default: -25.",
    )
    parser.add_argument(
        "--lat-max",
        type=float,
        default=25.0,
        help="Northern latitude limit. Default: 25.",
    )
    parser.add_argument(
        "--no-show",
        action="store_true",
        help="Do not open an interactive window after creating/saving the animation.",
    )
    parser.add_argument(
        "--no-cartopy",
        action="store_true",
        help="Use plain Matplotlib axes even when Cartopy is installed.",
    )
    return parser.parse_args()


def require_variable(ds: xr.Dataset, varname: str) -> xr.DataArray:
    if varname not in ds:
        available = ", ".join(ds.data_vars)
        raise KeyError(f"Variable {varname!r} not found. Available variables: {available}")
    field = ds[varname]
    if "time" not in field.dims:
        raise ValueError(f"Variable {varname!r} has no time dimension.")
    return field


def prepare_field(field: xr.DataArray, stride: int, frames: int | None) -> xr.DataArray:
    field = field.isel(time=slice(None, None, max(1, stride)))
    if frames is not None:
        field = field.isel(time=slice(0, frames))
    return field


def robust_symmetric_limits(values: np.ndarray, mask_values: np.ndarray) -> tuple[float, float]:
    mask3d = np.broadcast_to(mask_values <= 0, values.shape)
    sample = np.ma.masked_where(mask3d, values)
    abs_limit = float(np.nanpercentile(np.abs(sample.compressed()), 98))
    if abs_limit == 0.0:
        abs_limit = 1.0
    return -abs_limit, abs_limit


def main() -> None:
    args = parse_args()
    input_path = args.input.expanduser().resolve()
    if not input_path.exists():
        raise FileNotFoundError(input_path)

    ds = xr.open_dataset(input_path)
    ssta = prepare_field(require_variable(ds, args.ssta_var), args.stride, args.frames)
    thermocline = prepare_field(
        require_variable(ds, args.thermocline_var), args.stride, args.frames
    )

    lon = ds["lon_p"].values
    lat = ds["lat_p"].values
    lon2d, lat2d = np.meshgrid(lon, lat)
    extent = [
        float(np.nanmin(lon)),
        float(np.nanmax(lon)),
        -20,
        20,
    ]

    mask_sst = ds.get("mask_sst", ds.get("mask_p"))
    mask_p = ds.get("mask_p", ds.get("mask_sst"))
    if mask_sst is None or mask_p is None:
        raise KeyError("Neither mask_sst nor mask_p was found in the input file.")
    mask_sst_values = mask_sst.values
    mask_p_values = mask_p.values

    h_vmin = args.thermocline_vmin
    h_vmax = args.thermocline_vmax
    if h_vmin is None or h_vmax is None:
        h_vmin, h_vmax = robust_symmetric_limits(thermocline.values, mask_p_values)

    ssta_levels = np.linspace(args.ssta_vmin, args.ssta_vmax, args.levels)
    thermocline_levels = np.linspace(h_vmin, h_vmax, args.levels)

    use_cartopy = ccrs is not None and not args.no_cartopy
    if use_cartopy:
        data_crs = ccrs.PlateCarree()
        map_crs = ccrs.PlateCarree(central_longitude=180)
        fig, axes = plt.subplots(
            nrows=2,
            ncols=1,
            figsize=(8, 5),
            subplot_kw={"projection": map_crs},
        )
    else:
        data_crs = None
        fig, axes = plt.subplots(nrows=2, ncols=1, figsize=(8, 5))
    fig.subplots_adjust(left=0.08, right=0.86, bottom=0.08, top=0.93, hspace=0.24)

    if ccrs is None and not args.no_cartopy:
        print("Cartopy was not found; using plain Matplotlib axes.")

    def decorate_axis(ax) -> None:
        if use_cartopy:
            ax.set_extent(extent, crs=data_crs)
            ax.add_feature(
                cfeature.LAND,
                facecolor="0.82",
                edgecolor="0.45",
                linewidth=0.4,
                zorder=2,
            )
            ax.coastlines(resolution="110m", linewidth=0.6, color="0.25", zorder=3)
            ax.gridlines(
                draw_labels=False,
                linewidth=0.3,
                color="0.4",
                alpha=0.4,
                linestyle="-",
            )
            lon_step = 20.0
            lat_step = 10.0
            lon_ticks = np.arange(
                np.ceil(extent[0] / lon_step) * lon_step,
                extent[1] + lon_step,
                lon_step,
            )
            lat_ticks = np.arange(
                np.ceil(extent[2] / lat_step) * lat_step,
                extent[3] + lat_step,
                lat_step,
            )
            ax.set_xticks(lon_ticks, crs=data_crs)
            ax.set_yticks(lat_ticks, crs=data_crs)
            ax.xaxis.set_major_formatter(LongitudeFormatter())
            ax.yaxis.set_major_formatter(LatitudeFormatter())
            ax.tick_params(labelsize=8)
            return

        ax.set_xlabel("Longitude")
        ax.set_ylabel("Latitude")
        ax.set_xlim(extent[0], extent[1])
        ax.set_ylim(extent[2], extent[3])
        ax.grid(True, linewidth=0.3, alpha=0.4)

    def masked_frame(field: xr.DataArray, mask_values: np.ndarray, frame_index: int):
        values = field.isel(time=frame_index).values
        return np.ma.masked_where(mask_values <= 0, values)

    def draw_panel(ax, values, mask_values, levels, cmap, title):
        contour_kwargs = {"transform": data_crs} if use_cartopy else {}
        contour = ax.contourf(
            lon2d,
            lat2d,
            values,
            levels=levels,
            cmap=cmap,
            extend="both",
            zorder=1,
            **contour_kwargs,
        )
        ax.contour(
            lon2d,
            lat2d,
            mask_values,
            levels=[0.5],
            colors="0.25",
            linewidths=0.5,
            zorder=4,
            **contour_kwargs,
        )
        decorate_axis(ax)
        ax.set_title(title, fontsize=11)
        return contour

    def draw(frame_index: int):
        for ax in axes:
            ax.clear()

        date_label = np.datetime_as_string(ssta.time.values[frame_index], unit="D")
        draw_panel(
            axes[0],
            masked_frame(ssta, mask_sst_values, frame_index),
            mask_sst_values,
            ssta_levels,
            "RdBu_r",
            f"SST anomaly ({args.ssta_var})  {date_label}",
        )
        draw_panel(
            axes[1],
            masked_frame(thermocline, mask_p_values, frame_index),
            mask_p_values,
            thermocline_levels,
            "RdBu_r",
            f"Thermocline depth ({args.thermocline_var})  {date_label}",
        )
        return []

    first_ssta = draw_panel(
        axes[0],
        masked_frame(ssta, mask_sst_values, 0),
        mask_sst_values,
        ssta_levels,
        "RdBu_r",
        f"SST anomaly ({args.ssta_var})",
    )
    first_h = draw_panel(
        axes[1],
        masked_frame(thermocline, mask_p_values, 0),
        mask_p_values,
        thermocline_levels,
        "RdBu_r",
        f"Thermocline depth ({args.thermocline_var})",
    )

    cbar_ssta = fig.colorbar(first_ssta, ax=axes[0], orientation="vertical", pad=0.02,shrink=0.7)
    cbar_ssta.set_label("SSTA (degC)")
    cbar_h = fig.colorbar(first_h, ax=axes[1], orientation="vertical", pad=0.02,shrink=0.7)
    cbar_h.set_label("Thermocline depth (m)")

    animation = FuncAnimation(
        fig,
        draw,
        frames=ssta.sizes["time"],
        interval=args.interval,
        blit=False,
    )

    if args.output is not None:
        output_path = args.output.expanduser()
        if not output_path.is_absolute():
            output_path = Path(__file__).resolve().parent / output_path
        output_path.parent.mkdir(parents=True, exist_ok=True)
        if output_path.suffix.lower() == ".gif":
            animation.save(output_path, writer=PillowWriter(fps=args.fps))
        else:
            animation.save(output_path, fps=args.fps)
        print(f"Saved animation to {output_path}")

    if not args.no_show:
        plt.show()


if __name__ == "__main__":
    main()
