#!/usr/bin/env python3
"""Animate CGCM SST anomaly output on a lon-lat map.

Example:
    python animate_cgcm_sst.py
    python animate_cgcm_sst.py --stride 3 --output cgcm_ssta.gif
    python animate_cgcm_sst.py --var ssta --frames 120 --output cgcm_ssta.mp4
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
        description="Animate SST/SSTA from a CGCM NetCDF output file."
    )
    parser.add_argument(
        "--input",
        type=Path,
        default=DEFAULT_INPUT,
        help=f"Input NetCDF file. Default: {DEFAULT_INPUT}",
    )
    parser.add_argument(
        "--var",
        default="ssta",
        help="Variable to plot. Default: ssta. If 'sst' is requested but absent, ssta is used.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Optional animation output path, e.g. cgcm_ssta.gif or cgcm_ssta.mp4.",
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
        "--vmin",
        type=float,
        default=-2.0,
        help="Lower color limit. Default: -3.",
    )
    parser.add_argument(
        "--vmax",
        type=float,
        default=2.0,
        help="Upper color limit. Default: 3.",
    )
    parser.add_argument(
        "--levels",
        type=int,
        default=25,
        help="Number of contour levels. Default: 25.",
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


def get_plot_variable(ds: xr.Dataset, requested: str) -> str:
    if requested in ds:
        return requested
    if requested.lower() == "sst" and "ssta" in ds:
        print("Variable 'sst' was not found; using 'ssta' from this output file.")
        return "ssta"
    available = ", ".join(ds.data_vars)
    raise KeyError(f"Variable {requested!r} not found. Available variables: {available}")


def prepare_field(ds: xr.Dataset, varname: str, stride: int, frames: int | None) -> xr.DataArray:
    field = ds[varname]
    if "time" not in field.dims:
        raise ValueError(f"Variable {varname!r} has no time dimension.")

    field = field.isel(time=slice(None, None, max(1, stride)))
    if frames is not None:
        field = field.isel(time=slice(0, frames))
    return field


def main() -> None:
    args = parse_args()
    input_path = args.input.expanduser().resolve()
    if not input_path.exists():
        raise FileNotFoundError(input_path)

    ds = xr.open_dataset(input_path)
    varname = get_plot_variable(ds, args.var)
    field = prepare_field(ds, varname, args.stride, args.frames)

    lon = ds["lon_p"].values
    lat = ds["lat_p"].values
    lon2d, lat2d = np.meshgrid(lon, lat)

    mask = ds.get("mask_sst", ds.get("mask_p"))
    mask_values = mask.values if mask is not None else np.ones((lat.size, lon.size))
    levels = np.linspace(args.vmin, args.vmax, args.levels)
    extent = [
        float(np.nanmin(lon)),
        float(np.nanmax(lon)),
        float(-30),
        float(30),
    ]

    use_cartopy = ccrs is not None and not args.no_cartopy
    if use_cartopy:
        data_crs = ccrs.PlateCarree()
        map_crs = ccrs.PlateCarree(central_longitude=180)
        fig, ax = plt.subplots(
            figsize=(9, 3),
            subplot_kw={"projection": map_crs},
        )
    else:
        data_crs = None
        fig, ax = plt.subplots(figsize=(9, 3))
    fig.subplots_adjust(left=0.08, right=0.88, bottom=0.14, top=0.9)

    if ccrs is None and not args.no_cartopy:
        print("Cartopy was not found; using plain Matplotlib axes.")

    def decorate_axis() -> None:
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
            ax.tick_params(labelsize=9)
            return

        ax.set_xlabel("Longitude")
        ax.set_ylabel("Latitude")
        ax.set_xlim(extent[0], extent[1])
        ax.set_ylim(extent[2], extent[3])
        ax.grid(True, linewidth=0.3, alpha=0.4)

    def draw(frame_index: int):
        ax.clear()
        values = field.isel(time=frame_index).values
        values = np.ma.masked_where(mask_values <= 0, values)

        contour_kwargs = {"transform": data_crs} if use_cartopy else {}
        ax.contourf(
            lon2d,
            lat2d,
            values,
            levels=levels,
            cmap="RdBu_r",
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
        decorate_axis()
        ax.set_title(f"{varname.upper()}  {np.datetime_as_string(field.time.values[frame_index], unit='D')}")
        return []

    contour_kwargs = {"transform": data_crs} if use_cartopy else {}
    first = ax.contourf(
        lon2d,
        lat2d,
        np.ma.masked_where(mask_values <= 0, field.isel(time=0).values),
        levels=levels,
        cmap="RdBu_r",
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
    decorate_axis()
    cbar = fig.colorbar(first, ax=ax, orientation="vertical", pad=0.02,shrink=0.7)
    cbar.set_label(f"{varname} (degC)")

    animation = FuncAnimation(
        fig,
        draw,
        frames=field.sizes["time"],
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
