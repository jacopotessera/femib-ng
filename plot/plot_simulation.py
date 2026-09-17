#!/bin/python

from pathlib import Path

import h5py
import matplotlib
import typer

from plot_utils import calc_plot_data, sorted_timesteps
from views import ring, ring_minimal

# Whatever GUI backend matplotlib picked by default (Qt/Tk/...) routes every
# frame through its toolkit's event loop and widget machinery for no reason
# here -- this script always renders to a file, never a window.
matplotlib.use("Agg", force=True)

VIEWS = {
    "ring": ring.render,
    "ring-minimal": ring_minimal.render,
}


def main(
    input_path: Path = typer.Argument(
        ..., help="Simulation .h5 file (from femib::write::save_sim/save_plot_data)."
    ),
    output_path: Path | None = typer.Argument(
        None,
        help="Where to write the render. Defaults to <input_path> with a "
        ".png (single timestep) or .gif (multiple timesteps) suffix.",
    ),
    view: str | None = typer.Option(
        None,
        "--type",
        help="Which registered view to render with. Defaults to the file's "
        "own sim_type attribute (written by save_sim), falling back "
        "to 'ring' if that's unset.",
    ),
    fps: int = typer.Option(
        15, "--fps", min=1, help="Frames per second for animated (.gif) output."
    ),
    stride: int = typer.Option(
        1,
        "--stride",
        min=1,
        help="Render every Nth simulation timestep. GIF frame delays can't go "
        "below 10ms (100fps), so past that ceiling this is the only way to "
        "make the animation play faster.",
    ),
):
    """Render a femib-ng simulation .h5 file to a .png or .gif."""
    if not input_path.exists():
        typer.echo(f"error: {input_path} does not exist")
        raise typer.Exit(code=1)

    with h5py.File(input_path, "r") as f:
        groups = sorted_timesteps(f)
        if not groups:
            typer.echo(f"error: no timestep data in {input_path}")
            raise typer.Exit(code=1)
        data = calc_plot_data(groups)
        stored_view = f.attrs.get("sim_type", "")

    if stride > 1:
        data = {k: v[::stride] for k, v in data.items()}

    view_name = view or stored_view or "ring"
    if view_name not in VIEWS:
        available = ", ".join(sorted(VIEWS))
        typer.echo(f"error: unknown plot type {view_name!r} (available: {available})")
        raise typer.Exit(code=1)

    if output_path is None:
        suffix = ".gif" if len(data["T"]) > 1 else ".png"
        output_path = input_path.with_suffix(suffix)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    typer.echo(f"Rendering {input_path} as '{view_name}' -> {output_path}")
    VIEWS[view_name](data, output_path, fps)
    typer.echo(f"Wrote {output_path}")


if __name__ == "__main__":
    typer.run(main)
