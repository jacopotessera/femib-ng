# "ring-minimal" view: just the structure outline, nothing else. A second,
# deliberately sparse view to prove the --type override actually switches
# layouts, not a claim that this is the best minimal view.

import matplotlib.pyplot as pyplot
import numpy
import typer

from .common import STRUCT_COLOR, overlay_structure, render_animation_gif


def _snapshot(data, path):
    sx, sy = data["SX"][0], data["SY"][0]
    if len(sx) == 0:
        typer.echo(f"warning: nothing to plot for {path}")
        return

    fig, ax = pyplot.subplots(figsize=(5, 4.6))
    fig.set_tight_layout(True)
    overlay_structure(ax, sx, sy)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_aspect("equal", adjustable="box")

    title = path.name
    if not numpy.isnan(data["AREA"][0]):
        title += "  (area={:.4g}, aspect={:.3f})".format(
            data["AREA"][0], data["ASPECT"][0]
        )
    ax.set_title(title)

    fig.savefig(path, dpi=150)
    pyplot.close(fig)


def _animation(data, path, fps):
    has_struct = any(len(sx) for sx in data["SX"])
    if not has_struct:
        typer.echo(f"warning: nothing to plot for {path}")
        return

    fig, ax = pyplot.subplots(figsize=(5, 4.6), dpi=80)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_aspect("equal", adjustable="box")
    ax.set_title(path.name)
    (struct_line,) = ax.plot([], [], color=STRUCT_COLOR, linewidth=2, zorder=5)

    def update(i):
        sx, sy = data["SX"][i], data["SY"][i]
        xs_ = numpy.append(sx, sx[0]) if len(sx) else sx
        ys_ = numpy.append(sy, sy[0]) if len(sy) else sy
        struct_line.set_data(xs_, ys_)
        ax.xaxis.label.set_text("timestep {}".format(data["T"][i]))

    last_t = data["T"][-1]
    ax.set_xlabel(f"timestep {last_t}")
    fig.tight_layout()

    dynamic_artists = [(struct_line, ax), (ax.xaxis.label, ax)]
    render_animation_gif(fig, dynamic_artists, update, len(data["T"]), path, fps)
    pyplot.close(fig)


def render(data, path, fps=15):
    if len(data["T"]) == 1:
        _snapshot(data, path)
    else:
        _animation(data, path, fps)
