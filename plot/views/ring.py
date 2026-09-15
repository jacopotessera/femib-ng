# "ring" view: velocity/pressure fields + a closed-polygon structure overlay
# + area/aspect-ratio time series.

import matplotlib
import matplotlib.colors
import matplotlib.pyplot as pyplot
import numpy
import typer

from .common import STRUCT_COLOR, overlay_structure, render_animation_gif


def _snapshot(data, path):
    x, y = data["X"][0], data["Y"][0]
    u, v = data["U"][0], data["V"][0]
    q = data["P"][0]
    sx, sy = data["SX"][0], data["SY"][0]
    has_struct = len(sx) > 0

    panels = []
    if len(x) and len(u):
        panels.append("velocity")
    if len(x) and len(q):
        panels.append("pressure")
    if not panels and len(x):
        panels.append("points")
    if not panels and not has_struct:
        typer.echo(f"warning: nothing to plot for {path}")
        return
    if not panels:
        panels.append("structure only")

    ncols = len(panels)
    fig, axes = pyplot.subplots(1, ncols, figsize=(5 * ncols, 4.6), squeeze=False)
    axes = axes[0]
    fig.set_tight_layout(True)

    for ax, kind in zip(axes, panels, strict=True):
        if kind == "velocity":
            mag = numpy.sqrt(u**2 + v**2)
            quiv = ax.quiver(
                x, y, u, v, mag, pivot="tail", units="xy", cmap=pyplot.cm.viridis
            )
            fig.colorbar(quiv, ax=ax, label="|u|")
            ax.set_title("velocity")
        elif kind == "pressure":
            try:
                cs = ax.tricontourf(x, y, q, levels=14, cmap=pyplot.cm.magma)
                fig.colorbar(cs, ax=ax, label="p")
            except (ValueError, RuntimeError):
                sc = ax.scatter(x, y, c=q, cmap=pyplot.cm.magma)
                fig.colorbar(sc, ax=ax, label="p")
            ax.set_title("pressure")
        elif kind == "points":
            ax.scatter(x, y)
            ax.set_title("points")
        else:
            ax.set_title("structure")
        if has_struct:
            overlay_structure(ax, sx, sy)
        # adjustable='box' keeps equal aspect by resizing the drawn box, not
        # by rescaling the data limits -- axis('equal') does the latter
        # (adjustable='datalim'), which would silently zoom the domain based
        # on whatever this frame's data extent happens to be.
        ax.set_aspect("equal", adjustable="box")

    title = path.name
    if has_struct and not numpy.isnan(data["AREA"][0]):
        title += "  (area={:.4g}, aspect={:.3f})".format(
            data["AREA"][0], data["ASPECT"][0]
        )
    fig.suptitle(title)

    fig.savefig(path, dpi=150)
    pyplot.close(fig)


def _animation(data, path):
    # Velocity/pressure aren't necessarily saved every step (e.g. only at
    # checkpoints); forward-fill each gap with the last checkpoint's field
    # so the quiver/pressure panels stay populated every frame instead of
    # flashing blank in between.
    last_x_positions = last_y_positions = numpy.array([])
    last_u_velocities = last_v_velocities = numpy.array([])
    last_pressures = numpy.array([])
    for i in range(len(data["T"])):
        if len(data["U"][i]):
            (
                last_x_positions,
                last_y_positions,
                last_u_velocities,
                last_v_velocities,
                last_pressures,
            ) = (
                data["X"][i],
                data["Y"][i],
                data["U"][i],
                data["V"][i],
                data["P"][i],
            )
        elif len(last_u_velocities):
            data["X"][i], data["Y"][i], data["U"][i], data["V"][i], data["P"][i] = (
                last_x_positions,
                last_y_positions,
                last_u_velocities,
                last_v_velocities,
                last_pressures,
            )

    has_u = any(len(u) for u in data["U"])
    has_q = any(len(q) for q in data["P"])
    has_struct = any(len(sx) for sx in data["SX"])

    top_panels = []
    if has_u:
        top_panels.append("velocity")
    if has_q:
        top_panels.append("pressure")
    if not top_panels and not has_struct:
        typer.echo(f"warning: nothing to plot for {path}")
        return

    nrows = 2 if has_struct else 1
    ncols = max(len(top_panels), 2 if has_struct else 1)
    # dpi=80 (matplotlib's default is ~100): this is read at a glance, not
    # zoomed into, and nearly every per-frame cost (drawing the
    # quiver/scatter, GIF quantization, GIF frame-diffing/encoding) scales
    # with pixel count.
    fig, axes = pyplot.subplots(
        nrows, ncols, figsize=(5 * ncols, 4.6 * nrows), squeeze=False, dpi=80
    )
    # NOT set_tight_layout(True): that installs a layout engine that
    # recomputes the tight bounding box of EVERY text artist (every tick
    # label, every axis, on every subplot) on every single draw -- with the
    # per-frame 'timestep N' xlabel update below changing text each frame,
    # that turned into ~120s of a ~280s render doing nothing but bbox math.
    # A single one-shot tight_layout() call once, right before rendering
    # starts, gets the same spacing without paying for it every frame.
    top_axes = axes[0]
    for ax in top_axes[len(top_panels) :]:
        ax.axis("off")

    # Fixed color scales across the whole animation (not per-frame) so
    # color is comparable frame-to-frame -- otherwise a quiet frame and a
    # vigorous frame would both autoscale to "full brightness", hiding
    # exactly the amplitude change a non-stationary run is meant to show.
    mags = [
        numpy.sqrt(u**2 + v**2)
        for u, v in zip(data["U"], data["V"], strict=True)
        if len(u)
    ]
    vmin, vmax = min((m.min() for m in mags), default=0.0), max(
        (m.max() for m in mags), default=1.0
    )
    if vmax <= vmin:
        vmax = vmin + 1e-9
    vcmap, vnorm = pyplot.cm.viridis, matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)

    # quiver's own scale=None default autoscales the arrow-to-data-unit
    # ratio from EACH CALL's own u,v -- which, called once per frame, makes
    # the same physical speed draw a different-length arrow in a quiet
    # frame than in a vigorous one. Fix scale globally instead, from the
    # same vmax used for color, so arrow length is comparable frame to
    # frame just like color already is: the fastest frame's arrow spans
    # about 3 grid cells (data["X"] values are a full min-to-max grid, so
    # the smallest positive gap between them is the grid spacing) -- long
    # enough to actually read as an arrow rather than a speck, at the cost
    # of some overlap between neighboring arrows on the fastest frames.
    # quiver_width is the shaft width, also in data ('xy') units since
    # units='xy' below applies to the whole arrow, not just its length;
    # quiver's own default width is tuned for the 'width'(axes-relative)
    # unit system, so left alone here it renders as a near-invisible hairline.
    xs = numpy.unique(numpy.concatenate([x for x in data["X"] if len(x)]))
    grid_spacing = numpy.min(numpy.diff(numpy.sort(xs))) if len(xs) > 1 else 0.04
    quiver_scale = vmax / (3 * grid_spacing) if vmax > 0 else 1.0
    quiver_width = 0.15 * grid_spacing

    all_p = [p for p in data["P"] if len(p)]
    pmin, pmax = min((p.min() for p in all_p), default=0.0), max(
        (p.max() for p in all_p), default=1.0
    )
    if pmax <= pmin:
        pmax = pmin + 1e-9
    pcmap, pnorm = pyplot.cm.magma, matplotlib.colors.Normalize(vmin=pmin, vmax=pmax)

    vel_ax = pres_ax = None
    idx = 0
    if has_u:
        vel_ax = top_axes[idx]
        idx += 1
        sm = pyplot.cm.ScalarMappable(cmap=vcmap, norm=vnorm)
        sm.set_array([])
        fig.colorbar(sm, ax=vel_ax, label="|u|")
    if has_q:
        pres_ax = top_axes[idx]
        idx += 1
        sm2 = pyplot.cm.ScalarMappable(cmap=pcmap, norm=pnorm)
        sm2.set_array([])
        fig.colorbar(sm2, ax=pres_ax, label="p")

    area_ax = aspect_ax = area_marker = aspect_marker = None
    if has_struct:
        area_ax, aspect_ax = axes[1][0], axes[1][1]
        for ax in axes[1][2:]:
            ax.axis("off")
        timesteps = numpy.array(data["T"])
        areas = numpy.array(data["AREA"], dtype=float)
        aspect_ratios = numpy.array(data["ASPECT"], dtype=float)
        area_ax.plot(timesteps, areas, color=STRUCT_COLOR)
        area_ax.set_title("structure area")
        area_ax.set_xlabel("timestep")
        (area_marker,) = area_ax.plot([], [], "o", color="black")
        aspect_ax.plot(timesteps, aspect_ratios, color=STRUCT_COLOR)
        aspect_ax.set_title("structure aspect ratio")
        aspect_ax.set_xlabel("timestep")
        (aspect_marker,) = aspect_ax.plot([], [], "o", color="black")

    # The sampling grid (x,y) is the same fixed regular grid on every
    # field-carrying frame (finite_element_space::plot always resamples the
    # same box at the same delta), so it only needs to be read once; each
    # frame then just pushes new U/V/color data into the SAME artists via
    # set_UVC/set_array/set_data, which is far cheaper than recreating them.
    field_idx = next((i for i in range(len(data["T"])) if len(data["X"][i])), None)
    x0, y0 = (
        (data["X"][field_idx], data["Y"][field_idx])
        if field_idx is not None
        else (numpy.array([]), numpy.array([]))
    )
    zeros0 = numpy.zeros_like(x0)

    vel_quiv = vel_struct_line = None
    if vel_ax is not None:
        vel_ax.set_title("velocity")
        vel_ax.set_xlim(0, 1)
        vel_ax.set_ylim(0, 1)
        # adjustable='box' keeps this exact xlim/ylim -- axis('equal')
        # (adjustable='datalim') would instead stretch the domain to fit
        # each frame's quiver-arrow extent, making it drift frame to frame
        # even though the fluid mesh itself never moves.
        vel_ax.set_aspect("equal", adjustable="box")
        if len(x0):
            vel_quiv = vel_ax.quiver(
                x0,
                y0,
                zeros0,
                zeros0,
                zeros0,
                pivot="tail",
                units="xy",
                cmap=vcmap,
                norm=vnorm,
                scale=quiver_scale,
                scale_units="xy",
                width=quiver_width,
            )
        (vel_struct_line,) = vel_ax.plot(
            [], [], color=STRUCT_COLOR, linewidth=2, zorder=5
        )

    pres_sc = pres_struct_line = None
    if pres_ax is not None:
        pres_ax.set_title("pressure")
        pres_ax.set_xlim(0, 1)
        pres_ax.set_ylim(0, 1)
        pres_ax.set_aspect("equal", adjustable="box")
        if len(x0):
            pres_sc = pres_ax.scatter(x0, y0, c=zeros0, cmap=pcmap, norm=pnorm, s=8)
        (pres_struct_line,) = pres_ax.plot(
            [], [], color=STRUCT_COLOR, linewidth=2, zorder=5
        )

    def update(i):
        u, v = data["U"][i], data["V"][i]
        q = data["P"][i]
        sx, sy = data["SX"][i], data["SY"][i]
        t = data["T"][i]
        xs_ = numpy.append(sx, sx[0]) if len(sx) else sx
        ys_ = numpy.append(sy, sy[0]) if len(sy) else sy

        if vel_quiv is not None and len(u):
            vel_quiv.set_UVC(u, v, numpy.sqrt(u**2 + v**2))
        if vel_struct_line is not None:
            vel_struct_line.set_data(xs_, ys_)
        if vel_ax is not None:
            # .xaxis.label.set_text(), not set_xlabel(): set_xlabel also
            # recomputes the label's on-axes position every call (to stay
            # clear of the tick labels below it) -- pure overhead here
            # since that position never actually needs to change frame to
            # frame, only the text content does.
            vel_ax.xaxis.label.set_text(f"timestep {t}")

        if pres_sc is not None and len(q):
            pres_sc.set_array(q)
        if pres_struct_line is not None:
            pres_struct_line.set_data(xs_, ys_)
        if pres_ax is not None:
            pres_ax.xaxis.label.set_text(f"timestep {t}")

        if area_ax is not None and not numpy.isnan(data["AREA"][i]):
            area_marker.set_data([t], [data["AREA"][i]])
            aspect_marker.set_data([t], [data["ASPECT"][i]])

    # Set the widest label text ("timestep <last>") before the one-shot
    # tight_layout() below, so spacing is computed for the actual worst case
    # up front rather than shifting slightly as the number of digits in the
    # timestep grows over the animation.
    last_t = data["T"][-1]
    if vel_ax is not None:
        vel_ax.set_xlabel(f"timestep {last_t}")
    if pres_ax is not None:
        pres_ax.set_xlabel(f"timestep {last_t}")
    fig.tight_layout()

    # (artist, axes) pairs, not bare artists: a Text label's own .axes
    # attribute isn't reliably set the way a plotted Line2D/PathCollection's
    # is, so draw_artist needs to be called via the axes we already know
    # each one belongs to.
    dynamic_artists = [
        (a, ax)
        for a, ax in (
            (vel_quiv, vel_ax),
            (vel_struct_line, vel_ax),
            (pres_sc, pres_ax),
            (pres_struct_line, pres_ax),
            (area_marker, area_ax),
            (aspect_marker, aspect_ax),
            (vel_ax.xaxis.label if vel_ax is not None else None, vel_ax),
            (pres_ax.xaxis.label if pres_ax is not None else None, pres_ax),
        )
        if a is not None
    ]

    render_animation_gif(fig, dynamic_artists, update, len(data["T"]), path)
    pyplot.close(fig)


def render(data, path):
    if len(data["T"]) == 1:
        _snapshot(data, path)
    else:
        _animation(data, path)
