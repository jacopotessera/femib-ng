import os
import shutil
import subprocess
import tempfile

import numpy
from PIL import Image

# Structure color is deliberately distinct from both the velocity (viridis)
# and pressure (magma) colormaps so the material curve reads clearly against
# either field.
STRUCT_COLOR = "#d9750a"


def overlay_structure(ax, sx, sy):
    if len(sx) == 0:
        return
    xs = numpy.append(sx, sx[0])
    ys = numpy.append(sy, sy[0])
    ax.plot(xs, ys, color=STRUCT_COLOR, linewidth=2, zorder=5)


def render_animation_gif(fig, dynamic_artists, update_fn, n_frames, path, fps=15):
    # Manual blitting instead of matplotlib's Animation.save(): that method
    # always does a FULL canvas redraw per output frame (recomputing every
    # tick label's text layout, every axis, every colorbar) regardless of
    # any blit setting -- blit only ever applies to the on-screen path, never
    # to file output. Measured via cProfile at ~120-150s of a ~170s render on
    # the ring animation, almost all of it in text/bbox layout code, not
    # actual pixel drawing. None of the static decoration (axes, ticks,
    # colorbars, any background line plots) changes frame to frame, so it's
    # rendered once here, cached as a bitmap, and each frame only draws the
    # handful of artists that actually change on top of that background --
    # the standard matplotlib blitting pattern, driven by hand since
    # Animation.save() doesn't use it.
    if shutil.which("ffmpeg") is None:
        raise RuntimeError(
            "render_animation_gif requires ffmpeg on PATH: frames are "
            "encoded by piping a PNG sequence through it instead of "
            "buffering every frame in Python -- Pillow's own GIF writer "
            "(and imageio's, which just wraps it) builds a second full "
            "in-memory copy of every frame before writing any of them, "
            "which OOMs on long animations."
        )

    for a, _ in dynamic_artists:
        a.set_animated(True)
        # set_animated(True) alone does NOT make a plain canvas.draw() skip
        # an artist -- that skip only happens inside matplotlib's own blit
        # machinery, which this hand-rolled loop isn't using. Hide it
        # outright for the background capture instead (its pre-set state
        # would otherwise get baked into the cached background and show
        # through under every frame's actual content), then make it visible
        # again for the per-frame draws below.
        a.set_visible(False)

    canvas = fig.canvas
    canvas.draw()
    background = canvas.copy_from_bbox(fig.bbox)

    for a, _ in dynamic_artists:
        a.set_visible(True)

    # One frame written to disk and discarded at a time -- ffmpeg's
    # palettegen/paletteuse filters then read the PNG sequence directly
    # from disk (twice: once to build a shared palette, once to encode),
    # so peak memory is a single frame, not the whole animation.
    digits = len(str(n_frames - 1))
    with tempfile.TemporaryDirectory(prefix="femib_plot_frames_") as tmp_dir:
        pattern = os.path.join(tmp_dir, f"frame_%0{digits}d.png")
        for i in range(n_frames):
            update_fn(i)
            canvas.restore_region(background)
            for a, ax in dynamic_artists:
                ax.draw_artist(a)
            canvas.blit(fig.bbox)
            buf = numpy.asarray(canvas.buffer_rgba())
            # buf[:,:,:3], not .convert('RGB'): both drop the alpha channel,
            # but the numpy slice is a cheap view while .convert('RGB') is a
            # full pixel-by-pixel PIL color-mode conversion.
            Image.fromarray(buf[:, :, :3]).save(pattern % i)

        palette_path = os.path.join(tmp_dir, "palette.png")
        _run_ffmpeg(
            ["-framerate", str(fps), "-i", pattern, "-vf", "palettegen", palette_path]
        )
        _run_ffmpeg(
            [
                "-framerate",
                str(fps),
                "-i",
                pattern,
                "-i",
                palette_path,
                "-lavfi",
                "paletteuse",
                "-loop",
                "0",
                str(path),
            ]
        )


def _run_ffmpeg(args):
    result = subprocess.run(
        ["ffmpeg", "-y", "-loglevel", "error", *args],
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        raise RuntimeError(f"ffmpeg failed: {result.stderr}")
