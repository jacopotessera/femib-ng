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

    frames = []
    for i in range(n_frames):
        update_fn(i)
        canvas.restore_region(background)
        for a, ax in dynamic_artists:
            ax.draw_artist(a)
        canvas.blit(fig.bbox)
        buf = numpy.asarray(canvas.buffer_rgba())
        # buf[:,:,:3], not .convert('RGB'): both drop the alpha channel, but
        # the numpy slice is a cheap view while .convert('RGB') is a full
        # pixel-by-pixel PIL color-mode conversion.
        frames.append(Image.fromarray(buf[:, :, :3]))

    # Image.save(..., save_all=True) quantizes (24-bit RGB -> 8-bit palette)
    # each frame independently by default -- a full median-cut color search
    # per frame. Every frame here draws from the same handful of color
    # sources, so a palette built from one representative frame already
    # covers the rest well; reusing it turns every other frame's
    # quantization into a cheap nearest-color lookup instead of a fresh
    # search.
    palette_frame = frames[len(frames) // 2].quantize(colors=256)
    gif_frames = [f.quantize(palette=palette_frame) for f in frames]
    gif_frames[0].save(
        path, save_all=True, append_images=gif_frames[1:], duration=1000 // fps, loop=0
    )
