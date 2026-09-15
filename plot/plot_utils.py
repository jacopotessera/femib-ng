#!/bin/python

import numpy


def _timestep_number(group_name):
    # Groups are named "timestep_<time>"
    return int(group_name.rsplit("_", 1)[-1])


def sorted_timesteps(source):
    if hasattr(source, "keys"):
        names = [k for k in source.keys() if k.startswith("timestep_")]
        names.sort(key=_timestep_number)
        return [source[n] for n in names]
    return list(source)


def polygon_area(sx, sy):
    # Shoelace formula: signed area of the closed polygon (sx[i],sy[i]).
    # Structures are stored in a consistent orientation (see e.g.
    # femib::ib::build_ring's increasing-theta construction), so this is
    # positive in practice; abs() guards against a sign flip regardless.
    if len(sx) < 3:
        return numpy.nan
    return 0.5 * numpy.abs(numpy.sum(sx * numpy.roll(sy, -1) - numpy.roll(sx, -1) * sy))


def structure_aspect_ratio(sx, sy):
    # Bounding-box aspect ratio (width/height), matching ib_ns_demo.cpp's own
    # aspect_of() so the plotted time series matches the numbers logged to
    # stderr during the run.
    if len(sx) < 2:
        return numpy.nan
    width, height = sx.max() - sx.min(), sy.max() - sy.min()
    if height <= 0:
        return numpy.nan
    return width / height


def calc_plot_data(timesteps):
    groups = sorted_timesteps(timesteps)

    timestep_numbers, x_positions, y_positions = [], [], []
    u_velocities, v_velocities, pressures = [], [], []
    structure_x, structure_y, areas, aspect_ratios = [], [], [], []
    for grp in groups:
        name = grp.name.rsplit("/", 1)[-1]
        timestep_numbers.append(_timestep_number(name))

        if "x" in grp:
            x = grp["x"][:]
            x_positions.append(x[:, 0])
            y_positions.append(x[:, 1])
        else:
            x_positions.append(numpy.array([]))
            y_positions.append(numpy.array([]))

        if "u" in grp:
            u = grp["u"][:]
            u_velocities.append(u[:, 0])
            v_velocities.append(u[:, 1])
        else:
            u_velocities.append(numpy.array([]))
            v_velocities.append(numpy.array([]))

        pressures.append(grp["q"][:, 0] if "q" in grp else numpy.array([]))

        if "X" in grp:
            s = grp["X"][:]
            sx, sy = s[:, 0], s[:, 1]
            structure_x.append(sx)
            structure_y.append(sy)
            areas.append(polygon_area(sx, sy))
            aspect_ratios.append(structure_aspect_ratio(sx, sy))
        else:
            structure_x.append(numpy.array([]))
            structure_y.append(numpy.array([]))
            areas.append(numpy.nan)
            aspect_ratios.append(numpy.nan)

    return {
        "T": timestep_numbers,
        "X": x_positions,
        "Y": y_positions,
        "U": u_velocities,
        "V": v_velocities,
        "P": pressures,
        "SX": structure_x,
        "SY": structure_y,
        "AREA": areas,
        "ASPECT": aspect_ratios,
    }
