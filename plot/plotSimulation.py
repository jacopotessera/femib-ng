#!/bin/python

from plotUtils import parse_id, parse_input, calcPlotData, sortedTimesteps
import sys, os, glob, numpy, h5py, matplotlib, matplotlib.pyplot as pyplot
import matplotlib.colors
from matplotlib.animation import FuncAnimation
from PIL import Image

# TODO use typer
# TODO add instructions
class PlotSimulation():
	def __init__(self,directory="."):
		self.directory = directory
		self.sims = self._find_sims()

	def _find_sims(self):
		sims = []
		for path in sorted(glob.glob(os.path.join(self.directory,"*.h5"))):
			try:
				with h5py.File(path,"r") as f:
					sim_name = f.attrs.get("sim_name",os.path.basename(path))
					timesteps = sortedTimesteps(f)
					tmax = len(timesteps) - 1 if timesteps else -1
			except OSError as e:
				print("skipping "+path+": "+str(e))
				continue
			sims.append({"id": path, "sim_name": sim_name, "tmax": tmax})
		return sims

	def do(self):
		op = parse_input(self.sims)
		if op["op"] == "list":
			self.list()
		elif op["op"] == "plot":
			self.plot(op["id_"])
		elif op["op"] == "save":
			self.save(op["id_"])
		elif op["op"] == "delete":
			self.delete(op["id_"])
		else:
			print("invalid op :|")

	def list(self):
		print("Available Simulations:")
		for n,sim in enumerate(self.sims):
			steps = sim["tmax"]+1
			print("\t* ["+str(n+1)+"] "+str(sim["sim_name"])+" ("+sim["id"]+") - "+str(steps)+" timestep(s)")
		print("")
		print("plot with:")
		print("\t$ python3 plotSimulation.py n")
		sys.exit()

	def delete(self,ids):
		for id_ in ids:
			print("Deleting Simulation file "+id_+" ...")
			os.remove(id_)
		sys.exit()

	# TODO this is "make gif"
	def save(self,ids):
		# Force the non-interactive Agg backend for saving -- whatever GUI
		# backend matplotlib picked by default (Qt/Tk/...) routes every
		# frame through its toolkit's event loop and widget machinery for no
		# reason when nothing is ever shown on screen. Safe to switch here
		# since no figure has been created yet in this process.
		matplotlib.use('Agg',force=True)
		for id_ in ids:
			print("Saving Simulation "+id_+" ...")
			self.plot(id_,True)
		sys.exit()

	# TODO this is open plot window
	def plot(self,id_,save=False):
		print("Plotting Simulation "+id_+" ...")
		with h5py.File(id_,"r") as f:
			groups = sortedTimesteps(f)
			if not groups:
				print("No timestep data in "+id_)
				return
			data = calcPlotData(groups)

		if len(data["T"]) == 1:
			self.plot_snapshot(data,id_,save)
		else:
			self.plot_animation(data,id_,save)

	# Structure color is deliberately distinct from both the velocity
	# (viridis) and pressure (magma) colormaps so the material curve reads
	# clearly against either field.
	STRUCT_COLOR = '#d9750a'

	def _overlay_structure(self,ax,sx,sy):
		if len(sx) == 0:
			return
		xs = numpy.append(sx,sx[0])
		ys = numpy.append(sy,sy[0])
		ax.plot(xs,ys,color=self.STRUCT_COLOR,linewidth=2,zorder=5)

	def plot_snapshot(self,data,id_,save=False):
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
			print("Nothing to plot for "+id_)
			return
		if not panels:
			panels.append("structure only")

		ncols = len(panels)
		fig, axes = pyplot.subplots(1,ncols,figsize=(5*ncols,4.6),squeeze=False)
		axes = axes[0]
		fig.set_tight_layout(True)

		for ax,kind in zip(axes,panels):
			if kind=="velocity":
				mag = numpy.sqrt(u**2+v**2)
				quiv = ax.quiver(x,y,u,v,mag,pivot='tail',units='xy',cmap=pyplot.cm.viridis)
				fig.colorbar(quiv,ax=ax,label='|u|')
				ax.set_title("velocity")
			elif kind=="pressure":
				try:
					cs = ax.tricontourf(x,y,q,levels=14,cmap=pyplot.cm.magma)
					fig.colorbar(cs,ax=ax,label='p')
				except (ValueError,RuntimeError):
					sc = ax.scatter(x,y,c=q,cmap=pyplot.cm.magma)
					fig.colorbar(sc,ax=ax,label='p')
				ax.set_title("pressure")
			elif kind=="points":
				ax.scatter(x,y)
				ax.set_title("points")
			else:
				ax.set_title("structure")
			if has_struct:
				self._overlay_structure(ax,sx,sy)
			# adjustable='box' keeps equal aspect by resizing the drawn box,
			# not by rescaling the data limits -- axis('equal') does the
			# latter (adjustable='datalim'), which would silently zoom the
			# domain based on whatever this frame's data extent happens to be.
			ax.set_aspect('equal',adjustable='box')

		title = os.path.basename(id_)
		if has_struct and not numpy.isnan(data["AREA"][0]):
			title += "  (area={0:.4g}, aspect={1:.3f})".format(data["AREA"][0],data["ASPECT"][0])
		fig.suptitle(title)

		if save:
			os.makedirs('gifs',exist_ok=True)
			fig.savefig('gifs/Simulation_'+os.path.basename(id_)+'.png',dpi=150)
		else:
			pyplot.show()

	def plot_animation(self,data,id_,save=False):
		# Velocity/pressure are only saved every field_every steps (see
		# ib_ns_demo.cpp's dump_step); forward-fill each gap with the last
		# checkpoint's field so the quiver/pressure panels stay populated
		# every frame instead of flashing blank between checkpoints.
		last_X = last_Y = last_U = last_V = last_P = numpy.array([])
		for i in range(len(data["T"])):
			if len(data["U"][i]):
				last_X,last_Y,last_U,last_V,last_P = data["X"][i],data["Y"][i],data["U"][i],data["V"][i],data["P"][i]
			elif len(last_U):
				data["X"][i],data["Y"][i],data["U"][i],data["V"][i],data["P"][i] = last_X,last_Y,last_U,last_V,last_P

		has_u = any(len(u) for u in data["U"])
		has_q = any(len(q) for q in data["P"])
		has_struct = any(len(sx) for sx in data["SX"])

		top_panels = []
		if has_u: top_panels.append("velocity")
		if has_q: top_panels.append("pressure")
		if not top_panels and not has_struct:
			print("Nothing to plot for "+id_)
			return

		nrows = 2 if has_struct else 1
		ncols = max(len(top_panels),2 if has_struct else 1)
		# dpi=80 (matplotlib's default is ~100): this is read at a glance, not
		# zoomed into, and nearly every remaining per-frame cost (drawing the
		# quiver/scatter, GIF quantization, GIF frame-diffing/encoding) scales
		# with pixel count -- this used to be passed to anim.save() but was
		# lost when that call was replaced by the manual blit loop below.
		fig, axes = pyplot.subplots(nrows,ncols,figsize=(5*ncols,4.6*nrows),squeeze=False,dpi=80)
		# NOT set_tight_layout(True): that installs a layout engine that
		# recomputes the tight bounding box of EVERY text artist (every tick
		# label, every axis, on every subplot) on every single draw -- with
		# the per-frame 'timestep N' xlabel update below changing text each
		# frame, that turned into ~120s of a ~280s render doing nothing but
		# bbox math. A single one-shot tight_layout() call once, right
		# before the animation starts, gets the same spacing without paying
		# for it 361 more times.
		top_axes = axes[0]
		for ax in top_axes[len(top_panels):]:
			ax.axis('off')

		# Fixed color scales across the whole animation (not per-frame) so
		# color is comparable frame-to-frame -- otherwise a quiet frame and a
		# vigorous frame would both autoscale to "full brightness", hiding
		# exactly the amplitude change a non-stationary run is meant to show.
		mags = [numpy.sqrt(u**2+v**2) for u,v in zip(data["U"],data["V"]) if len(u)]
		vmin, vmax = min((m.min() for m in mags), default=0.0), max((m.max() for m in mags), default=1.0)
		if vmax <= vmin: vmax = vmin+1e-9
		vcmap, vnorm = pyplot.cm.viridis, matplotlib.colors.Normalize(vmin=vmin,vmax=vmax)

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
		grid_spacing = numpy.min(numpy.diff(numpy.sort(xs))) if len(xs)>1 else 0.04
		quiver_scale = vmax/(3*grid_spacing) if vmax>0 else 1.0
		quiver_width = 0.15*grid_spacing

		all_p = [p for p in data["P"] if len(p)]
		pmin, pmax = min((p.min() for p in all_p), default=0.0), max((p.max() for p in all_p), default=1.0)
		if pmax <= pmin: pmax = pmin+1e-9
		pcmap, pnorm = pyplot.cm.magma, matplotlib.colors.Normalize(vmin=pmin,vmax=pmax)

		vel_ax = pres_ax = None
		idx = 0
		if has_u:
			vel_ax = top_axes[idx]; idx += 1
			sm = pyplot.cm.ScalarMappable(cmap=vcmap,norm=vnorm); sm.set_array([])
			fig.colorbar(sm,ax=vel_ax,label='|u|')
		if has_q:
			pres_ax = top_axes[idx]; idx += 1
			sm2 = pyplot.cm.ScalarMappable(cmap=pcmap,norm=pnorm); sm2.set_array([])
			fig.colorbar(sm2,ax=pres_ax,label='p')

		area_ax = aspect_ax = area_marker = aspect_marker = None
		if has_struct:
			area_ax, aspect_ax = axes[1][0], axes[1][1]
			for ax in axes[1][2:]:
				ax.axis('off')
			T = numpy.array(data["T"])
			AREA = numpy.array(data["AREA"],dtype=float)
			ASPECT = numpy.array(data["ASPECT"],dtype=float)
			area_ax.plot(T,AREA,color=self.STRUCT_COLOR)
			area_ax.set_title("structure area")
			area_ax.set_xlabel("timestep")
			(area_marker,) = area_ax.plot([],[],'o',color='black')
			aspect_ax.plot(T,ASPECT,color=self.STRUCT_COLOR)
			aspect_ax.set_title("structure aspect ratio")
			aspect_ax.set_xlabel("timestep")
			(aspect_marker,) = aspect_ax.plot([],[],'o',color='black')

		# The old update() called ax.cla() then rebuilt the quiver/scatter/
		# structure-line from scratch every frame -- most of the per-frame
		# cost isn't the actual pixel rendering, it's Python-level artist
		# construction (quiver in particular recomputes arrow-head geometry
		# for every point) and cla() tearing down and rebuilding axis
		# ticks/spines. The sampling grid (x,y) is the same fixed regular
		# grid on every field-carrying frame (finite_element_space::plot
		# always resamples the same box at the same delta), so it only needs
		# to be read once; each frame then just pushes new U/V/color data
		# into the SAME artists via set_UVC/set_array/set_data, which is far
		# cheaper than recreating them.
		field_idx = next((i for i in range(len(data["T"])) if len(data["X"][i])), None)
		x0,y0 = (data["X"][field_idx],data["Y"][field_idx]) if field_idx is not None else (numpy.array([]),numpy.array([]))
		zeros0 = numpy.zeros_like(x0)

		vel_quiv = vel_struct_line = None
		if vel_ax is not None:
			vel_ax.set_title("velocity")
			vel_ax.set_xlim(0,1); vel_ax.set_ylim(0,1)
			# adjustable='box' keeps this exact xlim/ylim -- axis('equal')
			# (adjustable='datalim') would instead stretch the domain to fit
			# each frame's quiver-arrow extent, making it drift frame to
			# frame even though the fluid mesh itself never moves.
			vel_ax.set_aspect('equal',adjustable='box')
			if len(x0):
				vel_quiv = vel_ax.quiver(x0,y0,zeros0,zeros0,zeros0,pivot='tail',units='xy',
				                          cmap=vcmap,norm=vnorm,scale=quiver_scale,
				                          scale_units='xy',width=quiver_width)
			(vel_struct_line,) = vel_ax.plot([],[],color=self.STRUCT_COLOR,linewidth=2,zorder=5)

		pres_sc = pres_struct_line = None
		if pres_ax is not None:
			pres_ax.set_title("pressure")
			pres_ax.set_xlim(0,1); pres_ax.set_ylim(0,1)
			pres_ax.set_aspect('equal',adjustable='box')
			if len(x0):
				pres_sc = pres_ax.scatter(x0,y0,c=zeros0,cmap=pcmap,norm=pnorm,s=8)
			(pres_struct_line,) = pres_ax.plot([],[],color=self.STRUCT_COLOR,linewidth=2,zorder=5)

		def update(i):
			u, v = data["U"][i], data["V"][i]
			q = data["P"][i]
			sx, sy = data["SX"][i], data["SY"][i]
			t = data["T"][i]
			xs = numpy.append(sx,sx[0]) if len(sx) else sx
			ys = numpy.append(sy,sy[0]) if len(sy) else sy

			if vel_quiv is not None and len(u):
				vel_quiv.set_UVC(u,v,numpy.sqrt(u**2+v**2))
			if vel_struct_line is not None:
				vel_struct_line.set_data(xs,ys)
			if vel_ax is not None:
				# .xaxis.label.set_text(), not set_xlabel(): set_xlabel also
				# recomputes the label's on-axes position every call (to stay
				# clear of the tick labels below it) -- pure overhead here
				# since that position never actually needs to change frame
				# to frame, only the text content does.
				vel_ax.xaxis.label.set_text('timestep {0}'.format(t))

			if pres_sc is not None and len(q):
				pres_sc.set_array(q)
			if pres_struct_line is not None:
				pres_struct_line.set_data(xs,ys)
			if pres_ax is not None:
				pres_ax.xaxis.label.set_text('timestep {0}'.format(t))

			if area_ax is not None and not numpy.isnan(data["AREA"][i]):
				area_marker.set_data([t],[data["AREA"][i]])
				aspect_marker.set_data([t],[data["ASPECT"][i]])

			return []

		# Set the widest label text ("timestep <last>") before the one-shot
		# tight_layout() below, so spacing is computed for the actual worst
		# case up front rather than shifting slightly as the number of
		# digits in the timestep grows over the animation.
		last_t = data["T"][-1]
		if vel_ax is not None: vel_ax.set_xlabel('timestep {0}'.format(last_t))
		if pres_ax is not None: pres_ax.set_xlabel('timestep {0}'.format(last_t))
		fig.tight_layout()

		if not save:
			anim = FuncAnimation(fig, update, frames=numpy.arange(0,len(data["T"])), interval=1)
			pyplot.show()
			return

		# Manual blitting instead of anim.save(): matplotlib's Animation.save()
		# always does a FULL canvas redraw per output frame -- recomputing
		# every tick label's text layout, every axis, both colorbars, on
		# every single frame -- regardless of any blit setting on
		# FuncAnimation (blit only ever applies to the on-screen/interactive
		# path, never to file output). Measured via cProfile at ~120-150s of
		# a ~170s render for this exact animation, almost all of it in
		# matplotlib's text/bbox layout code, not actual pixel drawing. None
		# of that static decoration (axes, ticks, colorbars, the area/aspect
		# line plots) changes frame to frame, so it's rendered once here,
		# cached as a bitmap, and each frame only draws the handful of
		# artists that actually change (quiver, pressure scatter, structure
		# line, the two markers, the two xlabels) on top of that cached
		# background -- the standard matplotlib blitting pattern, just
		# driven by hand since Animation.save() doesn't use it.
		os.makedirs('gifs',exist_ok=True)
		# (artist, axes) pairs, not bare artists: a Text label's own .axes
		# attribute isn't reliably set the way a plotted Line2D/PathCollection's
		# is, so draw_artist needs to be called via the axes we already know
		# each one belongs to.
		dynamic_artists = [(a,ax) for a,ax in (
			(vel_quiv,vel_ax), (vel_struct_line,vel_ax),
			(pres_sc,pres_ax), (pres_struct_line,pres_ax),
			(area_marker,area_ax), (aspect_marker,aspect_ax),
			(vel_ax.xaxis.label if vel_ax is not None else None,vel_ax),
			(pres_ax.xaxis.label if pres_ax is not None else None,pres_ax),
		) if a is not None]
		# set_animated(True) alone does NOT make a plain canvas.draw() skip
		# these artists -- that skip only happens inside matplotlib's own
		# blit machinery, which this hand-rolled loop isn't using. Hide them
		# outright for the background capture instead (their pre-set state,
		# e.g. the "timestep <last>" label, would otherwise get baked into
		# the cached background and show through under every frame's actual
		# label), then make them visible again for the per-frame draws below.
		for a,_ in dynamic_artists:
			a.set_animated(True)
			a.set_visible(False)

		canvas = fig.canvas
		canvas.draw()
		background = canvas.copy_from_bbox(fig.bbox)

		for a,_ in dynamic_artists:
			a.set_visible(True)

		frames = []
		for i in range(len(data["T"])):
			update(i)
			canvas.restore_region(background)
			for a,ax in dynamic_artists:
				ax.draw_artist(a)
			canvas.blit(fig.bbox)
			buf = numpy.asarray(canvas.buffer_rgba())
			# buf[:,:,:3], not .convert('RGB'): both drop the alpha channel,
			# but the numpy slice is a cheap view while .convert('RGB') is a
			# full pixel-by-pixel PIL color-mode conversion -- measured via
			# cProfile as a real cost at 361 calls, one per frame.
			frames.append(Image.fromarray(buf[:,:,:3]))

		# Image.save(..., save_all=True) quantizes (24-bit RGB -> 8-bit
		# palette) each frame independently by default -- a full median-cut
		# color search per frame, measured via cProfile at ~15s of a ~28s
		# save, the single largest cost in the whole script by far. Every
		# frame here draws from the same handful of color sources (two fixed
		# colormaps, a white background, one overlay color), so a palette
		# built from one representative frame already covers the rest well;
		# reusing it turns every other frame's quantization into a cheap
		# nearest-color lookup instead of a fresh search.
		palette_frame = frames[len(frames)//2].quantize(colors=256)
		gif_frames = [f.quantize(palette=palette_frame) for f in frames]

		# 15fps: no explicit fps was ever set on the old anim.save() call
		# (PIL reports the resulting GIF's frame duration as unset, meaning
		# viewers fell back to their own default, commonly ~10fps) -- this
		# picks an explicit, reasonable pace instead of leaving it implicit.
		gif_frames[0].save('gifs/Simulation_'+os.path.basename(id_)+'.gif',
		                save_all=True, append_images=gif_frames[1:],
		                duration=1000//15, loop=0)

if __name__ == '__main__':
	directory = '.'
	plotSimulation = PlotSimulation(directory)
	plotSimulation.do()
