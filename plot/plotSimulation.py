#!/bin/python

from plotUtils import parse_id, parse_input, calcPlotData, sortedTimesteps
import sys, os, glob, numpy, h5py, matplotlib, matplotlib.pyplot as pyplot
from matplotlib.animation import FuncAnimation

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

	def plot_snapshot(self,data,id_,save=False):
		fig, ax = pyplot.subplots(1,1)
		fig.set_tight_layout(True)

		x, y = data["X"][0], data["Y"][0]
		u, v = data["U"][0], data["V"][0]
		q = data["P"][0]

		if len(x) and len(u):
			ax.quiver(x,y,u,v,pivot='tail',units='xy')
			ax.set_title("velocity field: "+os.path.basename(id_))
		elif len(x) and len(q):
			try:
				cs = ax.tricontourf(x,y,q,levels=14,cmap=pyplot.cm.magma)
				fig.colorbar(cs,ax=ax)
			except (ValueError,RuntimeError):
				sc = ax.scatter(x,y,c=q,cmap=pyplot.cm.magma)
				fig.colorbar(sc,ax=ax)
			ax.set_title("scalar field: "+os.path.basename(id_))
		elif len(x):
			ax.scatter(x,y)
			ax.set_title("points: "+os.path.basename(id_))
		else:
			print("Nothing to plot for "+id_)
			pyplot.close(fig)
			return

		ax.axis('equal')

		if save:
			os.makedirs('gifs',exist_ok=True)
			fig.savefig('gifs/Simulation_'+os.path.basename(id_)+'.png',dpi=150)
		else:
			pyplot.show()

	def plot_animation(self,data,id_,save=False):
		fig, ax = pyplot.subplots(1,1)
		fig.set_tight_layout(True)

		def update(i):
			ax.cla()
			x, y = data["X"][i], data["Y"][i]
			u, v = data["U"][i], data["V"][i]
			if len(x):
				ax.set_xlim([x.min(),x.max()])
				ax.set_ylim([y.min(),y.max()])
			label = 'timestep {0}'.format(data["T"][i])
			ax.set_xlabel(label)
			if len(u):
				ax.quiver(x,y,u,v,pivot='tail',units='xy')
			ax.axis('equal')
			return ax

		anim = FuncAnimation(fig, update, frames=numpy.arange(0,len(data["T"])), interval=1)

		if save:
			os.makedirs('gifs',exist_ok=True)
			anim.save('gifs/Simulation_'+os.path.basename(id_)+'.gif', dpi=150, writer='imagemagick')
		else:
			pyplot.show()

if __name__ == '__main__':
	directory = '.'
	plotSimulation = PlotSimulation(directory)
	plotSimulation.do()
