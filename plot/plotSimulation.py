#!/bin/python

#
#	plotSimulation.py
#

from plotUtils import parse_id, parse_input, calcPlotData
import sys, pymongo, numpy, matplotlib, matplotlib.pyplot as pyplot
import itertools
from matplotlib.animation import FuncAnimation

class PlotSimulation():

	def __init__(self,db_name):
		self.db = pymongo.MongoClient('localhost', 27017)[db_name]
		self.simCollection = self.db['sim']
		#self.timestepCollection = self.db['timestep']
		self.plotDataCollection = self.db['plot_data']
		self.sims = self.simCollection.find({},{"id": 1}).sort([("id", pymongo.DESCENDING)])
		#self.plotData = {}
		#self.plotConfig = {"ffw": 40, "steps": 20, "area0": numpy.pi*0.6*0.6, "x_min": -1, "x_max": 1, "y_min": -1, "y_max": 1}
	
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
		n = 1
		for n,sim in enumerate(self.sims):
			id_ = str(sim["id"])
			try:
				TMAX = self.plotDataCollection.find({"id" : id_},{"time" : 1}).sort([("time", pymongo.DESCENDING)]).limit(1)[0]["time"]
			except:
				TMAX = 0
			print("\t* ["+str(n+1)+"] "+ id_ +" - "+str(TMAX+1)+" timesteps")
		print("")
		print("plot with:")
		print("\t$ make plot ID=n")
		sys.exit()

	def delete(self,ids):
		for id_ in ids:
			print("Deleting Simulation with id "+id_+" ...")
			self.simCollection.delete_many({"id" : id_})
			#self.timestepCollection.delete_many({"id" : id_})
			self.plotDataCollection.delete_many({"id" : id_})
		sys.exit()

	def save(self,ids):
		for id_ in ids:
			print("Saving Simulation with id "+id_+" ...")
			self.plot(id_,True)
		sys.exit()

	def plot(self,id_,save=False):

		#sim = self.simCollection.find({"id" : id_},{"S":1,"S.mesh" : 1,"parameters" : 1,"full": 1})
		#parameters = sim[0]["parameters"]
		#deltat = sim[0]["parameters"]["deltat"]
		TMAX = self.plotDataCollection.find({"id" : id_},{"time" : 1, "x" : 1}).sort([("time", pymongo.DESCENDING)]).limit(1)[0]["time"]
		timesteps = self.plotDataCollection.find({"id" : id_},{"time" : 1, "x" : 1,"u": 1,"q": 1})
		timesteps = list(timesteps)
		if not save:
			print("Plotting Simulation with id "+id_+" ...")
			#print("full:",sim[0]["full"])
			#print(parameters)
			print(TMAX)

			print(timesteps)
		#grid_x, grid_y = numpy.meshgrid(	numpy.arange(self.plotConfig["x_min"],self.plotConfig["x_max"]+0.1,0.1),
		#									numpy.arange(self.plotConfig["y_min"],self.plotConfig["y_max"]+0.1,0.1))

		#self.plotData = calcPlotData(timesteps,parameters,TMAX,self.plotConfig)

		fig, ax = pyplot.subplots(1,1)
		fig.set_tight_layout(True)
		#print('fig size: {0} DPI, size in inches {1}'.format(fig.get_dpi(), fig.get_size_inches()))


		#grid_x = self.plotData["W"][0]["x"]
		#grid_y = self.plotData["W"][0]["y"]
		grid = list(map(lambda x: x[0], timesteps[0]["u"]))
		grid_x = list(set(map( lambda x: x[0], grid)))
		grid_y = list(set(map( lambda x: x[1], grid)))
		print(grid_x)
		print(grid_y)
		grid_x, grid_y = numpy.meshgrid(grid_x, grid_y)
		print(grid_x)
		print(grid_y)
		#grid_x, grid_y = numpy.meshgrid(	numpy.arange(0,1+0.1,0.5),
		#	numpy.arange(0,1+0.1,0.5))
		#print(grid_x)
		#print(grid_y)		
		
		
		
		def update(i):
			uu = list(map(lambda x : x[1][0],timesteps[i]["u"] ))
			vv = list(map(lambda x : x[1][1],timesteps[i]["u"] ))
			ax.cla()
			#ax[0][0].plot(self.plotData["X"][0][i],self.plotData["X"][1][i],".-")
			ax.set_xlim([0,1])
			ax.set_ylim([0,1])
			#label = 'timestep {0},time {3:.2f}, area {1:.4f}, ratio {2:.4f}'.format(ffw*i,area/area0,_min/_max,ffw*i*0.01)
			label = 'time {0:.2f}'.format(i)	
			ax.set_xlabel(label)
			ax.quiver(grid_x,grid_y,uu,vv,pivot='tail',units='xy',scale=2)
			#circle = pyplot.Circle((0,0), 0.8, color='r',fill=False)
			#ax[0][0].add_artist(circle)
			#circle = pyplot.Circle((0,0), 0.6, color='r',fill=False)
			#ax[0][0].add_artist(circle)
			ax.axis('equal')

			#ax[0][1].cla()
			#ax[0][1].set_xlim([-1,1])
			#ax[0][1].set_ylim([-1,1])
			#label = 'time {0:.2f}'.format(self.plotData["T"][i])
			#ax[0][1].set_xlabel(label)
			#ax[0][1].quiver(grid_x,grid_y,self.plotData["U"][0][i],self.plotData["U"][1][i],pivot='mid',width=0.005)
			#ax[0][1].axis('equal')

			#ax[0][2].cla()
			#ax[0][2].set_xlim([-1,1])
			#ax[0][2].set_ylim([-1,1])
			#label = 'time {0:.2f}'.format(self.plotData["T"][i])	
			#ax[0][2].set_xlabel(label)
			#matplotlib.rcParams['xtick.direction'] = 'out'
			#matplotlib.rcParams['ytick.direction'] = 'out'
			#matplotlib.rcParams['contour.negative_linestyle'] = 'dashed'

			#pNorm = list(itertools.chain.from_iterable(self.plotData["P"][i]))
			#vmin = min(pNorm)
			#vmax = max(pNorm)

			#norm = matplotlib.colors.Normalize(vmax=vmax,vmin=vmin)
			#ax[0][2].contourf(grid_x,grid_y,self.plotData["P"][i], 8, alpha=.75, norm=norm, cmap=pyplot.cm.magma)
			#C = ax[0][2].contour(grid_x,grid_y,self.plotData["P"][i], 8, colors='black', linewidth=.5)
			#ax[0][2].clabel(C, inline=1, fontsize=10)
			#pyplot.xticks([]), pyplot.yticks([])

			#cs = ax[0][2].contour(grid_x,grid_y,self.plotData["P"][i])
			#ax[0][2].clabel(cs, inline=1, fontsize=10)
			#ax[0][2].axis('equal')

			#ax[1][0].cla()
			#ax[1][0].plot(self.plotData["T"][:i],self.plotData["A"][0][:i],"-")
			#ax[1][0].set_xlim([0,self.plotData["T"][-1]])
			#ax[1][0].set_ylim([0.0,1])

			#ax[1][1].cla()
			#ax[1][1].plot(self.plotData["T"][:i],list(map(lambda x: x["e"],self.plotData["A"][1]))[:i],"-")
			#ax[1][1].set_xlim([0,self.plotData["T"][-1]])
			#ax[1][1].set_ylim([0.78,1.22])

			#ax[1][2].cla()
			#ax[1][2].plot(self.plotData["T"][:i],list(map(lambda x: x["a"]/(x["a"]+x["A"]),self.plotData["A"][1]))[:i],"-")
			#ax[1][2].plot(self.plotData["T"][:i],list(map(lambda x: x["A"]/(x["a"]+x["A"]),self.plotData["A"][1]))[:i],"-")
			#ax[1][2].set_xlim([0,self.plotData["T"][-1]])
			#ax[1][2].set_ylim([0.4,0.6])

			return ax

		anim = FuncAnimation(fig, update, frames=numpy.arange(0, len(timesteps)), interval=1)
		
		if save:
			anim.save('gifs/Simulation_'+id_+'.gif', dpi=600, writer='imagemagick')
		else:
			pyplot.show()

if __name__ == '__main__':
	db_name = 'femib_test'
	plotSimulation = PlotSimulation(db_name)
	plotSimulation.do()	
	