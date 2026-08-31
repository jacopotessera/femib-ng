#!/bin/python

import sys, numpy

def test():
	print("test")
	T = [
			None,
			"",
			"A",
			"1",
			"2,B",
			"1,2,1:2,3",
			"2:B",
			"5:3",
			"1,2,3:5:9",
			"5:5","1,2",
			"1:3",
			"1:4,6,8,10:15",
			"1,3:5,7"
		]
	for t in T:
		try:
			print(t,parse_id(t))
		except Exception as e:
			print(t)
			pass
	sims = [{"id": "a"},{"id": "b"},{"id": "c"},{"id": "d"},{"id": "e"}]
	inputs = [
				['plotUtils.py'],
				['plotUtils.py','1'],
				['plotUtils.py','1:2'],
				['plotUtils.py','A'],
				['plotUtils.py','1','op'],
				['plotUtils.py','1:2','op'],
				['plotUtils.py','A','op'],
				['plotUtils.py','A','op','wut'],
		]
	for i in inputs:
		sys.argv = i
		print(sys.argv,parse_input(sims))

def parse_id(i):
	ret = []
	try:
		if i is None or i == "":
			return []
		J = i.split(",")
		for j in J:
			K = j.split(":")
			if len(K)==1:
				ret.append(int(j))
			elif len(K)==2:
				k0 = int(K[0])
				k1 = int(K[1])+1
				if k0 < k1:
					for j in range(k0,k1):
						ret.append(j)
				else:
					raise Exception("invalid interval: " + j)
			else:
				raise Exception("invalid interval: " + j)
	except Exception as e:
		print(e)
		raise
	return list(set(ret))

def parse_input(sims):
	try:
		if len(sys.argv)==1:
			return {"op": "list"}
		elif len(sys.argv)==2:
			return {"op": "plot", "id_" : sims[int(sys.argv[1])-1]["id"]}
		elif len(sys.argv)==3:
			ids_ = parse_id(sys.argv[1])
			ids = []
			for id_ in ids_:
				ids.append(sims[id_-1]["id"])
			op = sys.argv[2]
			return {"op": op, "id_" : ids}
		else:
			raise Exception("cant parse input :|")
	except Exception as e:
		print(e)
		return {"op": "list"}

def _timestep_number(group_name):
	# Groups are named "timestep_<time>"
	return int(group_name.rsplit("_", 1)[-1])

def sortedTimesteps(source):
	if hasattr(source, "keys"):
		names = [k for k in source.keys() if k.startswith("timestep_")]
		names.sort(key=_timestep_number)
		return [source[n] for n in names]
	return list(source)

def calcPlotData(timesteps):
	groups = sortedTimesteps(timesteps)

	T, X, Y, U, V, P = [], [], [], [], [], []
	for grp in groups:
		name = grp.name.rsplit("/", 1)[-1]
		T.append(_timestep_number(name))

		if "x" in grp:
			x = grp["x"][:]
			X.append(x[:, 0])
			Y.append(x[:, 1])
		else:
			X.append(numpy.array([]))
			Y.append(numpy.array([]))

		if "u" in grp:
			u = grp["u"][:]
			U.append(u[:, 0])
			V.append(u[:, 1])
		else:
			U.append(numpy.array([]))
			V.append(numpy.array([]))

		P.append(grp["q"][:, 0] if "q" in grp else numpy.array([]))

	return {"T": T, "X": X, "Y": Y, "U": U, "V": V, "P": P}

if __name__ == '__main__':
	test()
