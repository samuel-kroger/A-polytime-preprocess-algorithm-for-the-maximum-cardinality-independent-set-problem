import networkx as nx
import gurobipy as gp
import time
import csv
import os
import itertools

def float_to_str(float):
	return '{:.2f}'.format(float)

class problem_instance(object):
	
	def __init__(self, input_filename):

		self.filename = input_file[7:]

		self.graph = self.read_graph(input_filename)
		self.n = len(self.graph.nodes())
		self.m = len(self.graph.edges())
		self.graph_density = float_to_str(nx.density(self.graph))

		self.MIP_time_limit = 3600

	def read_graph(self, fname):

		if fname.endswith('mtx'):
			edges = []

			with open(fname, 'r') as f:
				reader = csv.reader(f, delimiter=' ')
				edges = [row for row in reader if len(row) == 2]
				f.close()

			graph = nx.Graph()
			graph.add_edges_from(edges, nodetype=int)

		else:
			#graph = nx.read_edgelist(fname, nodetype=int, data=(("Type", str),))
			#graph = nx.read_edgelist(fname, nodetype=int)
			graph = nx.read_adjlist(fname)

		graph.remove_edges_from(nx.selfloop_edges(graph))
		graph.remove_nodes_from(list(nx.isolates(graph)))

		#print(nx.info(graph))

		return graph

	def compuete_maximum_independent_set_of_simplicial_nodes(self):
		self.number_of_simplicial_iterations += 1
		simplicials = find_simplicials(self.graph)		
		subgraph_of_simplicials = self.graph.subgraph(simplicials)

		simplicial_fixings = []
		for component in nx.connected_components(subgraph_of_simplicials):
			component = list(component)
			simplicial_fixings.append(component[0])

		return simplicial_fixings

	def compuete_maximum_independent_set(self, method, relax):
		self.lower_bound = 'NA'
		self.upper_bound = 'NA'
		self.total_time = 'NA'
		self.method_time = 0
		self.number_of_simplicial_iterations = 0
		self.number_of_simplicial_fixings = 0
		self.lower_bound = 'NA'
		self.upper_bound = 'NA'


		start_time = time.time()
		m = gp.Model()
		m.Params.timeLimit = self.MIP_time_limit
		m.Params.MIPGap = 0


		if method == 'recursive_simplicial':
			method_time_start = time.time()
			simplicial_fixings = 'while loop start'
			number_fixed = 0
			while simplicial_fixings != []:
				simplicial_fixings = self.compuete_maximum_independent_set_of_simplicial_nodes()
				for node in simplicial_fixings:
					number_fixed += 1
					#for neighbor in self.graph.neighbors(node):
					neighbors = list(self.graph.neighbors(node))

					self.graph.remove_nodes_from(neighbors)
					self.graph.remove_node(node)
			method_time_end = time.time()
			
			
			self.method_time = float_to_str(method_time_end - method_time_start)
			
			self.number_of_simplicial_fixings = number_fixed


		#m.Params.Cutoff = 10
		if relax:
			m._X = m.addVars(self.graph.nodes(), vtype=gp.GRB.CONTINUOUS)
		else:
			m._X = m.addVars(self.graph.nodes(), vtype=gp.GRB.BINARY)

		
		
		
		if method == 'one_step_simplicial':
			method_time_start = time.time()
			simplicial_fixings = self.compuete_maximum_independent_set_of_simplicial_nodes()
			for fixing in simplicial_fixings:
				m._X[fixing].lb = 1
				for neighbor in self.graph.neighbors(fixing):
					m._X[neighbor].ub = 0
			method_time_end = time.time()
			
			self.method_time = float_to_str(method_time_end - method_time_start)
			
			self.number_of_simplicial_fixings = len(simplicial_fixings)
		
		m.setObjective(gp.quicksum(m._X), gp.GRB.MAXIMIZE)
		m.addConstrs(m._X[i] + m._X[j] <= 1 for (i,j) in self.graph.edges())

		#if relax:
		#	t = m.relax()
		#	t.optimize()
		#	m.update()
		m.optimize()

		max_independent_set = []

		end_time = time.time()

		self.lower_bound = m.getAttr(gp.GRB.Attr.ObjVal)
		self.upper_bound = m.getAttr(gp.GRB.Attr.ObjBound)

		if method == 'recursive_simplicial':
			self.lower_bound += number_fixed
			self.upper_bound += number_fixed
		self.total_time = float_to_str(end_time - start_time)
		#self.method_time = float_to_str(method_time_end - method_time_start)

		#return len(max_independent_set), number_of_simplicial_fixings, float_to_str(end_time-start_time)

	def get_root_relaxtaion(self, method):
		self.compuete_maximum_independent_set(method, True)
		self.root_relaxation = self.lower_bound

	def write_results(self, method):
		if not os.path.exists('./results.csv'):
			with open('./results.csv', 'w') as csvfile:
				writer = csv.writer(csvfile, delimiter=',')
				writer.writerow(['Instance', 'method', 'lower bound', 'upper bound', 'root relaxation', 'total_time', 'method_time', 'number of simplicial nodes fixed', 'number of iterations'])

		with open('./results.csv', 'a') as csvfile:
			writer = csv.writer(csvfile, delimiter=',')
			writer.writerow([self.filename, method, self.lower_bound, self.upper_bound, self.root_relaxation, self.total_time, self.method_time, self.number_of_simplicial_fixings, self.number_of_simplicial_iterations])

	def make_row(self, method):
		instance.get_root_relaxtaion(method)
		if method == 'recursive_simplicial':
			self.graph = nx.read_adjlist('./data/' + self.filename)
		instance.compuete_maximum_independent_set(method, False)
		
		instance.write_results(method)

	def calculate_results(self):
		method = 'none'
		self.make_row(method)
		method = 'one_step_simplicial'
		self.make_row(method)
		method = 'recursive_simplicial'
		self.make_row(method)

	def write_graph_info(self):

		if not os.path.exists('./graph_info.csv'):
			with open('./graph_info.csv', 'w') as csvfile:
				writer = csv.writer(csvfile, delimiter=',')
				writer.writerow(['Instance', 'number of nodes', 'number of edges', 'density'])

		with open('./graph_info.csv', 'a') as csvfile:
			writer = csv.writer(csvfile, delimiter=',')
			writer.writerow([self.filename, self.n, self.m, self.graph_density])

def float_to_str(float):
	return '{:.2f}'.format(float)



def find_simplicials(graph):
			
	simplicial_dict = {}
	simplicial_counter = 0
	for vertex in graph.nodes():
		simplicial_dict[vertex] = False
		neighbors = list(graph.neighbors(vertex))
		
		is_simplical = True
		for pair in itertools.combinations(neighbors,2):

			if not graph.has_edge(pair[0],pair[1]):
				is_simplical = False
				break

		if is_simplical:
			simplicial_dict[vertex] = True
			simplicial_counter += 1

	return [v for v in graph.nodes() if simplicial_dict[v] == True]


input_files = []
skip_files = []

for file in os.listdir('./data'):
	if file in skip_files:
		continue

	input_files.append('./data/' + file)



#input_files = ['./data/CA-CondMat.txt']




for input_file in input_files:
	print(input_file)
	instance = problem_instance(input_file)
	
	print('Instance name: ' + instance.filename)
	
	#instance.write_graph_info()
	instance.calculate_results()
	

'''
for input_file in input_files:
	print(input_file)
	
	graph = nx.read_adjlist(input_file)
	#print(len(graph.nodes()))
	
	
		writer.writerow([input_file[6:], max_independent_set, number_of_simplicial_fixings, run_time])
		max_independent_set, number_of_simplicial_fixings, run_time = maximum_independent_set(graph, True)
		writer.writerow([input_file[6:], max_independent_set, number_of_simplicial_fixings, run_time])
'''