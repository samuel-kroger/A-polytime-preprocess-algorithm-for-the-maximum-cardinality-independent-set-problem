import networkx as nx
import gurobipy as gp
import matplotlib.pyplot as plt
import time
import csv
import os
import itertools
from pylab import Rectangle
import matplotlib.patches as mpatches
import random

def float_to_str(float):
	return '{:.2f}'.format(float)

class problem_instance(object):
	
	def __init__(self, input_filename):

		self.filename = input_filename

		self.graph = self.read_graph(self.filename)
		#self.graph = nx.path_graph(4)
		
		self.n = len(self.graph.nodes())
		self.m = len(self.graph.edges())
		self.graph_density = float_to_str(nx.density(self.graph))

		self.MIP_time_limit = 3* 3600
		self.MIP_print = 1

	def read_graph(self, fname):
		fname = './data/' + fname
		if fname.endswith('mtx'):
			edges = []

			with open(fname, 'r') as f:
				reader = csv.reader(f, delimiter=' ')
				edges = [row for row in reader if len(row) == 2]
				f.close()

			graph = nx.Graph()
			graph.add_edges_from(edges, nodetype=int)

		else:
			
			graph = nx.read_adjlist(fname)
			
			if fname.endswith('col'):
				
				for i in range(1, max([int(i) for i in graph.nodes()]) + 1):
					if str(i) not in list(graph.nodes()):
						graph.add_node(i)

			if fname.endswith('sam'):
				
				for i in range(1, max([int(i) for i in graph.nodes()]) + 1):
					if str(i) not in list(graph.nodes()):
						graph.add_node(i)
						print('happen')
						print(i)
				if fname.endswith('Erdos971.sam'):
					graph.add_node(491)
					graph.add_node(492)
				if fname.endswith('Erdos981.sam'):
					graph.add_node(484)
					graph.add_node(485)

		graph.remove_edges_from(nx.selfloop_edges(graph))
		
		graph = nx.convert_node_labels_to_integers(graph)
	


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

	def compuete_maximum_independent_set(self, method, relax, lp_fixing=False):
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
		m.Params.Output_Flag = self.MIP_print
		m.Params.timeLimit = self.MIP_time_limit
		m.Params.MIPGap = 0

		obj_value_tracker = 0
		
		
		if lp_fixing:
			method_time_start = time.time()
			self.number_of_simplicial_iterations += 1
			lp_fixings = []
			number_fixed = 0
			index = 0
			
			
			for vertex in self.relax_values:
				if vertex == 1:
					lp_fixings.append(index)
					obj_value_tracker += 1
				if vertex == 0:
					lp_fixings.append(index)
				index += 1

			self.graph.remove_nodes_from(lp_fixings)
			self.number_of_simplicial_fixings += len(lp_fixings)
			method_time_end = time.time()
			self.total_time += float_to_str(method_time_end - method_time_start)
		


			#return lp_fixings, number_fixed


		if method == 'recursive_simplicial':
			method_time_start = time.time()
			simplicial_fixings, number_fixed = self.recursive_fixing()
			simplicial_fixings = simplicial_fixings[1:]
			simplicial_fixings = [element for sublist in simplicial_fixings for element in sublist]
			method_time_end = time.time()
			
			
			self.method_time += method_time_end - method_time_start
			
			self.number_of_simplicial_fixings += number_fixed


		#m.Params.Cutoff = 10
		if relax:
			m._X = m.addVars(self.graph.nodes(), vtype=gp.GRB.CONTINUOUS, ub=1)
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
			self.method_time = method_time_end - method_time_start
			
			self.number_of_simplicial_fixings += len(simplicial_fixings)
		
		m.setObjective(gp.quicksum(m._X), gp.GRB.MAXIMIZE)
		m.addConstrs(m._X[i] + m._X[j] <= 1 for (i,j) in self.graph.edges())
		#m.addConstrs(m._X[i] <= 1 for i in self.graph.nodes()) 


		m.optimize()

		max_independent_set = []

		end_time = time.time()

		self.lower_bound = m.getAttr(gp.GRB.Attr.ObjVal) + obj_value_tracker
		self.upper_bound = m.getAttr(gp.GRB.Attr.ObjBound) + obj_value_tracker

		if method == 'recursive_simplicial':
			if not relax:
				self.lower_bound += len(simplicial_fixings)  
				self.upper_bound += len(simplicial_fixings) 
			

		if relax:
			#for v in m.getVars():
			#	print(f"{v.VarName} = {v.X}")
			self.relax_values = [v.x for v in m.getVars()]
		self.total_time = float_to_str(end_time - start_time)
		self.method_time = float_to_str(self.method_time)
	
	def recursive_fixing(self):
		simplicial_fixings = ['while loop start']
		number_fixed = 0
		while simplicial_fixings[-1] != []:
			simplicial_fixings.append(self.compuete_maximum_independent_set_of_simplicial_nodes())

			for node in simplicial_fixings[-1]:
				number_fixed += 1 + len(list(self.graph.neighbors(node)))
				neighbors = list(self.graph.neighbors(node))

				self.graph.remove_nodes_from(neighbors)
				self.graph.remove_node(node)
		return simplicial_fixings, number_fixed

	def get_root_relaxtaion(self, method):
		self.compuete_maximum_independent_set(method, True)
		self.root_relaxation = self.lower_bound

	def write_results(self, method):
		if not os.path.exists('./results.csv'):
			with open('./results.csv', 'w') as csvfile:
				writer = csv.writer(csvfile, delimiter=',')
				writer.writerow(['Instance', 'method', 'lower bound', 'upper bound', 'root relaxation', 'total_time', 'method_time', 'number nodes fixed', 'number of iterations'])

		with open('./results.csv', 'a') as csvfile:
			writer = csv.writer(csvfile, delimiter=',')
			writer.writerow([self.filename, method, self.lower_bound, self.upper_bound, self.root_relaxation, self.total_time, self.method_time, self.number_of_simplicial_fixings, self.number_of_simplicial_iterations])
			
	def make_row(self, method, lp_fixing=False):
		
		instance.get_root_relaxtaion('none')
		self.graph = self.read_graph(self.filename)
		instance.compuete_maximum_independent_set(method, False, lp_fixing)
		instance.write_results(method)
		self.graph = self.read_graph(self.filename)

	def calculate_results(self):
		#method = 'none'
		#self.make_row(method)
		#method = 'one_step_simplicial'
		#self.make_row(method)
		#method = 'recursive_simplicial'
		#self.make_row(method)
		method = 'LP fixing only'
		self.make_row(method, lp_fixing=True)
		method = 'recursive_simplicial'
		self.make_row(method, lp_fixing=True)

	def write_graph_info(self):

		if not os.path.exists('./graph_info.csv'):
			with open('./graph_info.csv', 'w') as csvfile:
				writer = csv.writer(csvfile, delimiter=',')
				writer.writerow(['Instance', 'number of nodes', 'number of edges', 'density'])

		with open('./graph_info.csv', 'a') as csvfile:
			writer = csv.writer(csvfile, delimiter=',')
			writer.writerow([self.filename, self.n, self.m, self.graph_density])

	def plot_fixings(self):
		self.number_of_simplicial_iterations = 0
		simplicial_fixings, number_of_simplicial_fixings = self.recursive_fixing()
		simplicial_fixings = simplicial_fixings[1:]
		self.graph = self.read_graph(self.filename)

		
		nodes = list(self.graph.nodes())
		removed_nodes = []
		pos = nx.spring_layout(self.graph, seed=5)  # positions for all nodes
		#plt.figure(figsize=(16, 6), dpi=80)
		fig, axes = plt.subplots(nrows=4, ncols=1, figsize=(10,15))
		ax = axes.flatten()
		#ax[0].title('super test')
		ax[0].set_axis_off()
		ax[1].set_axis_off()
		ax[2].set_axis_off()
		ax[3].set_axis_off()
		counter = 0

		#counter = 1
		for fixed_nodes in simplicial_fixings:
			neighbor_nodes = []
			for node in fixed_nodes:
				for potential_neighbor_node in self.graph.neighbors(node):
					if potential_neighbor_node in nodes:
						neighbor_nodes.append(potential_neighbor_node)
				
				nodes.remove(node)
				removed_nodes.append(node)
			#for node in neighbor_nodes:	
			#	if node in nodes:
			#		nodes.remove(node)

			options = {"edgecolors": "tab:gray", "node_size": 240, "alpha": 0.9}

			
			#ax[0].set_title('Computing independent set of simplicials')
			#ax[1].set_title('Removing the independent set of simplicials and their neighbors')

			#nx.draw_networkx_nodes(self.graph, pos, node_color="tab:blue", **options, ax=ax[0])
			if counter == 4:
				break
			nx.draw_networkx_nodes(nodes, pos, node_color="tab:blue", **options, ax=ax[counter])
			nx.draw_networkx_nodes(fixed_nodes, pos, node_color="tab:red", **options, ax=ax[counter])
			nx.draw_networkx_nodes(neighbor_nodes, pos, node_color="tab:olive", **options, ax=ax[counter])
			

			nx.draw_networkx_nodes(removed_nodes, pos, node_color="tab:olive", alpha=0.0, ax=ax[counter])
			subgraph = self.graph.subgraph(nodes + fixed_nodes)
			nx.draw_networkx_edges(subgraph, pos, width=2.0, alpha=.9, ax=ax[0])

			
			plt.tight_layout()
			#plt.axis("off",ax=ax[0])
			#plt.title(self.filename + ' ONE')
			
			#plt.savefig('./images/network_' + self.filename + '.png')
			#plt.show()
			
			for node in neighbor_nodes:
				if node in nodes:
					nodes.remove(node)
					removed_nodes.append(node)

			nx.draw_networkx_nodes(nodes, pos, node_color="tab:blue", **options, ax=ax[counter+1])
			#nx.draw_networkx_nodes(fixed_nodes, pos, node_color="tab:red", **options)
			#nx.draw_networkx_nodes(neighbor_nodes, pos, node_color="tab:olive", **options)
			nx.draw_networkx_nodes(removed_nodes, pos, node_color="tab:olive", alpha=0.0, ax=ax[counter+1])
			subgraph = self.graph.subgraph(nodes)
			nx.draw_networkx_edges(subgraph, pos, width=2.0, alpha=.9, ax=ax[2])
			nx.draw_networkx_edges(subgraph, pos, width=2.0, alpha=.9, ax=ax[counter+1])
			
			plt.tight_layout()
			#plt.axis("off",ax=ax[1])
			#plt.title(self.filename + ' TWO')

			#plt.show()
			counter += 2
		
		#plt.savefig('images/' + self.filename[:-4] + '_iteration_' + str(counter) + '.png')
		#plt.legend(['line1', 'line2', 'line3'], ['label1', 'label2', 'label3'])
		blue_patch = mpatches.Patch(color='tab:blue', label='Node')
		red_patch = mpatches.Patch(color='tab:red',  label='The red data')
		yellow_patch = mpatches.Patch(color='tab:olive', label='The red data')
		#plt.legend(handles=[blue_patch, red_patch, yellow_patch])
		plt.show()

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

def certain_graph():
	#graph = nx.cycle_graph(8)
	#graph.add_node(8)
	#graph.add_node(9)
	#edges= [(0,8), (7,8), (5,9), (4,9)]
	#graph.add_edges_from(edges)
	graph = nx.cycle_graph(5)
	graph.add_node(6)
	graph.add_edge(6,1)
	#graph.add_edge(6,3)
	
	#max(len(c) for c in nx.find_cliques(G))

	nx.draw(graph, with_labels=True)


	nx.write_adjlist(graph, './data/test.txt')
	instance = problem_instance('test.txt')
	instance.compuete_maximum_independent_set('recursive_simplicial', False)

	print(nx.is_chordal(graph))
	clique_number = max(len(c) for c in nx.find_cliques(graph))
	print(clique_number)
	print(nx.greedy_color(graph))
	plt.show()

input_files = []

for file in os.listdir('./data'):
	
	input_files.append(file)
	#instance = problem_instance(file)
	#instance.plot_fixings()
#input_files =  ['CA-CondMat.txt']
'''
#certain_graph()


i = 0
j = 0
while True:
	verticies = random.randint(200,300)
	probablity = random.randint(1,3)/300
	graph = nx.fast_gnp_random_graph(verticies, probablity)
	if len(graph.nodes()) == 0:
		continue
	nx.write_adjlist(graph, './data/test_graph.txt')
	instance = problem_instance('test_graph.txt')
	if nx.is_chordal(instance.graph):
		i += 1
		continue
	j += 1

	instance.compuete_maximum_independent_set('recursive_simplicial', False)
	if (i + j) % 100 == 0:
		print(i, ', ', j)
	if instance.lower_bound == 0 and instance.upper_bound == 0:
		nx.draw(instance.graph)
		plt.show()
		break
'''


for input_file in input_files:
	
	instance = problem_instance(input_file)
	print('Instance name: ' + instance.filename)
	
	instance.write_graph_info()
	instance.calculate_results()
	#instance.compuete_maximum_independent_set('none', True)
	#instance.compuete_maximum_independent_set('none', False, True)
	#instance.compuete_maximum_independent_set('none', True)
	#print(instance.relax_values)
	#print(len(instance.graph.nodes()))
	#instance.compuete_maximum_independent_set('recursive_simplicial', False)
	#print(len(instance.graph.nodes()))
	#print(instance.number_of_simplicial_fixings)
	#print(instance.number_of_simplicial_iterations)
