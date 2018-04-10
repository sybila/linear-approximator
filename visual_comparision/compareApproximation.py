import os
import argparse
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate

parser = argparse.ArgumentParser(description='Call approximator and visualise results.')
parser.add_argument('function_file', type=str, default='input.txt', help='file with input functions')
parser.add_argument('point_count', type=str, default=1000, help='number of max samples')
parser.add_argument('segment_count', type=str, default=10, help='number of chosen samples')
parser.add_argument('fast', type=str, default=False, help='use fast approximation')
args = parser.parse_args()

def mod(a, b):
	return a % b

COLOURS = "bgrcmykw"
FUNDICT = {'mod': mod, 'tanh': np.tanh, 'pow': pow}

def assignColours(names):
	"""
	Function assignColours.

	Assignes a colour to a function name.
	"""
	colour_names = {}
	for i in range(len(names)):
		colour_names[names[i]] = COLOURS[i]
	return colour_names

def visualise(original_functions, approximated_functions, var_name, var_range):
	"""
	Function visualise.

	Shows original and approximated functions in one plot.
	"""
	names = [i for i,j in original_functions]
	colour_names = assignColours(names)
	var_values = np.linspace(*var_range, num=1000)
	plt.figure()
	legend_names = []

	for (name, expression) in original_functions:
		fun_values = list(map(lambda value: eval(expression, {var_name: value}, FUNDICT), var_values))
		plt.plot(var_values, fun_values, color=colour_names[name], linestyle='solid')
		legend_names.append(name + "_orig")

	for (name, values) in approximated_functions:
		plt.plot([i for i,j in values], [j for i,j in values], color=colour_names[name], linestyle=':', linewidth=5)
		legend_names.append(name + "_aprox")

	plt.legend(legend_names)
	plt.title('Approximation vs. original functions')
	plt.show()

def readFunctions(file):
	"""
	Function readFunctions.

	Inports functions from input file.
	"""
	var_name, var_range, functions = "", "", []
	with open(file) as lines:
		for line in lines:
			if "var" in line:
				parts = line.split(" ")
				var_name = parts[1]
				var_range = list(map(float, parts[-1][1:-2].split("..")))
			else:
				parts = line.split("=")
				name = parts[0].split(" ")[1]
				fun = parts[1].rstrip()
				functions.append((name, fun))
	return var_name, var_range, functions

def readApproximations(output, var_name):
	"""
	Function readApproximations.

	Reads approximated functions from output of approximator.
	"""
	print(output[0])
	samples = []
	approximations = []
	for line in output[1:]:
		print(line)
		parts = line.split(": ")
		if var_name == parts[0]:
			samples = parts[1]
		else:
			approximations.append((parts[0], eval(parts[1])))
	return samples, approximations

if __name__ == '__main__':

	var_name, var_range, functions = readFunctions(args.function_file)
	output = os.popen("bin/approximator {0} {1} {2} {3}".\
		format(args.function_file, args.point_count, args.segment_count, args.fast)).readlines()
	samples, approximations = readApproximations(output, var_name)
	visualise(functions, approximations, var_name, var_range)