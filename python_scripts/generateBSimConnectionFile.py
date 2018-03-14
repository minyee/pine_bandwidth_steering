# This file generates the script for the connection files to be used in BookSim
import sys

def readSolutionMatrix(filename):
	f = open(filename, 'r')
	lines = f.readlines()
	cnt = 0
	size = int(lines[0])
	matrix = []
	row = 0
	for line in lines:
		if cnt == 0:
			cnt += 1
			continue
		else:
			words = line.split()
			arr = []
			for word in words:
				arr.append(int(word))
			matrix.append(arr)
	f.close()
	return matrix

def generateBSimScript(matrix, generalLatency, scriptFileName):
	f = open(scriptFileName, 'w'):
	f.writeline("lat %d" % generalLatency)
	size = len(matrix)
	for i in range(size):
		strBuilder = "%d " % i 
		for j in range(size):
			if matrix[i][j] > 0:
				strBuilder.append("%d (%d) " % (j, generalLatency/matrix[i][j]))
		f.writeline(strBuilder)
	f.close()
	return

solnFileName = sys.argv[1]
outputScriptFileName = sys.argv[2]
latency = 3 # default latency of 3 cycles
if len(sys.argv > 3):
	latency = int(sys.argv[3])
matrix = readSolutionMatrix(solnFileName)
generateBSimScript(matrx, latency, outputScriptFileName)
print "Successfully generated BookSim connection script"
