import random


bases = ["A", "C", "G", "T"]
n = [10, 30]
m = [1000, 5000]
ins = 10


for ni in n:
	for mi in m:
		for i in range(ins):
			f = open("./instancias/" + str(ni) + "-" + str(mi) + "-" + str(i+1) + ".txt", "w")

			for i in range(ni):
				f.write( ''.join(random.choice(bases) for _ in range(mi)))
				f.write("\n")

			f.close()