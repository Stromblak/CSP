import random


bases = ["A", "C", "G", "T"]
n = 10
m = 20


f = open("./instancias/" + str(n) + "-" + str(m) + ".txt", "w")

for i in range(n):
	f.write( ''.join(random.choice(bases) for _ in range(m)))
	f.write("\n")

f.close()