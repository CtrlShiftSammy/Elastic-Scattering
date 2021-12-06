from math import pi
import math


pi = 4 * math.atan(1.0)
data_file_c = open("Test/dat.txt", "w")
for i in range(11):
    data_file_c.write(str(math.sin(pi * i / 5)) + "\n")
data_file_c.close()