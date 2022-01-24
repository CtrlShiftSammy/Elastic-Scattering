from math import log10
import numpy as np

def linear_interpolate(r):
    global data
    m1 = 0
    m2 =0
    for i in range(0, data_lines - 1):
        if (r >= data[i, 0] and r <= data[i+1, 0]):
            y = data[i-1, 1] + (r - data[i-1, 0]) * (data[i, 1] - data[i-1, 1]) / (data[i, 0] - data[i-1, 0])
            return y
def interpolate(r):
    global data
    m1 = 0
    m2 =0
    for i in range(0, data_lines - 1):
        if (r >= data[i, 0] and r <= data[i+1, 0]):
            if i == 0:
                m1 = (data[i+1, 1] - data[i, 1]) / (data[i+1, 0] - data[i, 0])
                m2 = (data[i+2, 1] - data[i, 1]) / (data[i+2, 0] - data[i, 0])
            elif i == data_lines - 2:
                m1 = (data[i+1, 1] - data[i-1, 1]) / (data[i+1, 0] - data[i-1, 0])
                m2 = (data[i+1, 1] - data[i, 1]) / (data[i+1, 0] - data[i, 0])
            else:
                m1 = (data[i+1, 1] - data[i-1, 1]) / (data[i+1, 0] - data[i-1, 0])
                m2 = (data[i+2, 1] - data[i, 1]) / (data[i+2, 0] - data[i, 0])
            a = (m2 - m1) / (2 * (data[i+1, 0] - data[i, 0]))
            b = m1 - 2 * a * data[i, 0]
            c = ((data[i, 1] - b * data[i, 0] - a * data[i, 0] * data[i, 0]) * (data[i + 1, 0] - r) + (data[i+1, 1] - b * data[i + 1, 0] - a * data[i + 1, 0] * data[i + 1, 0]) * (r - data[i, 0])) / ((data[i + 1, 0] - r) + (r - data[i, 0]))
            y = a * r * r + b * r + c
            return y
def interpolate_logrho(r):
    global data
    m1 = 0
    m2 =0
    for i in range(0, data_lines - 1):
        if (r >= data[i, 0] and r <= data[i+1, 0]):
            if i == 0:
                m1 = (data[i+1, 3] - data[i, 3]) / (data[i+1, 0] - data[i, 0])
                m2 = (data[i+2, 3] - data[i, 3]) / (data[i+2, 0] - data[i, 0])
            elif i == data_lines - 2:
                m1 = (data[i+1, 3] - data[i-1, 3]) / (data[i+1, 0] - data[i-1, 0])
                m2 = (data[i+1, 3] - data[i, 3]) / (data[i+1, 0] - data[i, 0])
            else:
                m1 = (data[i+1, 3] - data[i-1, 3]) / (data[i+1, 0] - data[i-1, 0])
                m2 = (data[i+2, 3] - data[i, 3]) / (data[i+2, 0] - data[i, 0])
            a = (m2 - m1) / (2 * (data[i+1, 0] - data[i, 0]))
            b = m1 - 2 * a * data[i, 0]
            c = ((data[i, 3] - b * data[i, 0] - a * data[i, 0] * data[i, 0]) * (data[i + 1, 0] - r) + (data[i+1, 3] - b * data[i + 1, 0] - a * data[i + 1, 0] * data[i + 1, 0]) * (r - data[i, 0])) / ((data[i + 1, 0] - r) + (r - data[i, 0]))
            y = a * r * r + b * r + c
            return 10**y
def interpolate_logr(r):
    r = log10(r)
    global data
    m1 = 0
    m2 =0
    for i in range(0, data_lines - 1):
        if (r >= data[i, 2] and r <= data[i+1, 2]):
            if i == 0:
                m1 = (data[i+1, 1] - data[i, 1]) / (data[i+1, 2] - data[i, 2])
                m2 = (data[i+2, 1] - data[i, 1]) / (data[i+2, 2] - data[i, 2])
            elif i == data_lines - 2:
                m1 = (data[i+1, 1] - data[i-1, 1]) / (data[i+1, 2] - data[i-1, 2])
                m2 = (data[i+1, 1] - data[i, 1]) / (data[i+1, 2] - data[i, 2])
            else:
                m1 = (data[i+1, 1] - data[i-1, 1]) / (data[i+1, 2] - data[i-1, 2])
                m2 = (data[i+2, 1] - data[i, 1]) / (data[i+2, 2] - data[i, 2])
            a = (m2 - m1) / (2 * (data[i+1, 2] - data[i, 2]))
            b = m1 - 2 * a * data[i, 2]
            c = ((data[i, 1] - b * data[i, 2] - a * data[i, 2] * data[i, 2]) * (data[i + 1, 2] - r) + (data[i+1, 1] - b * data[i + 1, 2] - a * data[i + 1, 2] * data[i + 1, 2]) * (r - data[i, 2])) / ((data[i + 1, 2] - r) + (r - data[i, 2]))
            y = a * r * r + b * r + c
            return y
def interpolate_logr_logrho(r):
    r = log10(r)
    global data
    m1 = 0
    m2 =0
    for i in range(0, data_lines - 1):
        if (r >= data[i, 2] and r <= data[i+1, 2]):
            if i == 0:
                m1 = (data[i+1, 3] - data[i, 3]) / (data[i+1, 2] - data[i, 2])
                m2 = (data[i+2, 3] - data[i, 3]) / (data[i+2, 2] - data[i, 2])
            elif i == data_lines - 2:
                m1 = (data[i+1, 3] - data[i-1, 3]) / (data[i+1, 2] - data[i-1, 2])
                m2 = (data[i+1, 3] - data[i, 3]) / (data[i+1, 2] - data[i, 2])
            else:
                m1 = (data[i+1, 3] - data[i-1, 3]) / (data[i+1, 2] - data[i-1, 2])
                m2 = (data[i+2, 3] - data[i, 3]) / (data[i+2, 2] - data[i, 2])
            a = (m2 - m1) / (2 * (data[i+1, 2] - data[i, 2]))
            b = m1 - 2 * a * data[i, 2]
            c = ((data[i, 3] - b * data[i, 2] - a * data[i, 2] * data[i, 2]) * (data[i + 1, 2] - r) + (data[i+1, 3] - b * data[i + 1, 2] - a * data[i + 1, 2] * data[i + 1, 2]) * (r - data[i, 2])) / ((data[i + 1, 2] - r) + (r - data[i, 2]))
            y = a * r * r + b * r + c
            return 10**y

data_lines = 182
data_file = open("Iteration 07/density_xe+2.dat", "r")
data = np.empty(shape = (data_lines, 4))
for i in range(0, data_lines):
    str1, str2 = (data_file.readline()).strip().split("  ")
    data[i,0] = float(str1.strip())
    data[i,1] = float(str2.strip())
    data[i,2] = log10(data[i,0])
    data[i,3] = log10(data[i,1])
data_file.close()
smooth_data_file = open("Iteration 07/density_xe+2_smooth.dat", "w")
r = 0.0
while (r < 27.0):
    r += 0.0001
    #smooth_data_file.write(str(r) + "  " + str(linear_interpolate(r)) + "  " + str(interpolate(r)) + "  " + str(interpolate_logr(r)) + "  " + str(interpolate_logrho(r)) + "  " + str(interpolate_logr_logrho(r)) + "\n")
    smooth_data_file.write(str(r) + "  " + str(interpolate_logr_logrho(r)) + "\n")
smooth_data_file.close()