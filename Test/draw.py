import numpy as np

def data_smooth(r):
    global data
    m1 = 0
    m2 =0
    c = 0
    for i in range(0, 10):
        if (r >= i and r <= i+1):
            if i == 0:
                m1 = (data[i+1] - data[i]) / (i+1 - i)
                m2 = (data[i+2] - data[i]) / (i+2 - i)
            elif i == 9:
                m1 = (data[i+1] - data[i-1]) / (i+1 - (i-1))
                m2 = (data[i+1] - data[i]) / (i+1 - i)
            else:
                m1 = (data[i+1] - data[i-1]) / (i+1 - (i-1))
                m2 = (data[i+2] - data[i]) / (i+2 - i)
            #y = (m2 - m1) * r * r / (2 * (i+1 - i)) + (m1 * (i+1) - m2 * i) * r / (i+1 - i) + data[i] - ((m2 * i * i - m1 * i * i) / 2 + (m1 * i * (i+1) - m2 * i * i)) / (i+1 - i)
            a = (m2 - m1) / 2
            b = m1 - 2 * a * i
            c = ((data[i] - b * i - a * i * i) * (i + 1 - r) + (data[i+1] - b * (i+1) - a * (i+1) * (i+1)) * (r - i))
            y = a * r * r + b * r + c
            return y

data = np.empty(shape = (11))
data_file_c = open("Test/dat.txt", "r")
for i in range(11):
    data[i] = (data_file_c.readline()).strip()
data_file_c.close()

draw_file_c = open("Test/draw.txt", "w")
for i in range(100):
    draw_file_c.write(str(0.1 * i) + str(" ") + str(data_smooth(0.1 * i)) + "\n")
draw_file_c.close()