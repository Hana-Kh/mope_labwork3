import random
import numpy
from scipy.stats import t,f
import tkinter

import tkinter.messagebox

root = tkinter.Tk()

x1_min = -30
x1_max = 0

x2_min = -15
x2_max = 35

x3_min = -30
x3_max = -25

xm_min = (x1_min + x2_min + x3_min) / 3
xm_max = (x1_max + x2_max + x3_max) / 3
y_min = 200 + xm_min
y_max = 200 + xm_max

xn = [[-1, -1, -1],
      [-1, 1, 1],
      [1, -1, 1],
      [1, 1, -1]]

x = [[-30, -15, -30],
     [-30, 35, -25],
     [0, -15, -25],
     [0, 35, -30]]

m = 2
y = [[random.randint(int(y_min), int(y_max)) for i in range(m)] for j in range(4)]


def kohren(dispersion, m):
    fisher = fisher_t(0.95, 1, (m - 1) * 4)
    gt = fisher/(fisher+(m-1)-2)
    gp = max(dispersion) / sum(dispersion)
    return gp < gt

def student(dispersion_reproduction, m, y_mean, xn):
    tt = 0
    f3 = (m - 1) * 4
    prob = 0.95
    
    x_vec = [i*0.0001 for i in range(int(5/0.0001))]
    par = 0.5 + prob/0.1*0.05
    for i in x_vec:
        if abs(t.cdf(i, f3) - par) < 0.000005:
            tt = i
            break
        
    dispersion_statistic_mark = (dispersion_reproduction / (4 * m)) ** 0.5

    beta = [1 / 4 * sum(y_mean[j] for j in range(4))]
    for i in range(3):
        b = 0
        for j in range(4):
            b += y_mean[j] * xn[j][i]
        beta.append(1 / 4 * b)

    te = []
    for i in beta:
        te.append(abs(i) / dispersion_statistic_mark)

    return te[0] > tt, te[1] > tt, te[2] > tt, te[3] > tt


def normalized_multiplier(x, y_mean):
    mx1 = (x[0][0] + x[1][0] + x[2][0] + x[3][0]) / 4
    mx2 = (x[0][1] + x[1][1] + x[2][1] + x[3][1]) / 4
    mx3 = (x[0][2] + x[1][2] + x[2][2] + x[3][2]) / 4
    my = sum(y_mean) / 4
    	
    a11	= (x[0][0] ** 2 + x[1][0] ** 2 + x[2][0] ** 2 + x[3][0] ** 2) / 4
    a22	= (x[0][1] ** 2 + x[1][1] ** 2 + x[2][1] ** 2 + x[3][1] ** 2) / 4
    a33	= (x[0][2] ** 2 + x[1][2] ** 2 + x[2][2] ** 2 + x[3][2] ** 2) / 4
    a12	= (x[0][0] * x[0][1] + x[1][0] * x[1][1] + x[2][0] * x[2][1] + x[3][0] * x[3][1]) / 4
    a13	= (x[0][0] * x[0][2] + x[1][0] * x[1][2] + x[2][0] * x[2][2] + x[3][0] * x[3][2]) / 4
    a23	= (x[0][1] * x[0][2] + x[1][1] * x[1][2] + x[2][1] * x[2][2] + x[3][1] * x[3][2]) / 4
    	
    a1 = (x[0][0] * y_mean[0] + x[1][0] * y_mean[1] + x[2][0] * y_mean[2] + x[3][0] * y_mean[3]) / 4	
    a2 = (x[0][1] * y_mean[0] + x[1][1] * y_mean[1] + x[2][1] * y_mean[2] + x[3][1] * y_mean[3]) / 4	
    a3 = (x[0][2] * y_mean[0] + x[1][2] * y_mean[1] + x[2][2] * y_mean[2] + x[3][2] * y_mean[3]) / 4	
    
    a = numpy.array([[1, mx1, mx2, mx3],
                     [mx1, a11, a12, a13],
                     [mx2, a12, a22, a23],
                     [mx3, a13, a23, a33]])
    c = numpy.array([my, a1, a2, a3])
    b = numpy.linalg.solve(a, c)
    return b


def fisher_t(prob, d, f3):
    x_vec = [i*0.001 for i in range(int(10/0.001))]
    for i in x_vec:
        if abs(f.cdf(i, 4-d, f3)-prob) < 0.0001:
            return i


def fisher(m, d, y_mean, yo, dispersion_reproduction):

    dispersion_ad = 0
    for i in range(4):
        dispersion_ad += (yo[i] - y_mean[i]) ** 2
        
    dispersion_ad = dispersion_ad * m / (4 - d)

    fp = dispersion_ad / dispersion_reproduction
    f3 = (m - 1) * 4
    
    return fp < fisher_t(0.95, d, f3)


while True:
    while True:            
        y_mean = []
        for i in range(4):
            y_mean.append(sum(y[i]) / m)
        
        dispersion = []
        for i in range(len(y)):
            dispersion.append(0)
            for j in range(m):
                dispersion[i] += (y_mean[i] - y[i][j]) ** 2
            dispersion[i] /= m

        dispersion_reproduction = sum(dispersion) / 4

        if kohren(dispersion, m):
            break
        else:
            m += 1
            for i in range(4):
                y[i].append(random.randint(int(y_min), int(y_max)))

    k = student(dispersion_reproduction, m, y_mean, xn)
    d = sum(k)
    
    b = normalized_multiplier(x, y_mean)
    b = [b[i] * k[i] for i in range(4)]

    yo = []
    for i in range(4):
        yo.append(b[0] + b[1] * x[i][0] + b[2] * x[i][1] + b[3] * x[i][2])
    
    if d == 4:
        m += 1
        for i in range(4):
            y[i].append(random.randint(int(y_min), int(y_max)))
        
    elif fisher(m, d, y_mean, yo, dispersion_reproduction):
        break
    else:
        m += 1
        for i in range(4):
            y[i].append(random.randint(int(y_min), int(y_max)))


tkinter.Label(text="x1").grid()

tkinter.Label(text="x2").grid(row=0, column=1)
tkinter.Label(text="x3").grid(row=0, column=2)
for i in range(m):
    tkinter.Label(text="yi" + str(i + 1)).grid(row=0, column=i + 3)
for i in range(len(x)):
    for j in range(len(x[i])):
        tkinter.Label(text=x[i][j]).grid(row=i + 1, column=j)
for i in range(len(y)):
    for j in range(len(y[i])):
        tkinter.Label(text=(y[i][j])).grid(row=i + 1, column=j + 3)
tkinter.Label(text="Рівняння регресії:").grid(columnspan=m + 3)
text = "y = " + "{0:.2f}".format(b[0])
for i in range(3):
    if b[i + 1] != 0:
        text = text + " + {0:.2f}".format(b[i + 1]) + " * x" + str(i + 1)

tkinter.Label(text=text).grid(columnspan=m + 3)
tkinter.Label(text="Перевірка:").grid(columnspan=m + 3)

for i in range(4):
    tkinter.Label(text="yc" + str(i + 1) + " =" + "{0:.2f}".format(y_mean[i])).grid(columnspan=m + 3)
    tkinter.Label(text="y" + str(i + 1) + " = " + "{0:.2f}".format(yo[i])).grid(columnspan=m + 3)

root.mainloop()
