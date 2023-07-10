import numpy as np

p1 = np.array([[0, 10], [10, 0], [0, -10], [-10, 0]])
p2 = np.array([[0, 5], [5, 0], [0, -5], [-5, 0]])
p = p2 - p1
l = np.array([3.64005, 4.272, 6.57647, 6.18466])
a = np.zeros(3)
m = np.ones(shape=(3, 3))

for i in range(3):
    a[i] = l[i]**2 - p[i, 0]**2 - p[i, 1]**2
    m[i, 0] = 2 * p[i, 0]
    m[i, 1] = 2 * p[i, 1]


np.linalg.inv(m) @ a

x = np.linalg.inv(m.T@m) @ m.T@a
x.T