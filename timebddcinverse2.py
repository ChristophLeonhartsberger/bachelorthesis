import numpy as np
import matplotlib.pyplot as plt
import math

# Existing data
ndofs = [54189, 88101, 124986, 156141, 187104, 223692, 283458, 352743, 439650, 542166, 567690]
time1 = [5.696, 9.706, 13.947, 18.821, 23.371, 29.619, 39.164, 49.943, 65.468, 84.258, 89.947]
time2 = [10.449, 17.807, 29.035, 37.738, 49.090, 66.163, 96.115, 125.562, 170.525, 237.758, 252.718]

# Create the plot for Data Set 1
plt.figure()
plt.loglog(ndofs, time1, marker='o', linestyle='-', color='b', label='bddc')
plt.loglog(ndofs, time2, marker='o', linestyle='-.', color='r', label='inverse')
plt.loglog(ndofs, [1/5*time1[0] * (ndof / ndofs[0]) for ndof in ndofs], linestyle='--', color='g', label='O(n)')
#plt.loglog(ndofs, [1/3*time1[0] * ((ndof / ndofs[0])* math.log(ndof)) for ndof in ndofs], linestyle='--', color='c', label='O(n log n)')
plt.loglog(ndofs, [5*time1[0] * ((ndof / ndofs[0]) ** 2) for ndof in ndofs], linestyle='--', color='y', label='O(n^2)')
plt.xlabel('ndofs')
plt.ylabel('time (s)')
plt.grid(True)
plt.legend()

plt.show()