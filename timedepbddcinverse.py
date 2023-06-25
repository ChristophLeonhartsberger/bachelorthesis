import numpy as np
import matplotlib.pyplot as plt

# Existing data
ndofs = [54189, 124986, 567690]
time1 = [5.696, 13.947, 89.947]
time2 = [10.449, 29.035, 252.718]

# Create the plot for Data Set 1
plt.figure()
plt.loglog(ndofs, time1, marker='o', linestyle='-', color='b', label='bddc')
plt.loglog(ndofs, [1/5*time1[0] * (ndof / ndofs[0]) for ndof in ndofs], linestyle='--', color='g', label='O(n)')
plt.loglog(ndofs, [5*time1[0] * ((ndof / ndofs[0]) ** 2) for ndof in ndofs], linestyle='--', color='y', label='O(n^2)')
plt.xlabel('ndofs')
plt.ylabel('time (s)')
plt.grid(True)
plt.legend()

# Create the plot for Data Set 2
plt.figure()
plt.loglog(ndofs, time2, marker='o', linestyle='-', color='r', label='inverse')
plt.loglog(ndofs, [1/5*time2[0] * (ndof / ndofs[0]) for ndof in ndofs], linestyle='--', color='m', label='O(n)')
plt.loglog(ndofs, [5*time2[0] * ((ndof / ndofs[0]) ** 2) for ndof in ndofs], linestyle='--', color='c', label='O(n^2)')
plt.xlabel('ndofs')
plt.ylabel('time (s)')
plt.grid(True)
plt.legend()

# Show the plots
plt.show()
