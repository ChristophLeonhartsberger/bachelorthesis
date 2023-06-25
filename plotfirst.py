import matplotlib.pyplot as plt

ndofs_set1 = [9228, 14778, 20697, 25812, 30588, 35421, 43134, 53658, 66132, 80079, 84081]
time_set1 = [0.969, 1.473, 2.232, 3.061, 3.805, 4.955, 6.011, 8.238, 10.503, 12.983, 13.401]

ndofs_set2 = [54189, 88101, 124986, 156141, 187104, 223692, 283458, 352743, 439650, 542166, 567690]
time_set2 = [5.696, 9.706, 13.947, 18.821, 23.371, 29.619, 39.164, 49.943, 65.468, 84.258, 89.947]

# Create the plot for Data Set 1
plt.figure()
plt.loglog(ndofs_set1, time_set1, marker='o', linestyle='-', color='b', label='first order')
plt.loglog(ndofs_set1, [1/5*time_set1[0] * (ndof / ndofs_set1[0]) for ndof in ndofs_set1], linestyle='--', color='g', label='O(n)')
plt.loglog(ndofs_set1, [5*time_set1[0] * ((ndof / ndofs_set1[0]) ** 2) for ndof in ndofs_set1], linestyle='--', color='y', label='O(n^2)')
plt.xlabel('ndofs')
plt.ylabel('time (s)')
plt.grid(True)
plt.legend()

# Create the plot for Data Set 2
plt.figure()
plt.loglog(ndofs_set2, time_set2, marker='o', linestyle='-', color='r', label='second order')
plt.loglog(ndofs_set2, [1/5*time_set2[0] * (ndof / ndofs_set2[0]) for ndof in ndofs_set2], linestyle='--', color='m', label='O(n)')
plt.loglog(ndofs_set2, [5*time_set2[0] * ((ndof / ndofs_set2[0]) ** 2) for ndof in ndofs_set2], linestyle='--', color='c', label='O(n^2)')
plt.xlabel('ndofs')
plt.ylabel('time (s)')
plt.grid(True)
plt.legend()

# Show the plots
plt.show()