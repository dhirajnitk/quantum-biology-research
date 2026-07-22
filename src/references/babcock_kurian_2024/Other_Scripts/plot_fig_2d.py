from matplotlib import pyplot as plt
import numpy as np

# ASSUMING past spiral #20 (third 1)
#x_ax_real_num = np.array([1,2,3,4,5,6,7,8,9,10,11,12,13,20,30,40,50,60,70,80,90,100])

# ASSUMING past spiral #30 (val. 147)
x_ax = np.array([1,2,3,4,5,6,7,8,9,10,11,12,13,20,30,40,50,60,70,80])# ,90,100])

#y_ax_real_num = np.array([39,12,11,7,5,5,3,3,3,3,2,1,1,1,1,1,1,1,1,1,1,1])
#y_ax_z_proj_num = np.array([21,4,7,10,13,27,32,37,42,47,52,57,62,97,147,197,247,297,347,397,447,497])
#y_ax_real_approx = np.array([16,43,68,92,116,139,164,191,210,222,252,276,301,462,694,946,1170,1409,1656,1879,2135,2193])
#y_ax_z_proj_approx = np.array([19,204,290,14,9,10,667,35,60,100,271,283,304,27,29,62,1351,349,122,352,8290,1019])

y_ax_real_num = np.array([14.029, 21.273, 37.794, 51.511, 63.649, 83.938, 99.959, 116.53, 129.24, 155.31, 192.10, 217.507, 237.78, 345.10, 445.26, 500.55, 530.44, 547.633, 558.450, 565.72])
y_ax_z_proj_num = np.array([3.338, 5.519, 7.672, 10.332, 12.551, 15.184, 18.125, 20.967, 23.702, 26.325, 28.832, 31.223, 33.498, 45.922, 52.697, 53.275, 58.689, 58.512, 62.296, 62.087])
y_ax_real_approx = np.array([13.973, 27.902, 40.561, 52.582, 64.193, 75.505, 86.579, 97.665, 107.92, 118.61, 126.07, 134.18, 141.86, 190.46, 246.54, 282.54, 303.37, 315.64,323.41, 331.63])
y_ax_z_proj_approx = np.array([4.379, 4.718, 6.257, 7.393, 8.590, 9.257, 11.182, 12.047, 13.536, 31.046, 17.004, 15.425, 12.939, 24.211, 24.480, 22.828, 26.964, 28.595, 28.314, 27.474])

print(np.size(y_ax_real_num))
print(np.size(y_ax_z_proj_num))
print(np.size(y_ax_real_approx))
print(np.size(y_ax_z_proj_approx))
print(np.size(x_ax))

plt.xlim(0,104)
plt.xticks(ticks=np.linspace(0,104,9), labels=np.linspace(0,104,9))

plt.ylim(0.8,2200)
plt.yscale("log")

#plt.scatter(x_ax, y_ax_real_num, marker='x', color='b', linewidth=2, s=80, label="Exact diagonalization of realistic MT dipoles")
plt.scatter(x_ax, y_ax_real_num, marker='x', color='b', linewidth=2, s=80)
plt.plot(x_ax, y_ax_real_num, marker='x', color='b', label="Exact diagonalization of realistic MT dipoles")

#plt.scatter(x_ax, y_ax_real_approx, marker='x', color='c', linewidth=2, s=80, label="Approximate diagonalization of realistic MT dipoles")
plt.scatter(x_ax, y_ax_real_approx, marker='x', color='c', linewidth=2, s=80)
plt.plot(x_ax, y_ax_real_approx, marker='x', color='c', label="Approximate diagonalization of realistic MT dipoles")

#plt.scatter(x_ax, y_ax_z_proj_num, marker='x', color='r', linewidth=2, s=80, label="Exact diagonalization of z-projected MT dipoles")
plt.scatter(x_ax, y_ax_z_proj_num, marker='x', color='r', linewidth=2, s=80)
plt.plot(x_ax, y_ax_z_proj_num, marker='x', color='r', label="Exact diagonalization of z-projected MT dipoles")

#plt.scatter(x_ax, y_ax_z_proj_approx, marker='x', color='m', linewidth=2, s=80, label="Approximate diagonalization of z-projected MT dipoles")
plt.scatter(x_ax, y_ax_z_proj_approx, marker='x', color='m', linewidth=2, s=80)
plt.plot(x_ax, y_ax_z_proj_approx, marker='x', color='m', label="Approximate diagonalization of z-projected MT dipoles")

plt.legend()

plt.show()
