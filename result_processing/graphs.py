import numpy as np
import matplotlib.pyplot as plt

ind = input("Enter test number from 1, 2, 3, 4: \n")
test_data = np.genfromtxt("./results/res_test_" + ind +".txt", delimiter=' ')
if ind == "1":
    ro_l = 1.0; v_l = 0.0; p_l = 3.0
    ro_r = 1.0; v_r = 0.0; p_r = 1.0
    t = 0.18
elif ind == "2":
    # for test 2
    ro_l = 1.0; v_l = 1.0; p_l = 3.0
    ro_r = 1.0; v_r = -1.0; p_r = 1.0
    t = 0.1
elif ind == "3":
    # for test 3
    ro_l = 1.0; v_l = -0.1; p_l = 1.0
    ro_r = 1.0; v_r = 0.2; p_r = 1.0
    t = 0.1
elif ind == "4":
    ro_l = 1.0; v_l = 0.0; p_l = 1.0
    ro_r = 1.0; v_r = 0.0; p_r = 3.0
    t = 0.18


textstr = '\n'.join((
    r'Initial conditions:',
    r'$p_l=%.2f$, $p_r=%.2f$' % (p_l, p_r, ),
    r'$\rho_l=%.2f$, $\rho_r=%.2f$' % (ro_l, ro_r, ),
    r'$v_l=%.2f$, $v_r=%.2f$' % (v_l, v_r, ),
    r'$t=%.2f$' % (t, )))

# these are matplotlib.patch.Patch properties
props = dict(facecolor='white', alpha=0.3)

fig, (ax1, ax2, ax3) = plt.subplots(3, figsize=(8, 8))
fig.suptitle("Solution of test " + ind + ".")
ax1.plot(test_data[:, 0], test_data[:, 1], color="#FFC900")
# ax1.axvline(0, color='gray', linewidth=0.7, linestyle='dashed')
ax1.set_ylabel("Pressure")
ax2.plot(test_data[:, 0], test_data[:, 2], color="#B6FF00")
# ax2.axvline(0, color='gray', linewidth=0.7, linestyle='dashed')
ax2.set_ylabel("Density")
ax3.plot(test_data[:, 0], test_data[:, 3], color="#FF4A00")
# ax3.axvline(0, color='gray', linewidth=0.7, linestyle='dashed')
ax3.set_ylabel("Velocity")
ax3.set_xlabel("$x$")
fig.text(0.8, 1.45, textstr, transform=ax1.transAxes, fontsize=14, verticalalignment='top', bbox=props)
# plt.show()
plt.savefig("./results/res_test_" + ind + ".png")
