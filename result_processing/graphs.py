import numpy as np
import matplotlib.pyplot as plt

while (True):
    res_num = input("Choose a test!\nPossible input options: 1, 2, 3\n")

    if res_num == "1" or res_num == "2" or res_num == "3":
        filename = "./input/input" + res_num + ".txt" # works when launched from the math_task_2 folder
        
        try:
            with open(filename, 'r') as file:
                lines = file.read().splitlines()

            if len(lines) == 2:
                block1 = lines[0].split()
                block2 = lines[1].split()

                if len(block1) == 3 and len(block2) == 3:
                    rho_L, v_L, p_L = map(float, block1)
                    rho_R, v_R, p_R = map(float, block2)
                else:
                    print("File format is incorrect in terms of the number of parameters in each block.")
            else:
                print("File format is incorrect. It should contain 2 lines.")

        except FileNotFoundError:
            print("File not found.")
        break
    else:
        print("Invalid test number.")

output_data = np.genfromtxt("./solution/output" + res_num + ".txt", delimiter=' ')
exact_solution = np.genfromtxt("./solution/riemann_solution/exact-sol" + res_num + ".txt", delimiter=' ')

# these are matplotlib.patch.Patch properties
props = dict(facecolor='white', alpha=0.3)

fig, (ax1, ax2, ax3) = plt.subplots(3, figsize=(8, 8))
fig.suptitle("Solution № " + res_num)

ax1.plot(output_data[:, 0], output_data[:, 2], color="#D11428")
ax1.plot(exact_solution[:, 0], exact_solution[:, 2], '--', color="#666666", label="Exact Solution")
ax1.set_ylabel("Density")

ax2.plot(output_data[:, 0], output_data[:, 3], color="#534491")
ax2.plot(exact_solution[:, 0], exact_solution[:, 3], '--', color="#666666", label="Exact Solution")
ax2.set_ylabel("Velocity")

ax3.plot(output_data[:, 0], output_data[:, 1], color="#50C878")
ax3.plot(exact_solution[:, 0], exact_solution[:, 1], '--', color="#666666", label="Exact Solution")
ax3.set_ylabel("Pressure")

ax3.set_xlabel("$x$")

textstr = '\n'.join((
    r'Initial conditions:',
    r'$\rho_L=%.1f$, $v_L=%.1f$, $p_L=%.1f,$' % (rho_L, v_L, p_L, ),
    r'$\rho_R=%.1f$, $v_R=%.1f$, $p_R=%.1f.$' % (rho_R, v_R, p_R, )))

fig.text(0, 1.4, textstr, transform=ax1.transAxes, fontsize=9.5, verticalalignment='top')
plt.show()
# plt.savefig("./graphs/solution" + res_num + ".png")
