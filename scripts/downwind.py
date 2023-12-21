# -*- coding: utf-8 -*-
#
# Python script for solving the hyperbolic advection PDE with the (numerically unstable)
# downwind scheme. Since the method is unstable, it is not recommended to use it for
# any problem beyond these toy examples.
#

import numpy as np
import sympy as sp

sp.init_printing()

# Configuration
delta_x = 1/4
delta_t = 1/4
r = delta_t/delta_x
e = sp.Symbol('e')


# Boundary conditions
# u(x,0) = f(x)
def f(x):
    return e ** x


final_time = 1
required_steps = int(final_time/delta_t)
u = np.zeros((required_steps + 1, required_steps + 1)).astype(sp.Symbol)

# Apply initial conditions
u[0, :] = [f(j * delta_x) for j in range(0, required_steps + 1)]

if __name__ == '__main__':
    # Apply downwind scheme
    for k in range(1, required_steps + 1):
        for j in range(0, required_steps - k + 1):
            u[k, j] = (1 + r) * u[k - 1, j] - r * u[k - 1, j + 1]

    # Print results
    for k in range(0, required_steps + 1):
        print('k = ', k)
        for j in range(0, required_steps - k + 1):
            print('u[', k, ',', j, '] = ', sp.simplify(u[k, j]))
        print('\n')

    print(f"Which evaluates to {u[required_steps, 0].subs(e, np.e)}")
