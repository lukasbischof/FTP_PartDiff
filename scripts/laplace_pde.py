# Solves the following PDE using the finite difference method:
# with u(x,y): -(u_xx + u_yy) = f(x,y)
# Domain: 0 < x < 1, 0 < y < 1
# Homogeneous Dirichlet BC: u(x,0) = 0, u(x,1) = 0, u(0,y) = 0, u(1,y) = 0
#
# discrete laplace operator:
# 1/h^2 * [u(i+1,j) + u(i-1,j) + u(i,j+1) + u(i,j-1) - 4u(i,j)]

import numpy as np
import matplotlib.pyplot as plt

# Nuber of knots in each direction
n = 70

# Step size
h = 1 / (n - 1)
X = np.linspace(0, 1, n)
Y = np.linspace(0, 1, n)


def f(x, y):
    return np.sinh(2 * np.pi * y) * np.cos(2 * np.pi * x)


def build_A():
    # Build the matrix A for the inner knots
    # The columns of A correspond to the approximated values of `u` at the inner knots
    # The rows of A correspond to the approximated values of the laplace operator at the inner knots:
    #     A     |  u_1,1, u_2,1, ..., u_n-2,1, u_1,2, ..., u_n-2,2, ..., u_1,n-2, ..., u_n-2,n-2
    # ----------+-------------------------------------------------------------------------------
    # u_1,1     |  -4     1      ...  0        1           0        ...     0      ...        0
    # u_2,1     |  1     -4      ...  ?        0           ?        ...     ?      ...        ?
    # ...       |  ...   ...     ...  ...      ...         ...      ...     ...    ...      ...
    # u_n-2,1   |
    # ...       |
    # u_n-2,n-2 |
    #
    # The whole matrix has to be multiplied by 1/h^2 additionally, to account for the discretization of the
    # laplace operator

    A = np.zeros(((n - 2) ** 2, (n - 2) ** 2))
    for yi in range(0, (n - 2)):
        for xi in range(0, (n - 2)):
            current_row_index = xi + yi * (n - 2)
            north = (xi, yi - 1)
            east = (xi + 1, yi)
            south = (xi, yi + 1)
            west = (xi - 1, yi)

            A[current_row_index, current_row_index] = -4
            if north[1] >= 0:
                A[current_row_index, north[0] + north[1] * (n - 2)] = 1
            if east[0] < n - 2:
                A[current_row_index, east[0] + east[1] * (n - 2)] = 1
            if south[1] < n - 2:
                A[current_row_index, south[0] + south[1] * (n - 2)] = 1
            if west[0] >= 0:
                A[current_row_index, west[0] + west[1] * (n - 2)] = 1
    return A * -(1 / h ** 2)


if __name__ == '__main__':
    A = build_A()
    b = np.zeros((n - 2) ** 2)
    for yi in range(0, (n - 2)):
        for xi in range(0, (n - 2)):
            current_row_index = xi + yi * (n - 2)
            b[current_row_index] = f(X[xi + 1], Y[yi + 1])
    u = np.linalg.solve(A, b)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    xv, yv = np.meshgrid(X, Y)
    Z = np.zeros((n, n))
    Z[1:n - 1, 1:n - 1] = u.reshape(n - 2, n - 2)
    ax.plot_surface(xv, yv, Z)
    plt.show()
