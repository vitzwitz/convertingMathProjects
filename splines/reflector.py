"""
    Cubic Spline Interpolation - Natural Spline
    INPUT: X and Y are the vectors of given x-coordinates and y-coordinates respectively
    Author: Bri Miskovitz
    Converted from Matlab Project
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.axes as ax

def reflector(x,y,f,df):

        if type(x) != np.ndarray or type(y) != np.ndarray:
            x = np.asarray(x)
            y = np.asarray(y)
            f = np.asarray(f)
            df = np.asarray(df)



        # Stops if length(X) =/= length(Y)
        if len(x) != len(y):
            raise Warning("Vectors X and Y must be of the same length")

        # Number of points interpolating
        n = len(x)

        # Vector h with subintervals:
        # Step-size of x - COLUMN VECTOR
        h = []
        for j in range(n-2):
            h.append(x[j+1] - x[j])

        # Coefficient matrix A:
        # Create empty matrix
        A = np.zeros(n,n)

        # Natural Spline boundary conditions:
        A[0][0] = 1
        A[n-1][n-1] = 1

        for i in range(1,n-2):
            A[i][i-1] = h[i-1]
            # Diagonal elements
            A[i][i] = 2*(h[i-1] + h[i])
            A[i][i+1] = h[i]

        # Vector b - RHS vector - column vector
        b = []
        for i in range(1, n-2):
            b.append((3 / h[i]) * (y[i + 1] - y[i]) - (3 / h[i - 1]) * (y[i] - y[i - 1]))

        # Coefficient vector cj:  *** Needs : develop inverse matrix func ***
        cj = A/b

        # Coefficient vector bj:
        bj = []
        for i in range(n-2):
            bj.append((1 / h[i]) * (y[i + 1] - y[i]) - (1 / 3 * h[i]) * (2 * cj[i] + cj[i + 1]))


        # Coefficient vector dj:
        dj = []
        for i in range(n-2):
            dj.append((1/(3*h[i])) * (cj[i+1]-cj[i]))


        # Making a matrix P with all polynomials
        P = np.zeros(n-1,4)
        for i in range(n-2):
            P[i][1] = dj[i]
            P[i][2] = cj[i]
            P[i][3] = bj[i]
            P[i][4] = y[i]

        """
        ||**************************************************************||
        ||  Generating Shapes and Comparing them to Original Functions  ||
        ||**************************************************************||
        """
        # Data Points
        plt.plot(x,y, 'bo')
        resolution = 20

        # Issue : need to develop method to create functions using ndarrays and linspace
        s = np.zeros(n-2)
        for i in range(n-2):
            xs = np.linspace(x[i], x[i + 1], resolution)
            eq = y[i] + bj[i]*(xs-x[i]) + cj[i]*(xs-x[i])^2 + dj[i]*(x-x[i])^3
            np.append(s, eq)

            fig = plt.figure()
            plt.subplot(211)
            plt.plot(xs, s, 'k*', 'LineWidth', 2)
            plt.plot(xs, f, 'k*', 'LineWidth', 2)

            plt.legend('Data Points','Reflector','Function of Shape')
            plt.xlabel('x', fontsize=18)
            plt.ylabel('y', fontsize=16)
            plt.title('Comparing Reflector with Function of Shape using 20 equally spaced points', fontsize=20)
            jpg = "test" + str(i) + ".jpg"
            fig.savefig('jpg')