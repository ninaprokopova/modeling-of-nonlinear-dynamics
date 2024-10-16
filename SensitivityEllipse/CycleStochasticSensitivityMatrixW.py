import numpy as np

from model import Model
from model2 import Model2
from SensitivityEllipse.JacobiMatrixF import JacobiMatrixF


class CycleStochasticSensitivityMatrixW():
    def __init__(self):
        pass

    def w1(self, q1, q2, q3, b1, b2, b3, b4):
        """w1 находится из решения матричного уравнения W=BWB^T + Q, из предположения, что матрица Q - симметричная"""
        return -((q1 - b2 * b3 * q1 - b1 * b4 * q1 - (b4 ** 2) * q1 - b2 * b3 * (b4 ** 2) * q1 + b1 * (b4 ** 3) * q1 +
                  (b2 ** 2) * q2 - (b2 ** 3) * b3 * q2 + b1 * b2 ** 2 * b4 * q2 + 2 * b1 * b2 * q3 +
                  2 * b2 ** 2 * b3 * b4 * q3 - 2 * b1 * b2 * b4 ** 2 * q3) / (
                         -1 + b1 ** 2 + b2 * b3 + b1 ** 2 * b2 * b3 +
                         b2 ** 2 * b3 ** 2 - b2 ** 3 * b3 ** 3 + b1 * b4 - b1 ** 3 * b4 + 3 * b1 * b2 ** 2 * b3 ** 2 * b4 +
                         b4 ** 2 - b1 ** 2 * b4 ** 2 + b2 * b3 * b4 ** 2 - 3 * b1 ** 2 * b2 * b3 * b4 ** 2 - b1 * b4 ** 3 +
                         b1 ** 3 * b4 ** 3))

    def w2(self, q1, q2, q3, b1, b2, b3, b4):
        return -((-(b3 ** 2) * q1 + b2 * (b3 ** 3) * q1 - b1 * (b3 ** 2) * b4 * q1 - q2 + b1 ** 2 * q2
                  + b2 * b3 * q2 + b1 ** 2 * b2 * b3 * q2 + b1 * b4 * q2 - b1 ** 3 * b4 * q2
                  - 2 * b1 * b2 * b3 ** 2 * q3 - 2 * b3 * b4 * q3 + 2 * b1 ** 2 * b3 * b4 * q3) /
                  (1 - b1 ** 2 - b2 * b3 - b1 ** 2 * b2 * b3 - b2 ** 2 * b3 ** 2 + b2 ** 3 * b3 ** 3
                    - b1 * b4 + b1 ** 3 * b4 - 3 * b1 * b2 ** 2 * b3 ** 2 * b4 - b4 ** 2 + b1 ** 2 * b4 ** 2
                    - b2 * b3 * b4 ** 2 + 3 * b1 ** 2 * b2 * b3 * b4 ** 2 + b1 * b4 ** 3 - b1 ** 3 * b4 ** 3))

    def w3(self, q1, q2, q3, b1, b2, b3, b4):
        return -((-b1 * b3 * q1 - b2 * b3 ** 2 * b4 * q1 + b1 * b3 * b4 ** 2 * q1 - b1 * b2 ** 2 * b3 * q2 -
                  b2 * b4 * q2 + b1 ** 2 * b2 * b4 * q2 - q3 + b1 ** 2 * q3 + b2 ** 2 * b3 ** 2 * q3 + b4 ** 2 * q3 -
                  b1 ** 2 * b4 ** 2 * q3) / (
                         1 - b1 ** 2 - b2 * b3 - b1 ** 2 * b2 * b3 - b2 ** 2 * b3 ** 2 + b2 ** 3 * b3 ** 3 - b1 * b4 +
                         b1 ** 3 * b4 - 3 * b1 * b2 ** 2 * b3 ** 2 * b4 - b4 ** 2 + b1 ** 2 * b4 ** 2 - b2 * b3 * b4 ** 2 +
                         3 * b1 ** 2 * b2 * b3 * b4 ** 2 + b1 * b4 ** 3 - b1 ** 3 * b4 ** 3))

    def equals(self, x1, x2):
        return abs(x1 - x2) < 0.0001

    def get_cycle(self, model, x0, y0):
        x, y = model.transition_process(x0, y0, 50000-1)
        x_arr, y_arr = model.get_points(x, y, 1000)
        x_arr = list(map(lambda x: round(x, 8), x_arr))
        y_arr = list(map(lambda x: round(x, 8), y_arr))
        for i in range(1, len(x_arr)):
            if self.equals(x_arr[0], x_arr[i]) and self.equals(y_arr[0], y_arr[i]):
                x_arr = x_arr[0:i]
                y_arr = y_arr[0:i]
                break
        return x_arr, y_arr

    def get_B_matrix(self, jacobi_matrix_F_array):
        """Возвращает матрицу B = Fk * ... * F2 * F1"""
        B = jacobi_matrix_F_array[len(jacobi_matrix_F_array) - 1]
        for i in range(len(jacobi_matrix_F_array) - 2, -1, -1):
            B = np.dot(B, jacobi_matrix_F_array[i])
        return B

    def get_Q_matrix(self, model, jacobi_matrix_F_array):
        """Возвращает матрицу Q = Qk + Fk*Q(k-1)*Fk^T + ... + Fk*...*F2*Q1*F2^T*...*Fk^T
        Qk = [[1, 0], [0, 0]]"""
        length = len(jacobi_matrix_F_array)
        Q_array = model.get_Q_matrix() * length

        for i in range(length - 1):
            #   Q = [Q1, Q2, ... Qk]
            for j in range(i + 1, length):
                f_j = jacobi_matrix_F_array[j]
                f_j_transpose = np.transpose(jacobi_matrix_F_array[j])
                Q_array[i] = np.dot(f_j, Q_array[i])
                Q_array[i] = np.dot(Q_array[i], f_j_transpose)

        for i in range(1, length):
            Q_array[0][0][0] += Q_array[i][0][0]
            Q_array[0][0][1] += Q_array[i][0][1]
            Q_array[0][1][0] += Q_array[i][1][0]
            Q_array[0][1][1] += Q_array[i][1][1]
        return Q_array[0]


def main_model2():
    # построение W1 для тестового примерчика
    MU = 2.7
    SIGMA = 0.16
    X_0, Y_0 = 0.59, 0.69
    model = Model2(mu=MU, sigma=SIGMA)
    w_cycle_getter = CycleStochasticSensitivityMatrixW()
    x_cycle_arr, y_cycle_arr = w_cycle_getter.get_cycle(model, X_0, Y_0)
    jacobi_matrix_F_getter = JacobiMatrixF()
    jacobi_matrix_F_array = jacobi_matrix_F_getter.get_jacobi_matrix_F_array(model, x_cycle_arr, y_cycle_arr)

    print(len(x_cycle_arr))
    print('x =', x_cycle_arr[0])
    print('y =', y_cycle_arr[0])
    print('x_cycle_arr =', x_cycle_arr)
    print('y_cycle_arr =', y_cycle_arr)
    for i in range(len(jacobi_matrix_F_array)):
        print(f'f[{i}] =', jacobi_matrix_F_array[i])

    B = w_cycle_getter.get_B_matrix(jacobi_matrix_F_array)
    print('B =', B)
    Q = w_cycle_getter.get_Q_matrix(model, jacobi_matrix_F_array)
    print('Q =', Q)
    b1 = B[0][0]
    b2 = B[0][1]
    b3 = B[1][0]
    b4 = B[1][1]
    q1 = Q[0][0]
    q2 = Q[1][1]
    q3 = Q[0][1]
    w1 = w_cycle_getter.w1(q1, q2, q3, b1, b2, b3, b4)
    w2 = w_cycle_getter.w2(q1, q2, q3, b1, b2, b3, b4)
    w3 = w_cycle_getter.w3(q1, q2, q3, b1, b2, b3, b4)
    print(q1, q2, q3, b1, b2, b3, b4, sep=", ")
    print('w1 =', w1)
    print('w2 =', w2)
    print('w3 =', w3)

def main():
    MU = 6.
    I = 0.001
    X_0, Y_0 = 0.0, 0.0
    model = Model(_mu=MU, _I=I)
    w_cycle_getter = CycleStochasticSensitivityMatrixW()
    x_cycle_arr, y_cycle_arr = w_cycle_getter.get_cycle(model, X_0, Y_0)
    jacobi_matrix_F_getter = JacobiMatrixF()
    jacobi_matrix_F_array = jacobi_matrix_F_getter.get_jacobi_matrix_F_array(model, x_cycle_arr, y_cycle_arr)

    print(len(x_cycle_arr))
    print('x =', x_cycle_arr[0])
    print('y =', y_cycle_arr[0])
    print('x_cycle_arr =', x_cycle_arr)
    print('y_cycle_arr =', y_cycle_arr)
    for i in range(len(jacobi_matrix_F_array)):
        print(f'f[{i}] =', jacobi_matrix_F_array[i])

    B = w_cycle_getter.get_B_matrix(jacobi_matrix_F_array)
    print('B =', B)
    Q = w_cycle_getter.get_Q_matrix(model, jacobi_matrix_F_array)
    print('Q =', Q)
    b1 = B[0][0]
    b2 = B[0][1]
    b3 = B[1][0]
    b4 = B[1][1]
    q1 = Q[0][0]
    q2 = Q[1][1]
    q3 = Q[0][1]
    w1 = w_cycle_getter.w1(q1, q2, q3, b1, b2, b3, b4)
    w2 = w_cycle_getter.w2(q1, q2, q3, b1, b2, b3, b4)
    w3 = w_cycle_getter.w3(q1, q2, q3, b1, b2, b3, b4)
    print(q1, q2, q3, b1, b2, b3, b4, sep=", ")
    print('w1 =', w1)
    print('w2 =', w2)
    print('w3 =', w3)



if __name__ == '__main__':
    main_model2()
    #main()
