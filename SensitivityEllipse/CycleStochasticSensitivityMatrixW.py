import numpy as np

from Model import Model


class CycleStochasticSensitivityMatrixW():
    # Класс для нахождения матриц стохастической чувствительности цикла

    def __init__(self):
        pass

    # w1, w2, w3 - элементы матрицы стохастической чувствительности W_1 = [[w1, w3], [w3, w2]]
    # формулы для нахождения w1, w2, w3 получены в ВольфрамМатематика
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
        """Возвращает матрицу Q = Qk + Fk*Q(k-1)*Fk^T + ... + Fk*...*F2*Q1*F2^T*...*Fk^T"""
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

    def get_W_matrix(self, model: Model, previous_W_matrix: list[list[float]],
                     x: float, y: float) -> (float, float, float):
        """Возвращает элементы матрицы стохастической чувствительности w1, w2, w3. Матрица вычислется рекурсивно
        по формуле: w_i = jacobi_m_(i-1) * w_(i-1) * jacobi_m_(i-1)^T + q_i"""
        jacobi_matrix = model.get_Jacobi_matrix(x, y)
        jacobi_matrix_transpose = np.transpose(jacobi_matrix)
        q_matrix = model.get_Q_matrix()
        w_i_matrix = np.dot(jacobi_matrix, previous_W_matrix)
        w_i_matrix = np.dot(w_i_matrix, jacobi_matrix_transpose)
        w_i_matrix = np.add(w_i_matrix, q_matrix)
        return w_i_matrix[0][0][0], w_i_matrix[0][1][1], w_i_matrix[0][0][1]
