import math

from bifd import Bif_d
from model import Model
from data_writer import DataWriter
from EigValueVectors import get_eigenvalues, get_eigenvectors

class StochasticSensitivityMatrixW:
    """
    Матрица стохастической чувствительности W = ( w1  w3 )
                                                ( w3  w2 ) является симметрической.
    Класс содержит функции для нахождения элементов w1, w2, w3 матрицы стохастической чувствительности W,
    для нахождения собсвенных чисел и собственных векторов этой матрицы.
    """

    def __init__(self):
        pass

    def w1(self, mu: float, y: float, x: float) -> float:
        """Вычисляет значение w1 для равновесия (x,y)"""
        E = math.e
        return (1. + E**(mu * (1 - y)) * mu * x**1.5) / \
               (1.
                - 2.25 * E**(2. * mu * (1. - y)) * x**1.
                + E**(mu * (1. - y)) * mu * x**1.5
                + 2.25 * E**(3. * mu * (1. - y)) * mu * x**2.5
                - E**(2. * mu * (1. - y)) * mu**2 * x**3.
                - E**(3. * mu * (1. - y)) * mu**3 * x**4.5)

    def w2(self, mu: float, y: float, x: float) -> float:
        """Вычисляет значение w2 для равновесия (x,y)"""
        E = math.e
        return (1. + E**(mu * (1 - y)) * mu * x**1.5) / \
               (1.
                - 2.25 * E**(2. * mu * (1. - y)) * x**1.
                + E**(mu * (1. - y)) * mu * x**1.5
                + 2.25 * E**(3. * mu * (1. - y)) * mu * x**2.5
                - E**(2. * mu * (1. - y)) * mu**2 * x**3.
                - E**(3. * mu * (1. - y)) * mu**3 * x**4.5)

    def w3(self, mu: float, y: float, x: float) -> float:
        """Вычисляет значение w3 для равновесия (x,y)"""
        E = math.e
        return (1.5 * E**(mu * (1 - y)) * x**0.5) / \
               (1.
                - 2.25 * E**(2. * mu * (1. - y)) * x**1.
                + E**(mu * (1. - y)) * mu * x**1.5
                + 2.25 * E**(3. * mu * (1. - y)) * mu * x**2.5
                - E**(2. * mu * (1. - y)) * mu**2 * x**3.
                - E**(3. * mu * (1. - y)) * mu**3 * x**4.5)

    def get_w1_arr(self, mu_arr: list[float], x_arr: list[float], y_arr: list[float]) -> list[float]:
        """Формирует списко значений w1 для списка равновесий и параметра скорости роста"""
        if len(mu_arr) != len(x_arr) and len(mu_arr) != len(x_arr):
            raise Exception("Списки mu_arr, x_arr, y_arr имеют разную длину")
        w1_arr = []
        for i in range(len(mu_arr)):
            mu, x, y = mu_arr[i], x_arr[i], y_arr[i]
            w1 = self.w1(mu, x, y)
            w1_arr.append(w1)
        return w1_arr

    def get_w2_array(self, mu_arr: list[float], x_arr: list[float], y_arr: list[float]) -> list[float]:
        """Формирует список значений w2 для списка равновесий и параметра скорости роста"""
        if len(mu_arr) != len(x_arr) and len(mu_arr) != len(x_arr):
            raise Exception("Списки mu_arr, x_arr, y_arr имеют разную длину")
        w2_arr = []
        for i in range(len(mu_arr)):
            mu, x, y = mu_arr[i], x_arr[i], y_arr[i]
            w2 = self.w2(mu, x, y)
            w2_arr.append(w2)
        return w2_arr

    def get_w3_array(self, mu_arr: list[float], x_arr: list[float], y_arr: list[float]) -> list[float]:
        """Формирует список значений w3 для списка равновесий и параметра скорости роста"""
        if len(mu_arr) != len(x_arr) and len(mu_arr) != len(x_arr):
            raise Exception("Списки mu_arr, x_arr, y_arr имеют разную длину")
        w3_arr = []
        for i in range(len(mu_arr)):
            mu, x, y = mu_arr[i], x_arr[i], y_arr[i]
            w3 = self.w3(mu, x, y)
            w3_arr.append(w3)
        return w3_arr


def main():
    #    Данный код находит значения матрицы стохастической чувствительности w1, w2, w3
    #    для диапазона значений параметра скорости роста mu

    #mu_start, mu_end, step = 0.501, 0.995, 1e-3
    mu_start, mu_end, step = 0.204, 0.995, 1e-3
    start_iterations_n, iter_count = 2000, 1
    draw_plot = True
    x_start, y_start = 1.1, 1.1
    model = Model(_alpha=1.5, _I=0)
    bifd = Bif_d(model)
    mu_arr, x_arr, y_arr = bifd.get_data_stretching_attractor(mu_start, mu_end, step, start_iterations_n, iter_count,
                                                              x_start, y_start, draw_plot)
    w = StochasticSensitivityMatrixW()
    w1_arr = w.get_w1_arr(mu_arr, x_arr, y_arr)
    w2_arr = w.get_w2_array(mu_arr, x_arr, y_arr)
    w3_arr = w.get_w3_array(mu_arr, x_arr, y_arr)

    dw = DataWriter()
    dw.write_data(mu_arr, w1_arr, f'data\\w1.txt')
    dw.write_data(mu_arr, w2_arr, f'data\\w2.txt')
    dw.write_data(mu_arr, w3_arr, f'data\\w3.txt')



def main2():
    w = StochasticSensitivityMatrixW()
    mu = 0.6
    x, y = 0.001063, 0.001063
    w1 = w.w1(mu, x, y)
    w2 = w.w2(mu, x, y)
    w3 = w.w3(mu, x, y)
    print(f'w1 = {w1}')
    print(f'w2 = {w2}')
    print(f'w3 = {w3}')
    eta1, eta2 = get_eigenvalues(w1, w2, w3)
    print(eta1, eta2)
    u1, u2 = get_eigenvectors(eta1, eta2, w1, w2, w3)
    print(u1, u2)

if __name__ == '__main__':
    main2()
