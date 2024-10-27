import math
import numpy as np
import Model

from tqdm import tqdm
from DataWriter import DataWriter


class LyapunovExponent:
    # todo сделать описание функций и типы

    def __init__(self):
        pass

    def get_d(self, x, y):
        return math.sqrt(x ** 2 + y ** 2)

    # todo нужно приложить ссылку на алгоритм
    def get_main_lyapunov_exponent(self, model: Model.Model, x0, y0):
        start_iterations_n = 10 ** 3
        x0, y0 = model.transition_process(x0, y0, start_iterations_n)
        # print('mu: ', model.mu, 'x0, y0: ', x0, y0)
        iterations_n = 10 ** 6
        x_arr, y_arr = model.get_points(x0, y0, iterations_n)
        ln_p_arr = []

        DELTA = 10 ** (-6)
        R = DELTA / (math.sqrt(2))

        xvi, yvi = x_arr[0] + R, y_arr[0] + R
        for i in range(1, len(x_arr)):
            xvi_ = model.f(xvi, yvi)
            yvi_ = model.g(xvi)

            di = (xvi_ - x_arr[i], yvi_ - y_arr[i])
            di_len = self.get_d(di[0], di[1])

            if di_len == 0:
                di_len = 1e-20

            p_i = di_len / DELTA

            ln_p_arr.append(math.log(p_i))

            xvi = x_arr[i] + di[0] / di_len * DELTA
            yvi = y_arr[i] + di[1] / di_len * DELTA

        return round(np.sum(ln_p_arr) / len(ln_p_arr), 20)

    def get_lyapunov_exponent_data(self, model, mu_start, mu_end, step, x0, y0, extend_attractor=False):
        lyapunov_exponent_arr = []
        mu_arr = []

        range_ = np.arange(mu_start, mu_end, step)
        x_start = []
        y_start = []
        if not extend_attractor:
            x_start = np.ones(len(range_)) * x0
            y_start = np.ones(len(range_)) * y0
        else:
            x = x0
            y = y0
            for mu_ in range_:
                model.mu = round(mu_, 20)
                x, y = model.transition_process(x, y, 5 * (10 ** 4))
                x_start.append(x)
                y_start.append(y)

        counter = 0
        for mu_ in tqdm(range_):
            model.mu = round(mu_, 20)
            print(model.mu)
            x_ = x_start[counter]
            y_ = y_start[counter]
            x_, y_ = model.transition_process(x_, y_, 10000)
            l = self.get_main_lyapunov_exponent(model, x_, y_)
            counter += 1
            lyapunov_exponent_arr.append(l)
            mu_arr.append(mu_)
        return mu_arr, lyapunov_exponent_arr


# todo вынести код ниже в отдельный скрипт
if __name__ == '__main__':
    lyapunov_exp_getter = LyapunovExponent()
    alpha = 1.5
    I = 0.00
    model = Model.Model(_alpha=alpha, _mu=None, _I=I)
    mu_start, mu_end, step = 0.13017, 0.13023, 0.000001
    x0, y0 = 1.1, 1.1
    mu_arr, lyapunov_exponent_arr = lyapunov_exp_getter.get_lyapunov_exponent_data(model, mu_start, mu_end, step, x0, y0,
                                                                                   extend_attractor=False)
    data_writer = DataWriter()
    data_writer.write_data(mu_arr, lyapunov_exponent_arr, 'data/lp.txt')

    # model.mu = 2.5949
    # x0, y0 = 1.1, 1.1
    # print('lp1: ', lyapunov_p.get_main_lyapunov_exponent(model, x0, y0))
