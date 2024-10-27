import math
import numpy as np

from Model import Model
from tqdm import tqdm


class LyapunovExponent:

    def __init__(self):
        pass

    def get_d(self, x: float, y: float) -> float:
        return math.sqrt(x ** 2 + y ** 2)

    def get_main_lyapunov_exponent(self, model: Model, x0: float, y0: float) -> float:
        """Алгоритм нахождения показателя Ляпунова находится в Theory/LyapunovExponentAlgorithm"""
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

    def get_lyapunov_exponent_data(self, model: Model, mu_start: float, mu_end: float, step: float,
                                   x0: float, y0: float, extend_attractor: bool = False):
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
