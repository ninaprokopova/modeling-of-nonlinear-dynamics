from DataWriter import DataWriter
from Model import Model
from ModelIndicators.LyapunovExponent import LyapunovExponent


def get_lyapunov_exponent_for_point():
    # Скрипт для вычисления старшего показателя Ляпунова конкретного параметра скорости роста mu
    I = 0.00
    model = Model(I=I)
    model.mu = 2.5949
    x0, y0 = 1.1, 1.1
    lyapunov_exp_getter = LyapunovExponent()
    print('lp1: ', lyapunov_exp_getter.get_main_lyapunov_exponent(model, x0, y0))


def get_lyapunov_exponent_for_mu_arr():
    # Скрипт для вычисления старшего показателя Ляпунова для диапазона параметра скорости роста mu
    alpha = 1.5
    I = 0.00
    model = Model(alpha=alpha, mu=None, I=I)
    mu_start, mu_end, step = 0.13017, 0.13023, 0.000001
    x0, y0 = 1.1, 1.1
    lyapunov_exp_getter = LyapunovExponent()
    mu_arr, lyapunov_exponent_arr = lyapunov_exp_getter.get_lyapunov_exponent_data(model,
                                                                                   mu_start, mu_end, step,
                                                                                   x0, y0,
                                                                                   extend_attractor=False)
    data_writer = DataWriter()
    data_writer.write_data(mu_arr, lyapunov_exponent_arr, 'data/lp.txt')


if __name__ == '__main__':
    get_lyapunov_exponent_for_mu_arr()
    # get_lyapunov_exponent_for_point()


