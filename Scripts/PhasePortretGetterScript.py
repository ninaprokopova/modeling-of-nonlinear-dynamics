from Model import Model
from Model2 import Model2
from DataWriter import DataWriter


def main():
    """
    Скрипт для построения точек траектории после переходного процесса.
    По этим точкам можно построить фазовый портрет или временной ряд.
    """
    alpha = 1.5
    I = 0.00
    mu = 0.131
    x0, y0 = 1.1, 1.1
    points_count = 300000
    start_iterations_count = 3000
    eps = 0.001
    model = Model(_alpha=alpha, _mu=mu, _I=I, _eps=eps)
    x, y = model.transition_process(x0, y0, start_iterations_count)
    x_arr, y_arr = model.get_points(x, y, points_count)
    w = DataWriter()
    w.write_data(x_arr, y_arr, "data/phase_p.txt")


def main_model2():
    mu = 2.7
    sigma = 0.16
    x0, y0 = 0.59, 0.69
    points_count = 30000
    start_iterations_count = 10000
    eps = 0.006
    model = Model2(mu=mu, sigma=sigma, eps=eps)
    x, y = model.transition_process(x0, y0, start_iterations_count)
    x_arr, y_arr = model.get_points(x, y, points_count)
    w = DataWriter()
    w.write_data(x_arr, y_arr, "data/phase_p.txt")


if __name__ == '__main__':
    main()