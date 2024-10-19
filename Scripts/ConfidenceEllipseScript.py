import math

from SensitivityEllipse.EigenValueVectors import get_eigenvalues, get_eigenvectors
from SensitivityEllipse.ConfidenceEllipse import ConfidenceEllipse
from DataWriter import DataWriter


def main():
    #   Данный код для точки x,y и матрицы стохастической чувствительности W
    #   находит собственные числа, собственные значения и строит доверительный эллипс

    # доверительная вероятность
    P = 0.995
    Q_2 = -math.log(1 - P)
    # интенсивность шума
    EPS = 0.0001
    # Следующие 6 строчек нужно получить из скрипта CycleStohasticSensitivityMatrixWScript
    # Точка, для которой нужно построить эллипс, и её матрица стохастической чувствительности W
    X = 1.99336832
    Y = 7.63349383
    w1 = 153814.6850562113
    w2 = 4013601.1075127176
    w3 = 785711.8845088588
    eta1, eta2 = get_eigenvalues(w1, w2, w3)
    u1, u2 = get_eigenvectors(eta1, eta2, w1, w2, w3)
    confidence_ellipse = ConfidenceEllipse(q_2=Q_2, eps=EPS, eta1=eta1, eta2=eta2, u1=u1, u2=u2, xeq=X, yeq=Y)
    x_ellipse_arr, y_ellipse_arr = confidence_ellipse.get_confidence_ellipse()
    writer = DataWriter()
    path = f'../data/ellipse.txt'
    writer.write_data(x_ellipse_arr, y_ellipse_arr, path)

def main_many_ellipses():
    #   Данный код для точек из файла data\w_cycle.txt и  их матриц стохастической чувствительности W
    #   находит собственные числа, собственные значения и строит доверительные эллипсы
    #   Файл data\w_cycle.txt можно получить с помощью скрипта CycleStohasticSensitivityMatrxiWScript,
    #   вызвав функцию main()

    # доверительная вероятность
    P = 0.995
    Q_2 = -math.log(1 - P)
    # интенсивность шума
    EPS = 0.0005

    x_ellipses_arr = []
    y_ellipses_arr = []

    with open ("..\data\w_cycle.txt") as f:
        for line in f:
            X, Y, w1, w2, w3 = list(map(float, line.split()))
            eta1, eta2 = get_eigenvalues(w1, w2, w3)
            u1, u2 = get_eigenvectors(eta1, eta2, w1, w2, w3)
            confidence_ellipse = ConfidenceEllipse(q_2=Q_2, eps=EPS, eta1=eta1, eta2=eta2, u1=u1, u2=u2, xeq=X, yeq=Y)
            x_ellipse_arr, y_ellipse_arr = confidence_ellipse.get_confidence_ellipse()
            x_ellipses_arr += x_ellipse_arr
            y_ellipses_arr += y_ellipse_arr
    writer = DataWriter()
    path = f'../data/ellipse.txt'
    writer.write_data(x_ellipses_arr, y_ellipses_arr, path)


if __name__ == '__main__':
    #main()
    main_many_ellipses()