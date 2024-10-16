import math

from SensitivityEllipse.EigenValueVectors import get_eigenvalues, get_eigenvectors
from SensitivityEllipse.ConfidenceEllipse import ConfidenceEllipse
from data_writer import DataWriter


def main():
    #   Данный код для точки x,y и матрицы стохастической чувствительности W
    #   находит собственные числа, собственные значения и строит доверительный эллипс

    # доверительная вероятность
    P = 0.995
    Q_2 = -math.log(1 - P)
    # интенсивность шума
    EPS = 0.006
    # Следующие 6 строчек нужно получить из скрипта CycleStohasticSensitivityMatrixWScript
    # Точка, для которой нужно построить эллипс, и её матрица стохастической чувствительности W
    X = 0.67407407
    Y = 0.57777778
    w1 = 5.5170090195675145
    w2 = 9.907257038351842
    w3 = -5.61946207692496
    eta1, eta2 = get_eigenvalues(w1, w2, w3)
    u1, u2 = get_eigenvectors(eta1, eta2, w1, w2, w3)
    confidence_ellipse = ConfidenceEllipse(q_2=Q_2, eps=EPS, eta1=eta1, eta2=eta2, u1=u1, u2=u2, xeq=X, yeq=Y)
    x_ellipse_arr, y_ellipse_arr = confidence_ellipse.get_confidence_ellipse()
    writer = DataWriter()
    path = f'../data/ellipse.txt'
    writer.write_data(x_ellipse_arr, y_ellipse_arr, path)


if __name__ == '__main__':
    main()