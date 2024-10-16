import math

import numpy as np

from tqdm import tqdm
from data_writer import DataWriter
from EquilibriumStochasticSensitivityMatrixW import StochasticSensitivityMatrixW
from EigValueVectors import get_eigenvalues, get_eigenvectors


class ConfidenceEllipse:
    """
    Класс для построения доверительного эллипса точки. При инициализации нужно передать параметры необходимые для
    построения эллипса - доверительную вероятность, интенсивность шума, собственные числа матрицы стохастической
    чувствительности W, значения собсвенных векторов собственных чисел, координаты точки
    """

    def __init__(self, q_2, eps, eta1, eta2, u1, u2, xeq, yeq):
        """
        :param q: константа, задающая доверительную вероятность эллипса
        :param eps: параметр интенсивности шума
        :param eta1: собственное число1 матрицы стохастической чувствительности W
        :param eta2: собственное число1 матрицы стохастической чувствительности W
        :param u1x: первая координата собственного вектора собственного числа1 матрицы стохастической чувствительности W
        :param u1y: вторая координата собственного вектора собственного числа1 матрицы стохастической чувствительности W
        :param u2x: первая координата собственного вектора собственного числа2 матрицы стохастической чувствительности W
        :param u2y: вторая координата собственного вектора собственного число2 матрицы стохастической чувствительности W
        :param xeq: координата x равновесия, для которого строится доверительный эллипс
        :param yeq: координата y равновесия, для которого строится доверительный эллипс
        """
        self.q_2 = q_2
        self.eps = eps
        self.eta1 = eta1
        self.eta2 = eta2
        self.u1x = u1[0]
        self.u1y = u1[1]
        self.u2x = u2[0]
        self.u2y = u2[1]
        self.xeq = xeq
        self.yeq = yeq

    def get_ro(self, fi):
        """
        :param fi: угловая координата, принадлежащая доверительному эллипсу
        :return: радиальную кордината доверительного эллипса для данного значение угловой координаты fi
        """
        ro2 = (2 * self.q_2 * (self.eps ** 2) * self.eta1 * self.eta2) / (self.eta2 * (math.cos(fi) ** 2) + self.eta1 * (math.sin(fi) ** 2))
        ro = math.sqrt(ro2)
        return ro, -ro

    def get_z1_z2(self, ro, fi):
        """
        :param ro: радиальная координата доверительного эллипса
        :param fi: угловая координата доверильного эллипса
        :return: декартовы координаты доверительного эллипса в системе z1, z2
        """
        z1 = ro * math.cos(fi)
        z2 = ro * math.sin(fi)
        return z1, z2

    def get_x_y(self, z1, z2):
        """
        :param z1: первая координата в системе координат z1, z2
        :param z2: вторая координата в системе координат z1, z2
        :return: возвращает значение координат в системе координат x, y
        """
        x = -((-self.u1y * self.u2x * self.xeq + self.u1x * self.u2y * self.xeq + self.u2y * z1 - self.u1y * z2)
              / (self.u1y * self.u2x - self.u1x * self.u2y))
        y = -((-self.u1y * self.u2x * self.yeq + self.u1x * self.u2y * self.yeq - self.u2x * z1 + self.u1x * z2)
              / (self.u1y * self.u2x - self.u1x * self.u2y))
        return x, y

    def get_confidence_ellipse(self):
        """
        Перебирает значения угла fi от 0 до 180 и находит координаты доверительного эллипса для уже известных
        собственных значений и собственных векторов матрицы стохастической чувствительности W
        :return: список значений координат x,y доверительного эллипса
        """
        xarr = []
        yarr = []
        for fi in tqdm(np.arange(0, 180, 0.01)):
            ro1, ro2 = self.get_ro(fi)
            z1_1, z1_2 = self.get_z1_z2(ro1, fi)
            z2_1, z2_2 = self.get_z1_z2(ro2, fi)
            x1, y1 = self.get_x_y(z1_1, z1_2)
            x2, y2 = self.get_x_y(z2_1, z2_2)
            xarr.append(x1)
            yarr.append(y1)
            xarr.append(x2)
            yarr.append(y2)
        return xarr, yarr


def main():
    #   Данный код для равновесия x,y находит матрицу стохастической чувствительности W, её собственные числа,
    #   собственные значения и строит доверительный эллипс
    # доверительная вероятность
    P = 0.995
    Q_2 = -math.log(1 - P)
    # интенсивность шума
    EPS = 0.0002
    # точка, для которой стоится эллипс
    X_EQ = 0.001063
    Y_EQ = 0.001063
    # значения матрицы стохастической чувствительности
    w1 = 1.007993557803417
    w2 = 1.007993557803417
    w3 = 0.08976331836818677
    eta1, eta2 = get_eigenvalues(w1, w2, w3)
    u1, u2 = get_eigenvectors(eta1, eta2, w1, w2, w3)
    confidence_ellipse = ConfidenceEllipse(q_2=Q_2, eps=EPS, eta1=eta1, eta2=eta2, u1=u1, u2=u2, xeq=X_EQ, yeq=Y_EQ)
    x_ellipse_arr, y_ellipse_arr = confidence_ellipse.get_confidence_ellipse()
    writer = DataWriter()
    path = f'../data/ellipse.txt'
    writer.write_data(x_ellipse_arr, y_ellipse_arr, path)


if __name__ == '__main__':
    main()
