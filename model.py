import math
import random


class Model:

    def __init__(self, _alpha=1.5, _mu=None, _I=None, _eps=0):
        """
        Инициализация параметров модели
        :param _alpha: _alpha > 0 определяет силу Олли эффекта
        :param _mu: _mu > 0 характеризует скорость роста
        :param _I: определяем иммиграцию 0<=I<=1
        """
        self.alpha = _alpha
        self.mu = _mu
        self.I = _I
        self.eps = _eps

    def f(self, x, y):
        """Возвращает x_{t+1}=f(x_t, y_t)"""
        return (math.pow(x, self.alpha)) * (math.pow(math.e, self.mu * (1 - y))) + self.I \
               + self.eps * random.normalvariate(0, 1)

    def g(self, x):
        """Возвращает y_{t+1}=g(x_t, y_t)"""
        return x

    def get_dfdx(self, xi, yi):
        """Элемент с индексом 0,0 в матрице Якоби F"""
        f1 = 1.5 * math.sqrt(xi) * (math.e ** (self.mu * (1 - yi)))
        return round(f1, 7)

    def get_dfdy(self, xi, yi):
        """Элемент с индексом 0,1 в матрице Якоби F"""
        f2 = - self.mu * xi * math.sqrt(xi) * (math.e ** (self.mu * (1 - yi)))
        return round(f2, 7)

    def get_dgdx(self, xi, yi):
        """Элемент с индексом 1,0 в матрице Якоби F"""
        return 1.

    def get_dgdy(self, xi, yi):
        """Элемент с индексом 1,1 в матрице Якоби F"""
        return 0.

    def get_Q_matrix(self):
        """Что такое матрица Q??? - матрица из уравнений всяких из статьи Льва Борисовича и Ирины Адольфовны"""
        return [[[1., 0.], [0., 0.]]]

    def transition_process(self, x_start, y_start, iterations_n=1000):
        """
        Переходный процесс
        :param x_start: float, начальное x
        :param y_start: float, начальное y
        :param iterations_n: int, количество иттераций переходного процесса
        :return: float, float
        """
        x = x_start
        y = y_start
        eps = self.eps
        self.eps = 0
        for _ in range(iterations_n):
            x_old = x
            y_old = y
            x = self.f(x_old, y_old)
            y = self.g(x_old)
        self.eps = eps
        return x, y

    def get_points(self, x0, y0, iterations_n):
        """
        Возвращает массивы с точками x и y полученные за количество иттераций iterations_n
        :param x0: float
        :param y0: float
        :param iterations_n: int
        :return: [float], [float]
        """
        x = x0
        y = y0
        x_arr = []
        y_arr = []
        for _ in range(iterations_n):
            x_arr.append(x)
            y_arr.append(y)
            if x < 0:
                break
            x_old = x
            y_old = y
            x = self.f(x_old, y_old)
            y = self.g(x_old)
        return x_arr, y_arr
