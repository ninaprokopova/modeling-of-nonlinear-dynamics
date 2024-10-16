import random


class Model2:
    """Данная модель использовалась, чтобы отладить построение эллипсов для цикла"""

    def __init__(self, mu=None, sigma=None, eps=0):
        """
        Инициализация параметров модели
        """
        self.mu = mu
        self.sigma = sigma
        self.eps = eps

    def f(self, x, y):
        """
        Возвращает x_{t+1}=f(x_t, y_t)
        :param x: float, x_t
        :param y: float, y_t
        :return: float
        """
        return self.mu * x * (1 - x) + self.sigma * (y - x) + self.eps * random.normalvariate(0, 1)

    def g(self, x, y):
        """
        Возвращает y_{t+1}=g(x_t, y_t)
        :param x: float, x_t
        :return: float
        """
        return self.mu * y * (1 - y) + self.sigma * (x - y) + self.eps * random.normalvariate(0, 1)

    def get_dfdx(self, xi, yi):
        """Элемент с индексом 0,0 в матрице Якоби F"""
        f1 = self.mu * (1 - 2 * xi) - self.sigma
        return round(f1, 15)

    def get_dfdy(self, xi, yi):
        """Элемент с индексом 0,1 в матрице Якоби F"""
        f2 = self.sigma
        return round(f2, 15)

    def get_dgdx(self, xi, yi):
        """Элемент с индексом 1,0 в матрице Якоби F"""
        f2 = self.sigma
        return round(f2, 15)

    def get_dgdy(self, xi, yi):
        """Элемент с индексом 1,1 в матрице Якоби F"""
        f1 = self.mu * (1 - 2 * yi) - self.sigma
        return round(f1, 15)

    def get_Q_matrix(self):
        """Что такое матрица Q???"""
        return [[[1., 0.], [0., 1.]]]

    def get_Jacobi_matrix(self, xi, yi):
        return [[self.get_dfdx(xi, yi), self.get_dfdy(xi, yi)],
                [self.get_dgdx(xi, yi), self.get_dgdy(xi, yi)]]

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
            y = self.g(x_old, y_old)
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
            y = self.g(x_old, y_old)
        return x_arr, y_arr
