import math


class JacobiMatrixF:

    def get_jacobi_matrix_F(self, model, x, y):
        """Возвращает матрицу Якоби F для элемента цикла x, y и модели model (тут хранятся параметры)"""
        f1 = model.get_dfdx(x, y)
        f2 = model.get_dfdy(x, y)
        f3 = model.get_dgdx(x, y)
        f4 = model.get_dgdy(x, y)
        jacobi_matrix_f = [[f1, f2], [f3, f4]]
        return jacobi_matrix_f

    def get_jacobi_matrix_F_array(self, model, x_arr, y_arr):
        """Возвращает список из матриц Якоби для списка равновесий"""
        jacobi_matrix_F_array = []
        for i in range(len(x_arr)):
            jacobi_matrix_F = self.get_jacobi_matrix_F(model, x_arr[i], y_arr[i])
            jacobi_matrix_F_array.append(jacobi_matrix_F)
        return jacobi_matrix_F_array