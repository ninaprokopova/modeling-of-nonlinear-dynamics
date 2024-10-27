from Model import Model
from Model2 import Model2
from SensitivityEllipse.JacobiMatrixF import JacobiMatrixF
from SensitivityEllipse.CycleStochasticSensitivityMatrixW import CycleStochasticSensitivityMatrixW
from DataWriter import DataWriter


def main_model2():
    # построение W1 для тестового примера, чтобы убедиться, что код работает верно
    MU = 2.7
    SIGMA = 0.16
    X_0, Y_0 = 0.59, 0.69
    model = Model2(mu=MU, sigma=SIGMA)
    w_cycle_getter = CycleStochasticSensitivityMatrixW()
    x_cycle_arr, y_cycle_arr = w_cycle_getter.get_cycle(model, X_0, Y_0)
    jacobi_matrix_F_getter = JacobiMatrixF()
    jacobi_matrix_F_array = jacobi_matrix_F_getter.get_jacobi_matrix_F_array(model, x_cycle_arr, y_cycle_arr)

    print(f"Длина цикла = {len(x_cycle_arr)}")
    print(f"x_cycle_arr = {x_cycle_arr}")
    print(f"y_cycle_arr = {y_cycle_arr}")
    for i in range(len(jacobi_matrix_F_array)):
        print(f"f[{i}] = {jacobi_matrix_F_array[i]}")

    B = w_cycle_getter.get_B_matrix(jacobi_matrix_F_array)
    print(f'B = {B}')
    Q = w_cycle_getter.get_Q_matrix(model, jacobi_matrix_F_array)
    print(f'Q = {Q}')
    b1 = B[0][0]
    b2 = B[0][1]
    b3 = B[1][0]
    b4 = B[1][1]
    q1 = Q[0][0]
    q2 = Q[1][1]
    q3 = Q[0][1]
    w1 = w_cycle_getter.w1(q1, q2, q3, b1, b2, b3, b4)
    w2 = w_cycle_getter.w2(q1, q2, q3, b1, b2, b3, b4)
    w3 = w_cycle_getter.w3(q1, q2, q3, b1, b2, b3, b4)
    print("# Точка, для которой нужно построить эллипс, и её матрица стохастической чувствительности W")
    print(f"X = {x_cycle_arr[0]}")
    print(f"Y = {y_cycle_arr[0]}")
    print(f"w1 = {w1}")
    print(f"w2 = {w2}")
    print(f"w3 = {w3}")

    w_1 = [[w1, w3], [w3, w2]]
    x = x_cycle_arr[0]
    y = y_cycle_arr[0]
    w_2 = w_cycle_getter.get_W_matrix(model, w_1, x, y)
    print("# Точка, для которой нужно построить эллипс, и её матрица стохастической чувствительности W")
    print(f"X = {x_cycle_arr[1]}")
    print(f"Y = {y_cycle_arr[1]}")
    print(f"w1 = {w_2[0][0][0]}")
    print(f"w2 = {w_2[0][1][1]}")
    print(f"w3 = {w_2[0][0][1]}")

def main():
    """
    Скрипт находит матрицы стохастической чувствительности для циклов:
    1. Находит цикл
    2. Находит матрицы Якоби для элементов цикла
    3. Находит матрицы B и Q todo: стоит сделать ссылку на теорию
    4. Находит матрицу W_1
    5. Находит матрицы W_2..W_n, где n - длина цикла
    6. Записывает в файл точку и элементы матрицы W в формате (x y w1 w2 w3)
    """

    #  Параметры, которые нужно задать перед запуском скрипта
    MU = 0.131
    I = 0.00
    X_0, Y_0 = 1.1, 1.1
    model = Model(mu=MU, I=I)
    w_cycle_getter = CycleStochasticSensitivityMatrixW()
    x_cycle_arr, y_cycle_arr = w_cycle_getter.get_cycle(model, X_0, Y_0)
    jacobi_matrix_F_getter = JacobiMatrixF()
    jacobi_matrix_F_array = jacobi_matrix_F_getter.get_jacobi_matrix_F_array(model, x_cycle_arr, y_cycle_arr)

    print(f"Длина цикла = {len(x_cycle_arr)}")
    print(f"x_cycle_arr = {x_cycle_arr}")
    print(f"y_cycle_arr = {y_cycle_arr}")
    for i in range(len(jacobi_matrix_F_array)):
        print(f"f[{i}] = {jacobi_matrix_F_array[i]}")

    B = w_cycle_getter.get_B_matrix(jacobi_matrix_F_array)
    print(f'B = {B}')
    Q = w_cycle_getter.get_Q_matrix(model, jacobi_matrix_F_array)
    print(f'Q = {Q}')
    b1 = B[0][0]
    b2 = B[0][1]
    b3 = B[1][0]
    b4 = B[1][1]
    q1 = Q[0][0]
    q2 = Q[1][1]
    q3 = Q[0][1]
    w1 = w_cycle_getter.w1(q1, q2, q3, b1, b2, b3, b4)
    w2 = w_cycle_getter.w2(q1, q2, q3, b1, b2, b3, b4)
    w3 = w_cycle_getter.w3(q1, q2, q3, b1, b2, b3, b4)
    print("# Точка, для которой нужно построить эллипс, и её матрица стохастической чувствительности W")
    print(f"X = {x_cycle_arr[0]}")
    print(f"Y = {y_cycle_arr[0]}")
    print(f"w1 = {w1}")
    print(f"w2 = {w2}")
    print(f"w3 = {w3}")

    w_i = [[w1, w3], [w3, w2]]

    x_arr = [x_cycle_arr[0]]
    y_arr = [y_cycle_arr[0]]
    w1_arr = [w1]
    w2_arr = [w2]
    w3_arr = [w3]
    for i in range(len(x_cycle_arr) - 1):
        x = x_cycle_arr[i]
        y = y_cycle_arr[i]
        w1, w2, w3 = w_cycle_getter.get_W_matrix(model, w_i, x, y)
        print(f"{i+2} элемент цикла")
        print("# Точка, для которой нужно построить эллипс, и её матрица стохастической чувствительности W")
        print(f"X = {x_cycle_arr[i+1]}")
        print(f"Y = {y_cycle_arr[i+1]}")
        print(f"w1 = {w1}")
        print(f"w2 = {w2}")
        print(f"w3 = {w3}")
        w_i = [[w1, w3], [w3, w2]]
        x_arr.append(x_cycle_arr[i+1])
        y_arr.append(y_cycle_arr[i+1])
        w1_arr.append(w1)
        w2_arr.append(w2)
        w3_arr.append(w3)

    data_writer = DataWriter()
    data_writer.write_data_5arr(x_arr, y_arr, w1_arr, w2_arr, w3_arr, "..\data\w_cycle.txt")


if __name__ == '__main__':
    #main_model2()
    main()
