import math


def get_eigenvalues(w1, w2, w3):
    """Возвращает собственные числа эта1 и эта2 матрицы стохастической чувствительности W"""
    eta1 = (w1 + w2 + math.sqrt((w1 - w2) ** 2 + (2 * w3) ** 2)) / 2
    eta2 = (w1 + w2 - math.sqrt((w1 - w2) ** 2 + (2 * w3) ** 2)) / 2

    print(f'Eigenvalue_1 = {eta1}')
    print(f'Eigenvalue_2 = {eta2}')
    return eta1, eta2


def _get_u_value(w, eta):
    """
    Возвращает значение компоненты собственного вектора
    в особенной ситуации, описанной в теории
    """
    if abs(w - eta) < 0.00001:
        return 1
    return 0


def get_eigenvectors(eta1, eta2, w1, w2, w3):
    """
    Возвращает значения собственных векторов.
    При w3=0 возникает особенный случай,который обрабатывается отдельно
    """
    if abs(w3) < 0.00001:
        u11 = _get_u_value(w1, eta1)
        u12 = _get_u_value(w2, eta1)
        u21 = _get_u_value(w1, eta2)
        u22 = _get_u_value(w2, eta2)
        u1 = (u11, u12)
        u2 = (u21, u22)
    else:
        u1 = ((eta1 - w2) / w3, 1.)
        u2 = ((eta2 - w2) / w3, 1.)
    u1_normalized = get_normalized_vector(u1)
    u2_normalized = get_normalized_vector(u2)

    print(f'Eigenvalue_1 = {eta1}, eigenvector_1 = {u1_normalized}')
    print(f'Eigenvalue_2 = {eta2}, eigenvector_2 = {u2_normalized}')
    return u1_normalized, u2_normalized


def get_normalized_vector(v: tuple[float, float]) -> tuple[float, float]:
    v_length = math.sqrt(v[0] ** 2 + v[1] ** 2)
    v0 = v[0] / v_length
    v1 = v[1] / v_length
    return (v0, v1)