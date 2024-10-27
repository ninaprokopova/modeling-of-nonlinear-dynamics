import numpy as np

from tqdm import tqdm
from Model import Model
from DataWriter import DataWriter


class RotationNumber:
    """Класс для определения числа вращения системы при конкретном mu"""

    def get_rotation_number(self, x_arr):
        x_max = max(x_arr)
        x_min = min(x_arr)
        x_medium = (x_max + x_min) / 2
        counter = 0
        for i in range(len(x_arr) - 1):
            if x_arr[i] < x_medium <= x_arr[i + 1]:
                counter += 1
        return len(x_arr) / counter

    def get_rotation_number_array(self, model, mu_start, mu_end, step, x0, y0):
        mu_arr = []
        rot_num_arr = []
        for mu in tqdm(np.arange(mu_start, mu_end, step)):
            model.mu = mu
            x, y = model.transition_process(x0, y0, iterations_n=int(1e3))
            x_arr, _ = model.get_points(x, y, iterations_n=int(2e5))
            rot_num = self.get_rotation_number(x_arr)
            mu_arr.append(mu)
            rot_num_arr.append(rot_num)
        return mu_arr, rot_num_arr
