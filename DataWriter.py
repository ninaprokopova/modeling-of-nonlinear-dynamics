import os


class DataWriter:
    def __init__(self):
        None

    def write_data(self, x, y, file_path):
        float_precision = 6
        with open(file_path, 'w') as f:
            for i in range(len(x)):
                f.write(f'{x[i]:.{float_precision}f} {y[i]:.{float_precision}f}\n')
        print(f'the data is written to a file {file_path}')

    def write_data_3arr(self, x, y, z, file_path):
        float_precision = 6
        with open(file_path, 'w') as f:
            for i in range(len(x)):
                x_el = f'{x[i]:.{float_precision}f}'
                y_el = f'{y[i]:.{float_precision}f}'
                z_el = f'{z[i]:.{float_precision}f}'
                f.write(f'{x_el} {y_el} {z_el}\n')
        print(f'the data is written to a file {file_path}')

    def write_data_5arr(self, x, y, w1, w2, w3, file_path):
        float_precision = 6
        with open(file_path, 'w') as f:
            for i in range(len(x)):
                x_el = f'{x[i]:.{float_precision}f}'
                y_el = f'{y[i]:.{float_precision}f}'
                w1_el = f'{w1[i]:.{float_precision}f}'
                w2_el = f'{w2[i]:.{float_precision}f}'
                w3_el = f'{w3[i]:.{float_precision}f}'
                f.write(f'{x_el} {y_el} {w1_el} {w2_el} {w3_el}\n')
        print(f'the data is written to a file {file_path}')

    def clean_dir(self):
        dir = 'data'
        for f in os.listdir(dir):
            file_name, file_extension = os.path.splitext(f)
            if file_extension == ".txt":
                os.remove(os.path.join(dir, f))
