import numpy as np
import sys

from settings import Settings
from vars import Vars

class ModelParamsFormer(Settings, Vars):

    def __init__(self):
        self.form_clear_parameters()
        self.form_noize_graphs()
        self.form_measurment_noize_graphs()

    #### Getters =======================================================================================

    def get_clear_params(self):
        X = [self.dVe_clear, self.Fn_clear, self.Wdrn_clear]
        # return np.dot(self.H, X)[0]
        return X

    def get_noized_params(self):
        X = [self.dVe_noize, self.Fn_noize, self.Wdrn_noize]
        # return np.dot(self.H, X)[0]
        return X

    def get_measurment_noize_params(self):
        X = [self.dVe_mesuarment_noize, self.Fn_mesuarment_noize, self.Wdrn_mesuarment_noize]
        # return np.dot(self.H, X)[0]
        return X

    def get_egorushkin_params(self):
        X = []
        with open(self.egorushkin_params, 'r') as f:
            for v in f:
                X.append(float(v.strip()))
        return X

    #### Clear graphs ================================================================================== 
        
    # Чистые, незашумленные данные, их нужно сгенерировать 1 раз и сохранить
    # в файл
        
    # Если флаг 'generate_new_params' установлен (= True), то генерируется новая траектория, 
    # которая записывается в массивы и в файл clear_params.txt. Иначе, парсим файл clear_params.txt,
    # записываем данные отуда в массивы.
    def form_clear_parameters(self):
        if self.generate_new_params:
            # Generating new model parmas by rewriting 'clear_params.txt' file
            print('Generating new model params ...')
            with open(self.clear_params, 'w') as f:  
                self.dVe_clear.append(self.dVe_0)
                self.Fn_clear.append(self.Fn_0)
                self.Wdrn_clear.append(self.Wdrn_0)
                
                # Fill file witch initial conditions
                f.write("{} {} {}\n".format(self.dVe_0, self.Fn_0, self.Wdrn_0))
                for i in range(1, self.iter):
                    self.dVe_clear.append(self.dVe_clear[i - 1] - self.g * self.Fn_clear[i - 1])
                    self.Fn_clear.append(self.Fn_clear[i - 1] + self.dVe_clear[i - 1] / self.R + self.Wdrn_clear[i - 1])
                    self.Wdrn_clear.append(self.Wdrn_clear[i - 1])

                    f.write("{} {} {}\n".format(self.dVe_clear[i], self.Fn_clear[i], self.Wdrn_clear[i]))
            
            print('Model params generated successfully.')
        else:
            print('Getting model params from file ...')

            with open(self.clear_params, 'r') as f:
                params = f.read()
                if len(params) == 0:
                    print("Unable to read params from file: File is empty.")
                    print("Set 'generate_new_params' to 'True' in Gettings.py, to fill file.")
                    sys.exit(-1)
                
                lines = params.split('\n')
                for line in lines:
                    if line != '':
                        ve, fn, wdr = line.split()
                        self.dVe_clear.append(float(ve))
                        self.Fn_clear.append(float(fn))
                        self.Wdrn_clear.append(float(wdr))

            print('Model params getted successfully.')

    #### Noize graphs ==================================================================================

    # Зашумленные белым шумом графики. Флаг generate_new_params работает также, 
    # как и с чистыми графиками
    def form_noize_graphs(self):
        if self.generate_new_params:
            print('Generating new noized params ...')
            with open(self.noize_params, 'w') as f:        
                self.dVe_noize.append(self.dVe_0)
                self.Fn_noize.append(self.Fn_0)
                self.Wdrn_noize.append(self.Wdrn_0)
                
                f.write("{} {} {}\n".format(self.dVe_noize[0], self.Fn_noize[0], self.Wdrn_noize[0]))
                for i in range(1, self.iter):
                    self.X.append(np.dot(self.F , self.X[i - 1]) + np.dot(self.G, self.W[i - 1]))
                    self.dVe_noize.append(self.X[i][0][0])
                    self.Fn_noize.append(self.X[i][1][0])
                    self.Wdrn_noize.append(self.X[i][2][0])
                    
                    f.write("{} {} {}\n".format(self.dVe_noize[i], self.Fn_noize[i], self.Wdrn_noize[i]))
                
            print('Noized params generated successfully')
        else:
            print('Getting noized params from file ...')
            with open(self.noize_params, 'r') as f:
                params = f.read()
                if len(params) == 0:
                    print("Unable to read params from file: File is empty")
                    print("set 'generate_new_params' to 'True', to fill file")
                    sys.exit(-1)
                
                lines = params.split('\n')
                for line in lines:
                    if line != '':
                        ve, fn, wdr = line.split()
                        self.dVe_noize.append(float(ve))
                        self.Fn_noize.append(float(fn))
                        self.Wdrn_noize.append(float(wdr))
                    
            print('Noized params getted successfully')

    #### Measurement-noize graphs ======================================================================

    # Графикик с наложением измерительного шума. Фактически, это еще один белый шум, просто
    # меньшего порядка
    def form_measurment_noize_graphs(self):
        if self.generate_new_params:
            print('Generating new measurment-noized params ...')
            with open(self.meausure_noize_params, 'w') as f:        
                W_mesuarment_noize = np.random.normal(self.mu, self.sigma, self.iter)                
                for i in range(0, self.iter-1):

                    self.dVe_mesuarment_noize.append(self.dVe_noize[i]   + W_mesuarment_noize[i] * 1e-3)
                    self.Fn_mesuarment_noize.append(self.Fn_noize[i]     + W_mesuarment_noize[i] * 1e-5)
                    self.Wdrn_mesuarment_noize.append(self.Wdrn_noize[i] + W_mesuarment_noize[i] * 1e-8)
                
                    f.write("{} {} {}\n".format(self.dVe_mesuarment_noize[i], \
                        self.Fn_mesuarment_noize[i], \
                        self.Wdrn_mesuarment_noize[i]))
                
            print('Measurment-noized generated successfully')
        else:
            print('Getting noized params from file ...')
            with open(self.meausure_noize_params, 'r') as f:
                params = f.read()
                if len(params) == 0:
                    print("Unable to read params from file: File is empty")
                    print("set 'generate_new_params' to 'True', to fill file")
                    sys.exit(-1)
                
                lines = params.split('\n')
                for line in lines:
                    if line != '':
                        ve, fn, wdr = line.split()
                        self.dVe_mesuarment_noize.append(float(ve))
                        self.Fn_mesuarment_noize.append(float(fn))
                        self.Wdrn_mesuarment_noize.append(float(wdr))
                    
            print('Measurment-noized params getted successfully')
    