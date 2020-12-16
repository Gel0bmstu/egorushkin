import numpy as np
import math
import random

from vars import Vars
from settings import Settings

class FiltersPool(Settings, Vars):
    
    def __init__(self):
        pass

    def Q_picker(self, params):

        print('Try to pick best Q ...')

        saf_min = math.inf
        q_min = 0
        Q = 1e-12
        chosen_idx = 0

        X = [self.dVe_mesuarment_noize, self.Fn_mesuarment_noize, self.Wdrn_mesuarment_noize]

        for i in range(15):
            self.Q_peaker_legend.append('Q = {}'.format(Q))
            self.Q_line_width.append(0.5)
            self.Q_line_style.append('-')

            self.X_Q_picker.append([])
            saf = self.simple_kalman(params, Q)

            if saf < saf_min:
                q_min = Q
                saf_min = saf
                chosen_idx = i

            for v in self.X_simple_kalman:
                self.X_Q_picker[i].append(v)

            Q = Q / 1.5

        self.X_Q_picker.append([])
        for v in X[self.idx]:
            self.X_Q_picker[len(self.X_Q_picker)-1].append(v)

        self.Q_peaker_legend.append('Noized values')
        self.Q_line_width.append(1.5)
        self.Q_line_style.append('-.')
        self.Q_line_width[chosen_idx] = 1.5

        self.Q = q_min

        print('Picked Q is: {}'.format(q_min))
        print('Min SAF is: {}'.format(saf_min))

    # Обычный Калман
    def simple_kalman(self, noized_params, Q=None):
        """ Проверка фильтров на сходимость

        @param: noized_params [np.array([dVe, Fn, Wdr])]
        """
        print('Simple kalman filter ...')

        if Q is None:
            Q = self.Q


        H = np.array([[1, 0, 0]])
        # self.H = H

        # Выбираем фильтруемое значение с помощью вектора состояний H
        # FIXME: Тут приколы с numpy, т.к. H = [[1, 0, 0]], то результирующее Z
        # будет равно [[data, data, data ... ]], то есть обернуто в еще одни [].
        # Хз как фиксить, пока просто берем 1-ый элемент из получившегося массива,
        # то есть сам массив [data, data, data ... ]
        if len(noized_params) != 3:
            Z = noized_params
        else:
            Z = np.dot(H, noized_params)[0]

        # Коэффициент, характеризующий эфективность фильтра. Чем коэффициент ниже,
        # тем фильтр лучше фильтрует
        Saf = 0

        K = np.array([[0],
                      [0],
                      [0]])

        P_current = np.array([[0, 0, 0],
                              [0, 0, 0],
                              [0, 0, 0]])

        P = np.array([[1e3, 0,    0    ], 
                      [0,   1e-5, 0    ],
                      [0,   0,    1e-12]])


        # Для расчета Saf 
        X_noize = [self.dVe_noize, self.Fn_noize, self.Wdrn_noize]

        filtred_values = [np.array([[0], 
                                    [0],
                                    [0]])]

        # TODO: commnet all solutions
        for i in range(1, self.iter-1):
            xkd = np.dot(self.F, filtred_values[i - 1])
            P_current = np.dot(np.dot(self.F, P), np.transpose(self.F)) + Q
            K_num = np.dot(P_current, np.transpose(self.H))
            K_den = np.dot(np.dot(self.H, P_current), np.transpose(self.H)) + self.r
            K = K_num / K_den

            filtred_values.append(xkd + np.dot(K, (Z[i-1] - np.dot(self.H, xkd))))

            emkh = self.E - np.dot(K, self.H)

            P = np.dot(emkh, P_current)

        c = 0
        for v in filtred_values:
            self.X_simple_kalman[c] = v[self.idx][0]
            Saf += math.pow((X_noize[0][c] - v[0][0]), 2)
            c += 1
        
        return Saf

    # Адаптивный фильтр первого рода
    def first_order_adaptive_filter(self, noized_params):
        print('First-order adaptive filter ...')

        P = np.array([[1e3, 0,    0    ], 
                      [0,   1e-5, 0    ],
                      [0,   0,    1e-12]])

        P_current = np.array([[0, 0, 0],
                              [0, 0, 0],
                              [0, 0, 0]])

        K = np.array([[0],
                      [0],
                      [0]])

        P_current = np.array([[0, 0, 0],
                              [0, 0, 0],
                              [0, 0, 0]])

        Q = 3e-12        
        Q_mx = np.array([[0, 0, 0],
                         [0, 0, 0],
                         [0, 0, Q]])

        Tar = np.array([[1, 0, 0], 
                        [0, 1, 0],
                        [0, 0, 1]])

        C  = np.zeros(self.iter-1) 
        Nu = np.zeros(self.iter-1) 
        R  = np.zeros(self.iter-1) 

        Saf  = 0
        X_noize = [self.dVe_noize, self.Fn_noize, self.Wdrn_noize]
        filtred_values = [np.array([[0], 
                                    [0],
                                    [0]])]

        H = np.array([[1,0,0]])
        # H = self.H

        if len(noized_params) != 3:
            Z = noized_params
        else:
            Z = np.dot(H, noized_params)[0]

        for i in range(1, self.iter-1):
            xkd = np.dot(self.F, filtred_values[i - 1])
            P_current = np.dot(np.dot(self.F, P), np.transpose(self.F)) + Q_mx
            
            Nu[i] = Z[i] - np.dot(H, xkd)
            i_med = (i - 1) / i
            C[i]  = np.dot(i_med, C[i - 1]) + np.dot(np.dot((1 / i),  Nu), np.transpose(Nu))
            CmPc = C[i] - np.dot(np.dot(H, P_current), np.transpose(H))
            if  CmPc > 0:
                R[i] = CmPc
            else:
                R[i] = 0

            K_num = np.dot(P_current, np.transpose(H))
            K_den = np.dot(np.dot(H, P_current), np.transpose(H)) + R[i]
            K = np.divide(K_num, K_den)

            filtred_values.append(xkd + np.dot(K, Nu[i]))

            emkh = self.E - np.dot(K, H)

            P = np.dot(emkh, P_current)

            Tar = np.dot(Tar, np.dot(emkh, self.F))
            self.Tar_diagram[0][i] = Tar[0][0]
            self.Tar_diagram[1][i] = Tar[1][1]
            self.Tar_diagram[2][i] = Tar[2][2]

        c = 0
        for v in filtred_values:
            self.X_first_order[c] = v[self.idx][0]
            Saf += math.pow((X_noize[self.idx][c] - self.X_first_order[c]), 2)
            c += 1
        
        return Saf

    # Адаптивный фильтр второго рода
    def second_order_adaptive_filter(self, noized_params):
        print('Second-order adaptive filter ...')
        
        P = np.array([[1e3, 0,    0    ], 
                      [0,   1e-5, 0    ],
                      [0,   0,    1e-12]])

        P_current = np.array([[0, 0, 0],
                              [0, 0, 0],
                              [0, 0, 0]])

        K = np.array([[0, 0, 1]])
        Q = self.Q

        C  = np.zeros(self.iter-1) 
        Nu = np.zeros(self.iter-1) 
        R  = 5e-5

        Saf  = 0
        X_noize = [self.dVe_noize, self.Fn_noize, self.Wdrn_noize]
        filtred_values = [np.array([[0], 
                                    [0],
                                    [0]])]
        H = np.array([[1,0,0]])
        # H = self.H

        if len(noized_params) != 3:
            Z = noized_params
        else:
            Z = np.dot(H, noized_params)[0]

        for i in range(1, self.iter-1):
            xkd = np.dot(self.F, filtred_values[i - 1])
            Nu[i] = Z[i] - np.dot(H, xkd)[0][0]
            i_med = (i - 1) / i
            C[i]  = np.dot(i_med, C[i - 1]) + np.dot(np.dot((1 / i),  Nu), np.transpose(Nu))
            
            Q = np.dot(np.dot(K, C[i]), np.transpose(K))
            P_current = np.dot(np.dot(self.F, P), np.transpose(self.F)) + Q
            K_num = np.dot(P_current, np.transpose(H))
            K_den = np.dot(np.dot(H, P_current), np.transpose(H)) + R
            K = np.divide(K_num, K_den)

            filtred_values.append(xkd + np.dot(K, Nu[i]))

            emkh = self.E - np.dot(K, H)
            P = np.dot(emkh, P_current)

        c = 0
        for v in filtred_values:
            self.X_second_order[c] = v[self.idx][0]
            Saf += math.pow((X_noize[self.idx][c] - self.X_second_order[c]), 2)
            c += 1
        
        return Saf

    # Комбайн
    def direct_reverse_filter(self, noized_params):
        print('Direct-reverse combain filter')
        

        # local peaks
        # k = 450
        # median = np.median(noized_params)
        # for i in range(2, 10):
        #     for j in range(i * k -3, i * k + 3):
        #         noized_params[j] = median  * 0.01 * random.random() * random.choice([-1, 1])

        self.simple_kalman(noized_params)
        reverse = [np.array([[0], 
                             [0],
                             [0]])] * (self.iter-1)

        P = np.array([[1e3, 0,    0    ], 
                      [0,   1e-5, 0    ],
                      [0,   0,    1e-12]])

        P_current = np.array([[0, 0, 0],
                              [0, 0, 0],
                              [0, 0, 0]])

        K = np.array([1, 0, 0])


        if len(noized_params) != 3:
            Z = noized_params
        else:
            Z = np.dot(self.H, noized_params)[0]

        for i in range(self.iter -3, -1, -1):
            xkd = np.dot(self.F, reverse[i + 1])
            P_current = np.dot(np.dot(self.F, P), np.transpose(self.F)) + self.Q
            K_num = np.dot(P_current, np.transpose(self.H))
            K_den = np.dot(np.dot(self.H, P_current), np.transpose(self.H)) + self.r
            K = K_num / K_den

            reverse[i] = (xkd + np.dot(K, (Z[i+1] - np.dot(self.H, xkd))))

            emkh = self.E - np.dot(K, self.H)

            P = np.dot(emkh, P_current)

        c = 0
        for v in reverse:
            self.X_direct[c] = self.X_simple_kalman[c]
            self.X_reverse[c] = v[2][0]
            c += 1

        # index = 66
        for i in range(self.iter-1):
            self.X_combine[i] = ((self.X_direct[i] + self.X_reverse[i])/2)
        