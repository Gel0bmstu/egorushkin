import matplotlib.pyplot as plt
import numpy as np

class Settings:
    # ============================================================================================
    # Control flags
    generate_new_params = False  # Флаг перезаписи текущих параметров модели 
                                 # в файлах новыми, сгенерированными значениям

    plot_clear_graphs_flag     = False
    plot_noize_graphs_flag     = False
    plot_mesuarment_noize_flag = False

    plot_shodimost_flag                   = 0

    simple_kalman_flag                    = 0
    plot_first_order_filter_flag          = 0
    plot_second_order_filter_flag         = 0
    plot_all_adaprive_flag                = 0
    plot_combine_filter_flag              = 0
    plot_q_picker_flag                    = 1

    # ============================================================================================
    # Path to files
    clear_params          = './data/clear_params.txt'
    noize_params          = './data/noize_params.txt'
    meausure_noize_params = './data/measurment_noize_params.txt'
    egorushkin_params     = './data/egorushkin_params_test.txt'

    # ============================================================================================
    # Help functions