import numpy as np

class Vars:
    # ============================================================================================
    # Constants
    I     = np.eye(3) # Identity matrix
    R     = 6371000   # Earth radius               [m]
    g     = 9.81      # Gravitational acceleration [g/ms^2]
    U     = 7.27e-5   # Earth angular rotation     [rad/s]     
    t     = 5e3       # Simulation time            [sec]
    mu    = 0         # Expected value
    sigma = 1         # SKO
    b     = 5e-8      # Сorrelation parameter (beta)
    A     = 1e-4      # TODO: comment
    rate  = 1         # Model rate                 [Hz]
    Q     = 3.3e-13   # TODO: comment
    Q_a1  = 1e-14
    Q_a2  = 1e-14
    r     = 1e-2      # TODO: comment

    iter = int(t * 1/rate) # Model iterrations count [n]
    figure_index = 1

    E  = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
    # Ev  = np.array([[1],
    #                 [1],
    #                 [1]])
    Ev  = np.array([1,1,1])
    # ============================================================================================
    # Initial conditions
    dVe_0  = 0
    Fn_0   = 0
    Wdrn_0 = 0

    # ============================================================================================
    # Vars

    # "Вектор состояний" позволяет выбрать фильтруемую величину, тоесть
    # если на вох фильтра подаются все значения зашумленные (Ve, Fn, Wdr), 
    # а значение вектор H = (1,0,0), то отфильтруется только Ve
    H = np.array([[1, 0, 0]])

    # Noize data vars 
    W = np.random.normal(mu, sigma, iter) # White noize

    # TODO: comment
    F = np.array([[1,  -g*rate, 0],
                  [1/R/rate, 1, 1],
                  [0,   0, (1 - b) * 1/rate]])

    # Вектор значений, фактически - вектор векторов:
    #       | dVe |
    # X = [ | Fn  | ]
    #       | Wdr |
    X = [np.array([[0], 
                   [0],
                   [0]])]

    # Filters solution arrays
    X_simple_kalman   = [0] * (iter - 1)
    X_first_order     = [0] * (iter - 1)
    X_second_order    = [0] * (iter - 1)
    X_direct          = [0] * (iter - 1)
    X_reverse         = [0] * (iter - 1)
    X_combine         = [0] * (iter - 1)
    X_Q_picker        = []

    Q_peaker_legend = []
    Q_line_width    = []
    Q_line_style    = []

    # TODO: commnet
    G = np.array([[0], 
                  [0],
                  [A * pow(2 * b, 1/2)]])
                #   [1]])

    GQG=np.array([[Q, 0, 0],
                  [0, Q, 0],
                  [0, 0, Q]])

    F = np.array([[1,                      -g * 1 / rate, 0                  ],
                  [1 / R / rate,           1,             1 / rate           ],
                  [0,                      0,             (1 - b) * 1 / rate]])

    # Фильтры фильтруют только одно из 3-х значений модели (например,
    # только dVe). Для упрощения расчета различных параметров, например
    # коэффициента Saf, введем переменную idx, которая будет указывать на
    # нужный нам столбик в матрицах решений фильтров (X_kd и т.д.) 
    # idx = np.dot(np.array([0, 1, 2]), np.transpose(H))[0]
    idx = 2

    # Коэффициент, характеризующий эфективность фильтра. Чем коэффициент ниже,
    # тем фильтр лучше фильтрует
    Saf = 0

    Tar_diagram = ([0] * (iter - 1),
                   [0] * (iter - 1),
                   [0] * (iter - 1)) 

    # Clear params
    dVe_clear  = [] # Clear linear velocity
    Fn_clear   = [] # Clear angle
    Wdrn_clear = [] # Clear angular velocity

    # Noized params
    dVe_noize  = []
    Fn_noize   = []
    Wdrn_noize = []

    # Measurment noize params
    dVe_mesuarment_noize  = []
    Fn_mesuarment_noize   = []
    Wdrn_mesuarment_noize = []