from form_model_params import ModelParamsFormer

from filter import FiltersPool
from plotter import Plotter

if __name__ == "__main__":
    model_former = ModelParamsFormer()
    filters_pool = FiltersPool()
    
    noized_values = model_former.get_measurment_noize_params()
    # filters_pool.Q_picker(noized_values)

    # noized_values = model_former.get_egorushkin_params()
    filters_pool.simple_kalman(noized_values)
    filters_pool.first_order_adaptive_filter(noized_values)
    filters_pool.second_order_adaptive_filter(noized_values)
    filters_pool.direct_reverse_filter(noized_values)

    p = Plotter()
    p.run()