import matplotlib.pyplot as plt
import numpy as np

from vars import Vars
from settings import Settings

class Plotter(Settings, Vars):
    def __init__(self):
        pass

    def plot(self, data, xlabel, ylabel, title, legend=None, linewidth=None, linestyle=None):
        if legend is not None and len(data) != len(legend):
            print('Unable to plot graph, data length and legend width must be same')
            return
        
        if linewidth is not None and len(data) != len(legend):
            print('Unable to plot graph, data length and linewidth width must be same')
            return

        width = 20
        height = 12
        
        plt.figure(figsize=(width, height), num=self.figure_index)    
        
        c = 0
        for each in data:
            if linewidth is None:
                plt.plot(each, linewidth = 0.5, linestyle=linestyle[c])
                c += 1
            else:
                plt.plot(each, linewidth=linewidth[c], linestyle=linestyle[c])
                c += 1

        plt.suptitle(title, fontsize=20)
        plt.xlabel(xlabel, fontsize = 30)
        plt.ylabel(ylabel, fontsize = 30)
        plt.xticks(fontsize = 20)
        plt.yticks(fontsize = 20)
        if legend != None:
            plt.legend(legend, prop={'size': 20})
        plt.grid()

        self.figure_index += 1

    def run(self):
        X = [self.dVe_mesuarment_noize, self.Fn_mesuarment_noize, self.Wdrn_mesuarment_noize]

        if self.plot_clear_graphs_flag:
            self.plot([self.dVe_clear], "Time, s", "Ve, m/s", 'Clear')
            self.plot([self.Fn_clear], "Time, S", "Fn, rad", 'Clear')
            self.plot([self.Wdrn_clear], "Time, S", "Wdrn, rad/s", 'Clear')

        if self.plot_noize_graphs_flag:
            self.plot([self.dVe_noize], "Time, s", "Ve, m/s", 'Noize')
            self.plot([self.Fn_noize], "Time, S", "Fn, rad", 'Noize')
            self.plot([self.Wdrn_noize], "Time, S", "Wdrn, rad/s", 'Noize')
        
        if self.plot_mesuarment_noize_flag:
            self.plot([self.dVe_mesuarment_noize], "Time, s", "Ve, m/s", 'dVe Measure noize')
            self.plot([self.Fn_mesuarment_noize], "Time, s", "Fn, rad", 'Fn Measure noize')
            self.plot([self.Wdrn_mesuarment_noize], "Time, s", "Wdrn, rad/s", 'Wdrn Measure noize')
        
        if self.plot_shodimost_flag:
            self.plot(self.Tar_diagram, "Time, s", "Tar factors", 'Сходимость')

        if self.plot_q_picker_flag:
            self.plot(self.X_Q_picker,\
                "Time, s", "wrdn, rad", 'All Q picks',\
                legend=self.Q_peaker_legend,
                linewidth=self.Q_line_width,
                linestyle=self.Q_line_style
            )

        if self.simple_kalman_flag:
            self.plot(self.X_simple_kalman,\
                "Time, s", "wrdn, rad", 'Simple kalman',\
                legend=['Noized', 'Filtered'],\
                linewidth=[0.5, 1.5])

        if self.plot_first_order_filter_flag:
            self.plot([X[self.idx], self.X_first_order],\
                 "Time, s", "wrdn, rad", 'First order', ['Noized', 'Filtered'], [0.5, 1.5])

        if self.plot_second_order_filter_flag:
            self.plot([X[self.idx], self.X_second_order],\
                 "Time, s", "wrdn, rad", 'Second order', ['Noized', 'Filtered'], [0.5, 1.5])

        if self.plot_all_adaprive_flag:
            self.plot([X[self.idx], self.X_first_order, self.X_second_order],\
                 "Time, s", "wrdn, rad", 'All adaptive', \
                 ['Noized', 'First order', 'Second order'], \
                 [0.5, 1.5, 1.5])

        if self.plot_combine_filter_flag:
            self.plot([X[self.idx], self.X_direct, self.X_reverse, self.X_combine],\
                "Time, s", "wrdn, rad", 'Combine order', \
                legend=['Noized', 'Direct', 'Reverse', 'Combine'], \
                linewidth=[0.5, 1.5, 1.5, 1.5])

        plt.show()
