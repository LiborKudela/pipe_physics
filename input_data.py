from scipy.interpolate import interp1d

class Input_data():
    def __init__(self, t, tws, twr, ta, alpha_s=600, alpha_r=600, alpha_a=20):
        self.alpha_s = alpha_s
        self.alpha_r = alpha_r
        self.alpha_a = alpha_a
        self.tws = interp1d(t, tws, fill_value="extrapolate")
        self.twr = interp1d(t, twr, fill_value="extrapolate")
        self.ta = interp1d(t, ta, fill_value="extrapolate")
        self.t_stop = max(t)

    def T_ws(self, t):
        return float(self.tws(t))

    def T_wr(self, t):
        return float(self.twr(t))

    def T_a(self, t):
        return float(self.ta(t))

default_training_data = Input_data(
        [0, 5000, 1200000, 1205000, 2400000, 2405000, 7589000], # time
        [5, 10, 10, 10, 10, 10, 10],                            # t_ws
        [5, 5, 5, -10, -10, -10, -10],                          # t_wr
        [5, 5, 5, 5, 5, 0, 0])                                  # ta  

default_test_data = Input_data(
        [0, 200000, 600000, 1200000, 2400000, 4800000],  # time
        [5, 10, 10, 15, 30, 30],                         # t_ws
        [5, -10, -10, -15, -30, -30],                    # t_wr
        [5, -5, -5, 0, 0, 0])                            # ta  
