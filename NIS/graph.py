"""
The equivalent of 95% confidence level is a 0.05 percent error rate on the NIS. Being a Chi squared distribution
and having:
 - 2 degrees of freedom for LIDAR
 - 3 degrees of freedom for RADAR

Values taken from the table provided at the following link:

    https://en.wikipedia.org/wiki/Chi-squared_distribution#Table_of_χ2_values_vs_p-values
"""

import matplotlib.pyplot as plt

LIDAR_CONFIDENCE_LEVEL_95 = 5.99
RADAR_CONFIDENCE_LEVEL_95 = 7.82


def load_data(path):
    with open(path, 'r') as f:
        return [float(x) for x in f.read().split()]


def show_graph(data, confidence_level, out_name, nis_label):
    fig, ax = plt.subplots()

    print(out_name)
    print(data)
    ax.plot(data, ms=20, lw=2, alpha=0.7, label=nis_label)
    ax.axhline(y=confidence_level, color='r', linestyle='-', label="95%% (%s)" % confidence_level)
    ax.grid()
    ax.legend()

    plt.savefig("plot_%s.png" % out_name)


if __name__ == '__main__':
    lidar_data = load_data("NIS_lidar.txt")
    radar_data = load_data("NIS_radar.txt")

    show_graph(lidar_data, LIDAR_CONFIDENCE_LEVEL_95, "LIDAR", "NIS lidar")
    show_graph(radar_data[1:], RADAR_CONFIDENCE_LEVEL_95, "RADAR", "NIS radar")
    # the answer -> *