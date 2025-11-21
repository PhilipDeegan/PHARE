import os
import pandas
import subprocess
import numpy as np

import matplotlib.pyplot as plt

from pathlib import Path


def main():
    pandas.read_csv("tota.csv").astype(float).plot(kind="line", title="8 ranks")

    plt.show()


if __name__ == "__main__":
    main()
