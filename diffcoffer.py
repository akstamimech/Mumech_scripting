import Dfit.Dfit as Dfit
import os


def dfiter(files):
    res = Dfit.Dcov(fz = files, m = 20, tmin = 50, tmax = 2500, dt = 0.2, nitmax = 500, imgfmt = "pdf")
    res.run_Dfit()
    res.analysis(tc = 0.4)
    res.finite_size_correction(L = 30.6422, eta = 0.001, tc = 0.4)

def index_dat_files():
    dat_files = ["coords_kap0-COM.dat"]

    if not dat_files:
        print("No .dat files found in the directory.")
        return

    print("Found the following .dat files:")
    for dat_file in dat_files:
        print(dat_file)

    dfiter(dat_files)

# Replace 'your_directory_path' with the actual path to your directory
#directory_path = 'your_directory_path'
index_dat_files()
