import os
import h5py
import numpy as np
import pandas as pd

# InteractiveShell.ast_node_interactivity = "all"
# plt.style.use("default")

# f = h5py.File('./../phd_data.hdf5', 'a')
# f = h5py.File('./../data/phd_data_temp.hdf5', 'a')


number = 3
folder = "./../experimentos/20230720/test_"
fbg_n = "fbg_6"
date = "20230720"
# f.close()
data = []
files = (
    os.popen("ls " + folder + str(number) + "/" + fbg_n + "*.csv -r")
    .read()
    .split("\n")[:-1]
)
for i in files:
    print(i)
    data.append(pd.read_csv(i))

len(data)


f = h5py.File("phd_data.hdf5", "a")
# ff = f.require_group("accel_1/allan_variance")
# ff.attrs["calibrate_angular_coefficient"] = coef[0]
# ff.attrs["calibrate_linear_coefficient"] = coef[1]

# ff = f.require_group("accel_1/allan_variance/y_up_long")
# ff["y_calibrated"] = calibrate_data(ff["y"][:] / ff["tap"][:])


# fff = ff.require_group("test_of_pm_components")


# fff["sld_isolator_circulator/v_trans"] = pd.read_csv("./../accel_v1/20240111/circulator_analisys",sep='\t').iloc[:,1].to_numpy()
# fff["sld_isolator_circulator"].attrs["info"]= "Saída do transimpedância com o SLD ligado ao isolador e ao acoplador 99/1 e ao circulador."


# fff["sld_isolator_coupler/time"] = pd.read_csv("./../accel_v1/20240111/sld_isolator_acoplador.txt").iloc[:,0].to_numpy()
# fff["sld_isolator_coupler/v_trans"] = pd.read_csv("./../accel_v1/20240111/sld_isolator_acoplador.txt").iloc[:,1].to_numpy()
# fff["sld_isolator_coupler"].attrs["info"]= "Saída do transimpedância com o SLD ligado ao isolador e ao acoplador 99/1."


# fff["sld_isolator/time"] = pd.read_csv("./../accel_v1/20240111/sld_isolator.txt").iloc[:,0].to_numpy()
# fff["sld_isolator/v_trans"] = pd.read_csv("./../accel_v1/20240111/sld_isolator.txt").iloc[:,1].to_numpy()
# fff["sld_isolator"].attrs["info"]= "Saída do transimpedância com o SLD ligado ao isolador."


# fff["sld/time"] = pd.read_csv("./../accel_v1/20240111/sld.txt").iloc[:,0].to_numpy()
# fff["sld/v_trans"] = pd.read_csv("./../accel_v1/20240111/sld.txt").iloc[:,1].to_numpy()
# fff["sld"].attrs["info"]= "Saída do transimpedância com o SLD."

# f = h5py.File("./../phd_data.hdf5", "a")
# fff = f.require_group("accel_1/allan_variance")

# fff.create_dataset(
#     "y_up_long/time",
#     data=pd.read_csv("./../accel_v1/20240117/allan_variance_longo")
#     .iloc[:, 0]
#     .to_numpy(),
#     compression="gzip",
# )

# fff.create_dataset(
#     "y_up_long/y",
#     data=pd.read_csv("./../accel_v1/20240117/allan_variance_longo")
#     .iloc[:, 1]
#     .to_numpy(),
#     compression="gzip",
# )

# fff.create_dataset(
#     "y_up_long/tap",
#     data=pd.read_csv("./../accel_v1/20240117/allan_variance_longo")
#     .iloc[:, 2]
#     .to_numpy(),
#     compression="gzip",
# )

# fff["sld"].attrs["info"]= "Saída do transimpedância com o SLD ligado ao isolador."
# fff["sld/time"] = pd.read_csv("./../accel_v1/20240111/sld.txt").iloc[:,0].to_numpy()
# fff["sld/v_trans"] = pd.read_csv("./../accel_v1/20240111/sld.txt").iloc[:,1].to_numpy()


def read_all_data():
    _N = 0
    for i, count in zip(data, range(1, 1 + len(files))):
        _tt = f.require_group(
            "optical_mechanical_identification/test_fibers/"
            + fbg_n
            + "/test_00"
            + str(number)
            + "/"
            + fbg_n
            + "_"
            + "{:03d}".format(_N + count)
            + "/"
        )
        for k in [
            "actual_resolution",
            "ilx_current",
            "traction_N",
            "micrometer_position_um",
            "resolution_nm",
            "resolution_vbw",
            "initial_length_m",
            "room_temperature_C",
        ]:
            _tt.attrs[k] = data[count - 1][k][0]
        _tt.attrs["date"] = date
        for k in ["wavelength", "power_dbm"]:
            _tt[k] = data[count - 1][k]


read_all_data()
data = []
files = (
    os.popen("ls " + folder + str(number) + "/fonte*.csv -r").read().split("\n")[:-1]
)
for i in files:
    print(i)
    data.append(pd.read_csv(i))
len(data)


def read_all_data():
    _N = 0
    for i, count in zip(data, range(1, 1 + len(files))):
        _tt = f.require_group(
            "optical_mechanical_identification/test_fibers/"
            + fbg_n
            + "/test_00"
            + str(number)
            + "/fonte"
            + "{:03d}".format(_N + count)
            + "/"
        )
        for k in [
            "actual_resolution",
            "ilx_current",
            "resolution_nm",
            "resolution_vbw",
            "room_temperature_C",
        ]:
            _tt.attrs[k] = data[count - 1][k][0]
        _tt.attrs["date"] = date
        for k in ["wavelength", "power_dbm"]:
            _tt[k] = data[count - 1][k]


read_all_data()
f.close()
