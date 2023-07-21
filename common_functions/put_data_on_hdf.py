
import os
import h5py
import numpy as np
import pandas as pd

# InteractiveShell.ast_node_interactivity = "all"
# plt.style.use("default")

f = h5py.File('./../data/phd_data.hdf5', 'a')
# f = h5py.File('./../data/phd_data_temp.hdf5', 'a')


number= 3
folder = './../experimentos/20230720/test_'
fbg_n = 'fbg_6'
date = '20230720'
 # f.close()
data = []
files = os.popen(
    'ls '+folder+str(number)+'/'+fbg_n+'*.csv -r').read().split('\n')[:-1]
for i in files:
    print(i)
    data.append(pd.read_csv(i))
    
len(data)

# for i in range(len(data)):
#     print(data[i]['initial_length_m'][0])
#     data[i]['initial_length_m'] = 39.7e-3
#     data[i].to_csv(files[i])

def read_all_data():
    _N = 0
    for i, count in zip(data, range(1, 1 + len(files))):
        _tt = f.require_group("optical_mechanical_identification/test_fibers/" + fbg_n + "/test_00" + str(number) + "/" + fbg_n + "_" + "{:03d}".format(_N + count) + "/")
        for k in ['actual_resolution', 'ilx_current', 'traction_N', 'micrometer_position_um',
                  'resolution_nm', 'resolution_vbw', 'initial_length_m', 'room_temperature_C']:
            _tt.attrs[k] = data[count - 1][k][0]
        _tt.attrs['date'] = date
        for k in ['wavelength', 'power_dbm']:
            _tt[k] = data[count - 1][k]


read_all_data()
# f.close()


# READ SOURCE


data = []
files = os.popen(
    'ls '+folder+str(number)+'/fonte*.csv -r').read().split('\n')[:-1]
for i in files:
    print(i)
    data.append(pd.read_csv(i))
len(data)


def read_all_data():
    _N = 0
    for i, count in zip(data, range(1, 1 + len(files))):
        _tt = f.require_group(
            "optical_mechanical_identification/test_fibers/" + fbg_n + "/test_00" + str(number) + "/fonte" +
            "{:03d}".format(
                _N +
                count) +
            "/")
        for k in ['actual_resolution', 'ilx_current',
                  'resolution_nm', 'resolution_vbw', 'room_temperature_C']:
            _tt.attrs[k] = data[count - 1][k][0]
        _tt.attrs['date'] = date
        for k in ['wavelength', 'power_dbm']:
            _tt[k] = data[count - 1][k]


read_all_data()
f.close()


# f = h5py.File('./../data/phd_data.hdf5', 'a')
# ff = f['optical_mechanical_identification/test_fibers/fbg_2']
# fff = ff.require_group('test_001')
# # fff2 = fff.require_group('fbg_2')
# # fff3 = fff.require_group('fbg_3')
# # my_keys = list(ff.keys())

# for key in ff.keys():
#     print(key)
#     if key.startswith("fbg_2") or key.startswith("fonte_2"):
#         ff.move(key, 'test_001/' + key)
# #     elif key.startswith("fonte_3"):
#         ff.move(key, 'test_fibers/fbg_3/'+key)


