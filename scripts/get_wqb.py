import numpy as np
import os,sys
import  matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt

def get_wqb_simple(file_duck_dat):
    f = open(file_duck_dat,'r')
    data = []
    for line in f:
        a = line.split()
        data.append([float(a[1]), float(a[3]), float(a[5]), float(a[8])])
    f.close()
    data = np.array(data[1:])
    Work = data[:,3]
    Wqb_max = max(Work[400:])
    Wqb_min = min(Work[:400])
    Wqb_value = Wqb_max - Wqb_min
    return(Wqb_value, data, Wqb_min)


def get_Wqb_value_all(input_dir):
    file_list = []
    for fil in os.listdir(input_dir):
        if fil[-3:] == 'dat':
            file_list.append(fil)

    Wqb_values = []
    plt.figure(figsize = (7,7))
    for fil in file_list:
        Wqb_data = get_wqb_simple(fil)
        Wqb_values.append(Wqb_data[0])
        plt.plot(1*Wqb_data[1][:,0], Wqb_data[1][:,3]-Wqb_data[2])

    plt.xlabel('HB Distance (A)')
    plt.ylabel('Work (kcal/mol)')
    plt.savefig('wqb_plot.png')
    Wqb = min(Wqb_values)
    return(Wqb)

def main():
    if len(sys.argv) > 1:
        input_dir = sys.argv[1]
    else:
        input_dir = os.getcwd()
    print(get_Wqb_value_all(input_dir))

if __name__ == '__main__':
    main()

