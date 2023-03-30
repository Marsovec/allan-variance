import numpy as np
import matplotlib.pyplot as plt
import sys, os
from scipy.optimize import curve_fit
import argparse
from datetime import datetime as dt
from matplotlib import rcParams
rcParams["lines.markersize"] = 4

parser = argparse.ArgumentParser(prog="Allan Variance plotter",
                                 description="Plots Overlapping and Non-Overlapping Allan variance.")
parser.add_argument("filename", nargs="+", type=str, help="Path(s) to the data file(s).")
#parser.add_argument("-a", "--avars", default="nonoverlap", choices=["overlap", "nonoverlap", "both"], type=str, help="Choose which variance to calculate (default nonoverlap).")
parser.add_argument("-a", "--append", action="store_true", help="Represent all datapoints on one plot (separate by default).")
parser.add_argument("-l", "--legend", action="store_false", help="Do not draw legend (draw by default).")
args = parser.parse_args()

for filename in args.filename:
    if not os.path.isfile(filename):
        print(f"'{filename}' doesn't exist.")
        exit(1)

def f(x, k, n):
    return k*x+n
def real_f(x, k, n):
    return x**k+10**n-1

logged_xmin = 0
logged_xmax = 0
logged_ymin = 0
logged_ymax = 0
for filename in args.filename:
    data = np.loadtxt(filename, skiprows=2, usecols=(0,1), unpack=True)
    #plt.scatter(*data, label="data")
    #logged = np.log(data)
    logged = np.ma.log(data).filled(0)
    if args.append:
        plt.scatter(*logged, label=filename)
        logged_xmin = np.min([np.min(logged[0]), logged_xmin])
        logged_xmax = np.max([np.max(logged[0]), logged_xmax])
        logged_ymin = np.min([np.min(logged[1]), logged_ymin])
        logged_ymax = np.max([np.max(logged[1]), logged_ymax])

    if not args.append:
        plt.scatter(*logged, label="log-log data")
        popt, pcov = curve_fit(f, *logged)
        k_err = np.sqrt(np.diag(pcov))[0]
        n_err = np.sqrt(np.diag(pcov))[1]
        k = popt[0]
        n = popt[1]
        x = logged[0]
        y = f(x, *popt)
        plt.plot(x, y, label=f"k={round(k, 5)}"+r"$\pm$"+f"{round(k_err, 5)}"+"\n"+f"offset={round(n, 5)}"+r"$\pm$"+f"{round(n_err, 5)}", color="red")
        
        xticks = np.arange(np.min(logged[0]), np.max(logged[0]), np.e)
        xticks_labels = [str(int(i))+r"$e$" for i in xticks/np.e]
        plt.xticks(xticks, xticks_labels)
        
        yticks = np.arange(np.min(logged[1]), np.max(logged[1]), np.e)
        yticks_labels = [str(int(i))+r"$e$" for i in yticks/np.e]
        plt.yticks(yticks, yticks_labels)
        
        plt.xlabel("Group size ln(T)")
        plt.ylabel(r"Allan variance ln($\sigma^2$)")
        if filename[-15:] == "_nonoverlap.dat":
            plt.title("Non-overlapping Allan variance for\n"+filename[:-15])
        elif filename[-12:] == "_overlap.dat":
            plt.title("Overlapping Allan variance for\n"+filename[:-15])
        else:
            plt.title(filename)
        if args.legend:
            plt.legend()
        plt.savefig(filename+"_plot.png", dpi=300)
        print(f"Saved {filename+'_plot.png'}")
        plt.clf()
    elif args.append:
        with open(str(dt.now().strftime('%d.%m.%Y-%H:%M:%S'))+"_plot.txt", "a") as f:
            f.write(filename+"\n")
            
if args.append:
    xticks = np.arange(logged_xmin, logged_xmax, np.e)      
    xticks_labels = [str(int(i))+r"$e$" for i in xticks/np.e]           
    plt.xticks(xticks, xticks_labels)                                   
                                                                        
    yticks = np.arange(logged_ymin, logged_ymax, np.e)      
    yticks_labels = [str(int(i))+r"$e$" for i in yticks/np.e]           
    plt.yticks(yticks, yticks_labels)                                   
                                                                        
    plt.xlabel("Group size ln(T)")                                      
    plt.ylabel(r"Allan variance ln($\sigma^2$)")                        
    plt.title("Allan variance, appended together\n"+str(dt.now().strftime('%d.%m.%Y %H:%M:%S')))
    if args.legend:
        plt.legend()                                                     
    plt.savefig(str(dt.now().strftime('%d.%m.%Y-%H:%M:%S'))+"_plot.png", dpi=300)                          
    print(f"Saved {dt.now().strftime('%d.%m.%Y-%H:%M:%S')}_plot.png")                              
    plt.clf()                                                           
