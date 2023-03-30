import numpy as np
import argparse, os, time
import pandas as pd
from multiprocessing import Pool
from datetime import datetime as dt
pd.set_option('display.float_format', lambda x: '%.6e' % x)

def cluster_avgs_sum(seq, T):
    #seems to take about 16/T seconds on a 1 MB block
    N = seq.size #number of measurements
    K = int(N/T) #number of clusters in truncated seq
    N = K*T #calculate length of truncated seq
    #calculate per cluster averages
    #t = time.time()
    cluster_avgs = np.sum([seq[i:i+T] for i in range(0, N, T)], axis=1)/T
    result = np.sum(np.square([cluster_avgs[k+1]-cluster_avgs[k] for k in range(0, K-1)]))
    #print(time.time()-t,"s")
    return result
def avar_nonoverlap_calc(cluster_sum, K):
    return (1/(2*(K-1)))*cluster_sum
def avar_nonoverlap_sum(avar1, avar2, K1, K2, remainder=0):
    return 1/(2*(K1+K2-1))*(2*(K1-1)*avar1+2*(K2-1)*avar2+remainder)

parser = argparse.ArgumentParser(prog="Allan Variance calculator",
                                 description="Calculates Overlapping and Non-Overlapping Allan variance of a specified data block.")
parser.add_argument("filename", nargs="+", type=str, help="Path(s) to the data file(s).")
parser.add_argument("-b", "--blocks", default=1, type=float, help="Number of 1MB blocks to read concurrently aka chunk size (default 1).")
parser.add_argument("-a", "--avars", default="nonoverlap", choices=["overlap", "nonoverlap", "both"], type=str, help="Choose which variance to calculate (default nonoverlap).")
parser.add_argument("-T", "--tmin", default=0, type=int, help="Smallest cluster size exponent T=2^(tmin) to calculate with (default tmin=0).")
parser.add_argument("-s", "--sum", action="store_true", help="Whether to sum each chunk with the previous (default false).")
parser.add_argument("-c", "--count", default=0, type=int, help="Count of chunks to consider (default all).")
parser.add_argument("-o", "--offset", default=0, type=float, help="File offset in MB (default 0).")
parser.add_argument("-d", "--outdir", default="", type=str, help="Save avars to .dat files in a specified directory.")
args = parser.parse_args()
if args.tmin < 0:
    args.tmin = 0
if args.offset < 0:
    args.offset = 0
    
for filename in args.filename:
    if not os.path.isfile(filename):
        print(f"'{filename}' doesn't exist.")
        exit(1)
        
chunk_size = int(args.blocks*1e6*8)
print(f"#chunk size: {chunk_size}")
args.offset = int(args.offset*1e6*8)
print(f"#offset: {args.offset}")

if args.outdir:
    if not os.path.exists(args.outdir):
        os.mkdir(args.outdir)
    with open(os.path.join(args.outdir, "_INFO.txt"), "w") as f:
        f.write(f"Calculations started on {dt.now().strftime('%d.%m.%Y %H:%M:%S')}\n")
        f.write(f"chunk size: {chunk_size}\n")
        f.write(f"offset: {args.offset}\n")
        f.write(f"avars: {args.avars}\n")
        f.write("summed: "+str(args.sum))

for filename in args.filename:
    print(f"#filename: {filename}")    
    file_size = os.path.getsize(filename)*8 #in bits
    print(f"#file size: {file_size}")
    if chunk_size > file_size or chunk_size <= 0:
        chunk_size = file_size
    if args.offset+chunk_size > file_size:
        print("Offset (and chunk size) bigger than file size.")
        continue
    Ts = np.array(2**np.arange(args.tmin, np.log(chunk_size/2)/np.log(2)), dtype=int)
    #Ts = np.array(np.arange(args.tmin+1, chunk_size/2), dtype=int)
    Ks = np.array(chunk_size/Ts, dtype=int)
    avar_nonoverlap = {}
    for T in Ts:
        avar_nonoverlap[T] = 0
    avar_overlap = avar_nonoverlap.copy()
    
    #array of chunk offsets
    offsets = np.array([i for i in range(args.offset, file_size-chunk_size+1, chunk_size)])
    if args.count > 0 and file_size - (args.offset+args.count*chunk_size) >= chunk_size:
        offsets = offsets[:args.count]
    #print(f"#chunk positions: {', '.join([str(i) for i in offsets])} (bits)")
    
    last1 = avar_overlap.copy() #EXPERIMENTAL
    first2 = avar_overlap.copy() #EXPERIMENTAL
    #cluster_cum = avar_overlap.copy()
    for offset, i in zip(offsets, range(0, offsets.size)):
        print("#")
        print(f"#POSITION: {offset}")
        chunk = np.fromfile(filename, dtype=np.dtype("B"), count=int(chunk_size/8), offset=int(offset/8))
        chunk = np.unpackbits(chunk)
        if args.avars == "overlap" or args.avars == "both":
            print("#avar_overlap not done yet")
        if args.avars == "nonoverlap" or args.avars == "both":
            for T, K in zip(Ts, Ks):
                if args.sum: #EXPERIMENTAL
                    first2[T] = np.sum(chunk[:T])/T
                    if i == 0:
                        remainder = 0
                    else:
                        remainder = np.square(first2[T]-last1[T])
                    avar_nonoverlap[T] = avar_nonoverlap_sum(avar_nonoverlap[T], avar_nonoverlap_calc(cluster_avgs_sum(chunk, T), K), i*K, K, remainder)
                    last1[T] = np.sum(chunk[-T:])/T
                    #cluster_cum[T] += cluster_avgs_sum(chunk, T)
                    #avar_nonoverlap[T] = avar_nonoverlap_calc(cluster_cum[T], (i+1)*K)
                else:
                    avar_nonoverlap[T] = avar_nonoverlap_calc(cluster_avgs_sum(chunk, T), K)
        #format dataframe
        avar_nonoverlap_df = pd.DataFrame(avar_nonoverlap.values(), index=avar_nonoverlap.keys())
        avar_nonoverlap_df = "\n".join(avar_nonoverlap_df.to_string().split("\n")[1:])
        print(f"#T{(len(str(T)))*' '}avar_nonoverlap(T)")
        print(avar_nonoverlap_df)
        #save to a file, if necessary
        if args.outdir:
            if args.avars == "nonoverlap" or args.avars == "both":
                with open(os.path.join(args.outdir, filename.replace("/", "_")+"_offset="+str(offset)+"_nonoverlap.dat"), "w") as f:
                    f.write(f"#POSITION: {offset}\n")
                    f.write(f"#T{(len(str(T)))*' '}avar_nonoverlap(T)\n")
                    f.write(avar_nonoverlap_df)
    print("---")
