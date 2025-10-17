"""
This program is based on a one-pass algorithm for the calculation of an array of autocorrelations r[1], r[2], ... r[K].
The key feature of this algorithm is the circular array 'hold' which stores the (K + 1) most recent data points and the
associated index 'p' which points to the (rotating) head of the array.

Data is read from a text file in the format 1-data-point-per-line (with no blank lines).

NOTE: the constant K (maximum lag) MUST be smaller than the of data points in the text file, n. Moreover, if the
autocorrelations are to be statistically meaningful, K should be MUCH smaller than n.

Name            : acs.c  (AutoCorrelation Statistics)
Author          : Steve Park & Dave Geyer
Language        : ANSI C
Latest Revision : 10/02/1997
Compile with    : gcc -lm acs.c
Execute with    : acs.out < acs.dat
Translated by   : Philip Steele
Language        : Python 3.3
Latest Revision : 26/03/2014
Execute with    : python acs.py < acs.dat
Updated by      : Andrea Mazzucchi
Python Version  : 3.17
Latest Revision : 17/10/2025
Execute with    : python acs.py acs.dat
"""

import sys

import numpy as np

K = 50  # K is the maximum lag
SIZE = (K + 1)

_sum = 0.0  # sums x[i]
hold = []  # K + 1 most recent data points
p = 0  # points to the head of 'hold'
co_sum = [0 for _ in range(0, SIZE)]  # co_sum[j] sums x[i]  x[i+j]

file_path = None
if len(sys.argv) > 1:
    file_path = sys.argv[1]
else:
    print("Error: Please provide a file path as the first argument.")

i = 0  # data point index
with open(file_path, 'r') as file:
    for line in file:
        x = float(line)
        _sum += x
        hold.append(x)
        i += 1
        if i == SIZE:
            break

    for line in file:
        for j in range(0, SIZE):
            co_sum[j] += hold[p] * hold[(p + j) % SIZE]
        x = float(line)  # lines read in as string
        _sum += x
        hold[p] = x
        p = (p + 1) % SIZE
        i += 1

n = i  # the total number of data points

while i < n + SIZE:  # empty the circular array
    for j in range(0, SIZE):
        co_sum[j] += hold[p] * hold[(p + j) % SIZE]
    hold[p] = 0.0
    p = (p + 1) % SIZE
    i += 1

mean = _sum / n
for j in range(0, K + 1):
    co_sum[j] = (co_sum[j] / (n - j)) - (mean * mean)

print("for {0} data points".format(n))
print("the mean is ... \t{0:8.2f}".format(mean))
print("the std_dev is ... \t{0:8.2f}\n".format(np.sqrt(co_sum[0])))
print("  j (lag)   r[j] (autocorrelation)\n")
for j in range(1, SIZE):
    print("{0:3d}  {1:11.3f}".format(j, co_sum[j] / co_sum[0]))

# Output with acs.dat:
# for 500 data points
# the mean is ...       5.68
# the std_dev is ...    4.72

#   j (lag)   r[j] (autocorrelation)
#   1        0.966
#   2        0.931
#   3        0.898
#   4        0.862
#   5        0.831
#   6        0.799
#   7        0.764
#   8        0.730
#   9        0.700
#  10        0.668
#  11        0.641
#  12        0.612
#  13        0.583
#  14        0.554
#  15        0.524
#  16        0.494
#  17        0.465
#  18        0.437
#  19        0.410
#  20        0.382
#  21        0.358
#  22        0.334
#  23        0.308
#  24        0.277
#  25        0.247
#  26        0.221
#  27        0.193
#  28        0.168
#  29        0.148
#  30        0.131
#  31        0.119
#  32        0.105
#  33        0.089
#  34        0.076
#  35        0.069
#  36        0.064
#  37        0.062
#  38        0.062
#  39        0.061
#  40        0.056
#  41        0.054
#  42        0.060
#  43        0.062
#  44        0.067
#  45        0.068
#  46        0.059
#  47        0.055
#  48        0.053
#  49        0.052
#  50        0.050
