
import csv
from io import StringIO
import numpy as np

file = 'C:/Users/Mattia Mencagli/Documents/DATA/output_lizard_POST.csv'

print(f"Lizard .csv Header:\n0: NLOC; \n1: CCN; \n2: Token; \n3: Param; \n4: Number of lines; \
\n5: Detailed func name; \n6: File name; \n7: Func name; \n8: Func call; \n9: First line; \n10: Last line ")

toint = [0,1,2,3,4]
liz = []
buffer = ""

with open(file, 'r') as f:
    for line in f:
        line = line.rstrip('\n')  # Remove trailing newline for function which are "truncated"

        if line.endswith('\\'):  # If line ends with `\`, join with next line
            buffer += line[:-1]  # Remove `\` and add to buffer
        else:
            if buffer:  # If buffer is not empty, prepend it to the current line
                line = buffer + line
                buffer = ""  # Reset buffer

            # Parse the line with csv.reader
            reader = csv.reader(StringIO(line))
            try:
                tmp = next(reader)
                # Convert specified columns to int
                for i in toint:
                    tmp[i] = int(tmp[i])
                liz.append(tmp)
            except StopIteration:
                # Skip empty lines
                continue


total = len(liz)

CCNlim = 1
highCCN = []
summ = 0 
for l in liz:
    if l[1]>CCNlim:
        summ += 1
        highCCN.append([l[0],l[1],l[5]])
print(f"tot={total}; CCN>{CCNlim}? {summ}")

NLOCs = np.array([sublist[0] for sublist in liz])
NCCs  = np.array([sublist[1] for sublist in liz])

avg_NLOC = np.average(NLOCs)
avg_CCN = np.average(NCCs)
print(f"avg_NLOC:{avg_NLOC:.1f}; avg_CCN:{avg_CCN:.1f}")

for el in highCCN:
    if 'GeoLocalization::passive' in el[2]:
        print(el)
    elif 'GeoLocalization::filter_ned_passive' in el[2]:
        print(el)
    elif 'GeoLocalization::check_finale_passive' in el[2]:
        print(el)
    elif 'GeoLocalization::IIeq_passive_intersections' in el[2]:
        print(el)
