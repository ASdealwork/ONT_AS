import numpy as np
import argparse

# 读取文件
#data = np.loadtxt('chip1_PAM11408_PAM09650_length_num.txt', dtype=int)

parser = argparse.ArgumentParser(description='Process file for binning data.')
parser.add_argument('filename1', type=str, help='Name of the file to process')
parser.add_argument('filename2', type=str, help='Name of the file to process')
args = parser.parse_args()

# 读取文件
data = np.loadtxt(args.filename1, dtype=int)
#out_f2 = np.loadtxt(args.filename2, dtype=int)

bin_width = 8000
# 计算total_bins
max_value_second_column = np.max(data[:, 1])  # 获取第二列的最大值
total_bins = int(np.ceil(max_value_second_column / bin_width))  # 计算总的bin数

# 打印结果
print("Total bins:", total_bins)

# 计算每个bin的范围
bin_ranges = []
for i in range(total_bins):
    start = (bin_width) * i
    end = bin_width * (i + 1)
    bin_ranges.append((start, end))


# 计算每个bin范围内的总Bases数
bin_bases = []
for bin_range in bin_ranges:
    start, end = bin_range
    bin_indices = np.where((data[:, 1] >= start) & (data[:, 1] <= end))[0]
    bin_reads = data[bin_indices, 0]
    bin_lengths = data[bin_indices, 1]
    bin_total_bases = np.sum(bin_reads * bin_lengths)
    bin_bases.append(bin_total_bases / 1000000)

# 计算横轴的中心点（除以1000）
bin_centers = [(bin_range[0] + bin_range[1]) // 2 for bin_range in bin_ranges]

# 打印每个bin的范围和对应的Bases数
for i in range(total_bins):
    print(f"Bin {i+1}: Range {bin_ranges[i]}, Bases {bin_bases[i]}")

# 打印横轴的中心点（除以1000）
print("Bin Centers:", [center // 1000 for center in bin_centers])

with open(args.filename2, 'w') as f:
    for center, bases in zip(bin_centers, bin_bases):
        f.write(f"{center // 1000}\t{bases}\n")
