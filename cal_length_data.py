import numpy as np

# 读取文件
data = np.loadtxt('reads.sequence_length_num.txt', dtype=int)

# 设置bin的宽度和总bin数
bin_width = 8000
total_bins = 372

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

with open('reads.sequence_length_num_out.txt', 'w') as f:
    for center, bases in zip(bin_centers, bin_bases):
        f.write(f"{center // 1000}\t{bases}\n")
