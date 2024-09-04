# -*- coding: utf-8 -*-
# @Time    : 2024/6/25 15:34
# @Author  : jingjing.wu
# @FileName: static_ppi_plot.py
# @Software: PyCharm
# @Email: jingjing.wu@personalbio.cn
import pandas as pd
import igraph as ig
from igraph import Graph, plot
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.gridspec import GridSpec
import sys
import os

def get_subfolders(directory):
    # 列出目录中的所有文件和文件夹
    all_entries = os.listdir(directory)
    # 过滤出所有子文件夹
    subfolders = [entry for entry in all_entries if os.path.isdir(os.path.join(directory, entry))]
    return subfolders

def generate_filenames(directory):
    subfolders = get_subfolders(directory)
    network_filepaths = [os.path.join(directory,f"{subfolder}/",f"{subfolder}_network.xlsx") for subfolder in subfolders]
    attributes_filepaths = [os.path.join(directory,f"{subfolder}/", f"{subfolder}_attributes.xlsx") for subfolder in subfolders]
    subfolder_path = [os.path.join(directory,f"{subfolder}")for subfolder in subfolders]
    return network_filepaths, attributes_filepaths,subfolder_path

def plot_ppi(input_attributes_file, input_network_file, extra_info_file,output_dir):
    # 从文本文件中读取数据
    data = pd.read_excel(input_attributes_file,sheet_name='Sheet 1')
    ppi_data = pd.read_excel(input_network_file,sheet_name='Sheet 1')

    # 从额外信息文件中读取文件名部分
    extra_info = extra_info_file
    # with open(extra_info_file, 'r') as f:
    #     extra_info = f.read().strip()  # 读取并去除首尾空白

    # # 从DataFrame中提取边和权重
    # edges = list(zip(ppi_data['Node1'], ppi_data['Node2']))
    # scores = list(ppi_data['Score'])

    # 设置阈值
    threshold = 700

    # 获取所有节点
    nodes = list(set(ppi_data['Node1']).union(set(ppi_data['Node2'])))

    # 过滤 DataFrame 只保留 weight 大于阈值的边缘
    filtered_df = ppi_data[ppi_data['Score'] >= threshold]
    scores = list(filtered_df['Score'])

    # 创建空图并添加所有节点
    g = ig.Graph(directed=False)
    g.add_vertices(nodes)

    # 添加边缘
    edges = list(zip(filtered_df['Node1'], filtered_df['Node2']))
    g.add_edges(edges)
    g_simplified = g.simplify(multiple=True, loops=True)

    # 设置边缘权重
    g_simplified.es['weight'] = filtered_df['Node1'].tolist()

    # 计算节点的度
    degrees = g_simplified.degree()

    # 将节点度数添加为节点属性
    g_simplified.vs['degree'] = degrees

    # 根据score值调整边的粗细
    edge_weights = [score / max(scores) * 3 for score in scores]  # 归一化并调整粗细

    # 根据度数调整节点大小
    node_sizes = [(degree / max(degrees) * 20 + 10) * 2.2 for degree in degrees]  # 归一化并调整大小
    #node_sizes = degrees * 5000
    # 根据Regulation列调整节点颜色
    regulation_colors = {'up': '#f5780c', 'down': '#0586c0'}  # Morandi红和蓝
    node_colors = [regulation_colors[reg] for reg in data['Regulation']]

    # 可视化图
    layout = g_simplified.layout('fr')  # kk布局

    # 绘制图
    visual_style = {
        "vertex_size": node_sizes,
        "vertex_color": node_colors,
        "vertex_label": g_simplified.vs["name"],
        "vertex_frame_color": '#636363',  # 设置节点边的颜色
        "edge_width": edge_weights,
        "edge_color": '#636363',
        "layout": layout,
        "bbox": (1920, 1440),
        "margin": 10
    }

    # 使用Matplotlib绘制图例
    plt.figure(figsize=(19.2, 14.4))
    ig.plot(g_simplified, **visual_style, target=plt.gca())
    # 添加图例
    legend_elements = [Patch(facecolor='#f5780c', edgecolor='black', label='Up'),
                       Patch(facecolor='#0586c0', edgecolor='black', label='Down')]

    first_legend = plt.legend(handles=legend_elements, loc='upper right', fontsize=25, markerscale=4,labelspacing=0.4,handletextpad=0.2,title_fontsize=25)
    # 调整图例的整体大小和其他属性
    frame = first_legend.get_frame()
    frame.set_facecolor('white')  # 设置图例的背景色
    frame.set_linewidth(1.5)  # 设置图例的边框宽度
    frame.set_edgecolor('white')  # 设置图例的边框颜色
    # 设置图例的整体大小
    plt.setp(first_legend.get_texts(), fontsize='large')  # 设置图例文本的大小
    first_legend.set_title('Regulation', prop={'size': 'large'})  # 设置图例标题的大小
    # 调整图例的位置
    first_legend.set_bbox_to_anchor((1.1,0.7))  # 将图例放置在右上角

    # Add the legend manually to the Axes.
    plt.gca().add_artist(first_legend)

    def some_function(sample_size):
        s = [sample_size_1 * 10 for sample_size_1 in sample_size]
        return s

    # 添加节点大小图例
    for size in range(2, 8, 2):
        sample_size = [(size / max(degrees) * 20 + 10) * 5]
        plt.scatter([], [], c='k', alpha=0.4, s=some_function(sample_size), label=f'    {size}')
    nd_legend = plt.legend(columnspacing=2.5 ,scatterpoints=1, frameon=False, labelspacing=5 / max(degrees) + 1.5 , loc='lower right', fontsize=11,handletextpad=0.1,title_fontsize=25)

    nd_legend.set_title('Size', prop={'size': 'large'})

    nd_legend.set_bbox_to_anchor((1.05, 0.35))  # 将图例放置在右上角

    # 生成文件名
    base_filename = f"{extra_info}.network"
    output_png = f"{output_dir}/{base_filename}.png"
    output_pdf = f"{output_dir}/{base_filename}.pdf"
    output_svg = f"{output_dir}/{base_filename}.svg"

    # 保存和展示图
    plt.savefig(output_png)
    plt.savefig(output_pdf)
    plt.savefig(output_svg)

    return plt


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python static_ppi_plot.py <directory_path>")
        sys.exit(1)

    # 示例用法
    directory_path = sys.argv[1]
    subfolders = get_subfolders(directory_path)
    print(subfolders)
    network_files, attributes_files,subfolder_path = generate_filenames(directory_path)
    print("Network Files:")
    print(network_files)
    print("Attributes Files:")
    print(attributes_files)

    for i,j,k,l in zip(subfolders,network_files,attributes_files,subfolder_path):
        print(i,j,k,l)
        plot_ppi(k, j, i, l)

