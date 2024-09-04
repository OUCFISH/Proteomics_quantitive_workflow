# -*- coding: utf-8 -*-
# @Time    : 2024/6/24 14:50
# @Author  : jingjing.wu
# @FileName: ppi_plot.py
# @Software: PyCharm
# @Email: jingjing.wu@personalbio.cn

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from d3graph import d3graph
import numpy as np
from selenium import webdriver
from selenium.webdriver.chrome.service import Service
from PIL import Image
import time
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
    network_filepaths = [os.path.join(directory,f"{subfolder}/",f"{subfolder}_network_pro.xlsx") for subfolder in subfolders]
    attributes_filepaths = [os.path.join(directory,f"{subfolder}/", f"{subfolder}_attributes.xlsx") for subfolder in subfolders]
    subfolder_path = [os.path.join(directory,f"{subfolder}")for subfolder in subfolders]
    return network_filepaths, attributes_filepaths,subfolder_path


# def generate_size_legend_html(max_degree,target_max):
#     size_legend_html = ""

#     if max_degree <= 10:
#         sizes = [2, 4, 6]
#     else:
#         # 使用分割方法取其中的三个大小
#         step = max_degree // 4
#         sizes = [step, step * 2, step * 3]
#         # 计算当前最大值

#     current_max = max(sizes)
    
#     # 计算缩放因子
#     scale_factor = target_max / current_max
    
#     # 进行缩放
#     scaled_sizes = [size * scale_factor for size in sizes]

#     for size in scaled_sizes:
#         size_legend_html += f"""
#             <div style="background-color: black; width: {size * 2}px; height: {size * 2}px; display: inline-block; border-radius: {size}px; border: 1px solid black; margin-bottom: 5px;"></div>
#             <span style="display: inline-block; width: 20px; text-align: center;">{size}</span>
#             <br>
#         """

#     return size_legend_html
def generate_size_legend_html(max_degree, target_max):
    size_legend_html = ""

    if max_degree <= 10:
        sizes = [2, 4, 6]
    else:
        # 使用分割方法取其中的三个大小
        step = max_degree // 4
        sizes = [step, step * 2, step * 3]

    # 计算当前最大值
    current_max = max(sizes)
    
    # 计算缩放因子
    scale_factor = target_max / current_max
    
    # 进行缩放
    scaled_sizes = [size * scale_factor for size in sizes]

    for original_size, size in zip(sizes, scaled_sizes):
        size_legend_html += f"""
            <div style="background-color: black; width: {size * 2}px; height: {size * 2}px; display: inline-block; border-radius: {size}px; border: 1px solid black; margin-bottom: 5px;"></div>
            <span style="display: inline-block; width: 50px; text-align: center;">{original_size}</span>
            <br>
        """

    return size_legend_html

def ppi_plot(input_attributes_file, input_network_file, extra_info_file,output_dir):
    # 从文本文件中读取数据
    data = pd.read_excel(input_attributes_file)
    ppi_data = pd.read_excel(input_network_file)

    # 从额外信息文件中读取文件名部分
    extra_info = extra_info_file

    # 创建颜色映射
    # 假设 'Regulation' 列的值是 'up' 和 'down'
    up_color = '#f5780c'  # 红色
    down_color = '#0586c0'  # 蓝色

    # 根据 Regulation 值映射颜色
    data['Color'] = data['Regulation'].apply(lambda x: down_color if x == 'down' else up_color)

    # 提取边列表
    edges = ppi_data[['Node1', 'Node2', 'Score']]
    # 创建邻接矩阵
    nodes = pd.concat([edges['Node1'], edges['Node2']]).unique()
    adj_matrix = pd.DataFrame(0, index=nodes, columns=nodes)

    for _, row in edges.iterrows():
        adj_matrix.at[row['Node1'], row['Node2']] = row['Score']
        # adj_matrix.at[row['Node2'], row['Node1']] = row['Score']  # 无向图

    # 初始化 D3Graph
    d3 = d3graph(collision=0.2,charge= 200)

    # 构建图
    d3.graph(adj_matrix)

    # 创建节点颜色列表
    node_colors = []
    for node in nodes:
        if node in data['Protein'].values:
            color = data.loc[data['Protein'] == node, 'Color'].values[0]
        else:
            color = '#808080'  # 灰色
        node_colors.append(color)

    node_colors = data['Color'].values

    data_by_prot = data.sort_values(by='Protein')
    # #创建节点形状列表
    # markers = []
    # for node in nodes:
    #     if node in data['GeneName'].values:
    #         if data.loc[data['GeneName'] == node, 'Log2fc'].values[0] > 0:
    #             marker = 'rect'
    #         else:
    #             marker = 'circle'
    #     else:
    #         marker = 'square'  # 正方形
    #     markers.append(marker)
    # print(markers)
    # 设置节点属性
    # Set node properties
    tooltip = '\nGeneName: ' + data_by_prot['GeneName'] + '\nProtein id: ' + data_by_prot[
        'Protein'] + '\n-log10(Pvalue): ' + data_by_prot['-log10(Pvalue)'].values.astype(str) + '\nDegree: ' + \
              data_by_prot['Degree'].values.astype(str) + '\nLog2fc: ' + data_by_prot['Log2fc'].values.astype(str)
    tooltip = tooltip.values

    d3.set_node_properties(tooltip=tooltip, label=data_by_prot['Protein'], color=data_by_prot['Color'].values,
                           fontcolor='black', marker='circle', edge_color='#808080', size='degree')

    # 设置边的参数
    d3.set_edge_properties(directed=False, edge_distance=50)



    # 生成文件名
    base_filename = f"{extra_info}.network"
    output_html = f"{output_dir}/{base_filename}.html"
    # 将图形保存为HTML文件
    d3.show(filepath=output_html)
    #设置图例
    color_legend_html = f"""
        <div style="display: flex; align-items: center; justify-content: center; margin-bottom: 10px;">
            <div style="background-color: {up_color}; width: 30px; height: 30px; border-radius: 50%; margin-right: 10px;"></div>
            <span>Up&emsp;&nbsp;</span>
        </div>
        <div style="display: flex; align-items: center; justify-content: center; margin-bottom: 10px;">
            <div style="background-color: {down_color}; width: 30px; height: 30px; border-radius: 50%; margin-right: 10px;"></div>
            <span>Down</span>
        </div>
    """

    size_legend_html = ""

    # 将节点大小匹配最大度数  小于10 用2，4，6 大于10 用分割方法 取其中的三个
    # 示例使用

    max_value = data["Degree"].max()
    size_legend_html = generate_size_legend_html(max_value,20)
    # for size in [2, 4, 6, 8, 10]:  # 示例大小
    #     size_legend_html += f"""
    #         <div style="background-color: black; width: {size * 2}px; height: {size * 2}px; display: inline-block; border-radius: {size}px; border: 1px solid black; margin-bottom: 5px;"></div>
    #         <span style="display: inline-block; width: 20px; text-align: center;">{size}</span>
    #         <br>
    #     """
    with open(output_html, "a") as f:
        f.write(f"""
            <style>
                #legendContainer {{
                    position: absolute;
                    right: 20px;
                    top: 100px;
                    padding: 5px;
                    background-color: white;
                    border: 1px solid black;
                    border-radius: 5px;
                    box-shadow: 0px 0px 10px rgba(0, 0, 0, 0.1);
                    width: 100px;
                }}
                #sizeLegend, #colorLegend {{
                    margin-bottom: 20px;
                    text-align: center;
                }}
            </style>
            <div id="legendContainer">
                <div id="colorLegend">
                    <h3>Regulation</h3>
                    {color_legend_html}
                </div>
                <div id="sizeLegend">
                    <h3>Degree</h3>
                    {size_legend_html}
                </div>
            </div>
        """)

        #删除广告脚本
        # 读取文件内容
        text_to_delete = "<script async src='https://media.ethicalads.io/media/client/ethicalads.min.js'></script>"

        with open(output_html, 'r', encoding='utf-8') as file:
            lines = file.readlines()

        # 检查前20行并删除包含指定文本的行
        for i in range(min(20, len(lines))):
            if text_to_delete in lines[i]:
                lines.pop(i)
                break  # 假设每次只删除一行，找到并删除后即可退出循环
        # 将修改后的内容写回文件
        with open(output_html, 'w', encoding='utf-8') as file:
            file.writelines(lines)

    return

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python static_ppi_plot.py <directory_path>")
        sys.exit(1)

    # 示例用法
    directory_path = sys.argv[1]
    subfolders = get_subfolders(directory_path)
    network_files, attributes_files,subfolder_path = generate_filenames(directory_path)
    for i,j,k,l in zip(subfolders,network_files,attributes_files,subfolder_path):
        print('开始绘图')
        ppi_plot(k, j, i, l)

