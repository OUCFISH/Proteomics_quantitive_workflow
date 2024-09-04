# -*- coding: utf-8 -*-
# @Time    : 2024/8/19 15:20
# @Author  : chen.zy
# @FileName: create_qmd.py
# @Software: PyCharm
# @Email: chenzy@personalbio.cn


import os
import paddleocr
from paddleocr import PaddleOCR
import pandas as pd
import re
import argparse
import time



def create_pro_info():
    def ocr(pdf_file):
        texts = []
        PAGE_NUM = 1
        ocr = PaddleOCR(use_angle_cls=True, lang="ch", page_num=PAGE_NUM)
        result = ocr.ocr(pdf_file, cls=True)
        for idx in range(len(result)):
            res = result[idx]
            txts = [line[1][0] for line in res]
        texts += txts

        all_text = " ".join(texts)

        info_dict = {}
        info_dict["项目编号"] = os.path.basename(pdf_file).replace(".pdf", "")
        info_dict["合同号"] = re.findall("合同编号[：:;]{1,2}\s*(\S*)", all_text)[0] if re.findall("合同编号[：:;]{1,2}\s*(\S*)",
                                                                                           all_text) else ""
        projec_type = args.type
        if projec_type == "DIA":
            info_dict["项目名称"] = "Library-free DIA 蛋白质组学分析"
        elif projec_type == "TMT":
            info_dict["项目名称"] = "TMT 蛋白质组学分析"
        elif projec_type == "LFQ":
            info_dict["项目名称"] = "Label free 蛋白质组学分析"
        elif projec_type == "PQA":
            info_dict["项目名称"] = "Shotgun 蛋白鉴定"
        else:
            raise Exception("项目类型错误")
        info_dict["项目类别"] = "医学蛋白质组学分析"
        info_dict["项目类型"] = "标准分析"
        info_dict["完成日期"] = time.strftime("%Y/%m/%d", time.localtime())
        info_dict["客户姓名"] = re.findall("客户姓名[：:;]{1,2}\s*(\S*)", all_text)[0] if re.findall("客户姓名[：:;]{1,2}\s*(\S*)",
                                                                                            all_text) else ""
        info_dict["委托单位"] = re.findall("客户单位[：:;]{1,2}\s*(\S*)", all_text)[0] if re.findall("客户单位[：:;]{1,2}\s*(\S*)",
                                                                                            all_text) else ""
        return info_dict

    pdf_file = [fi for fi in os.listdir(r"./result/1_Sampleinfo") if fi.endswith(".pdf")][0]
    infos = ocr(os.path.join(r"./result/1_Sampleinfo", pdf_file))
    with open("./table/project_info.txt", "w", encoding="utf-8") as f:
        for k in ["项目编号", "合同号", "项目名称", "项目类别", "项目类型", "完成日期", "客户姓名", "委托单位"]:
            f.write(f"{k}\t{infos[k]}\n")

def create_sample_info():
    group_df = pd.read_excel("./temp/group.xlsx")
    groups = zip(group_df["SampleName"].tolist(), group_df["Group"].tolist())
    with open("./config/default.conf", "r", encoding="utf-8") as f:
        for line in f:
            if "organism_Scientific_name" in line:
                species = line.split("<-")[1].strip().split(" ")[0].replace("\"", "")
                break
    if args.type == 'TMT':
        xlsx = [fi for fi in os.listdir(r'./result/1_Sampleinfo') if '对照' in fi][0]
        try:
            df = pd.read_excel(rf'./result/1_Sampleinfo/{xlsx}')
        except:
            raise Exception("没有找到对照组信息")
        with open('./table/sample_info.txt', 'w', encoding='utf-8') as f:
            f.write('Sample Names\tGroup Names\tTMT Label\tSpecies\n')
            for sample, group in groups:
                f.write(f"{sample}\t{group}\t{df[df['报告内样本名称']==sample]['标记'].values[0]}\t{species}\n")
    else:
        with open("./table/sample_info.txt", "w", encoding="utf-8") as f:
            f.write("Sample Names\tGroup Names\tSpecies")
            for sample, group in groups:
                f.write(f"\n{sample}\t{group}\t{species}")

def create_diff_analyse_result():
    with open(r"./config/default.conf", "r", encoding="utf-8") as conf:
        for line in conf:
            if line.startswith("Fold_Change1"):
                FC = line.split(" ")[-1].strip()
            if line.startswith("pvalue"):
                pvalue = line.split(" ")[-1].strip()
                break
    contrast = pd.read_excel("./temp/contrast.xlsx")
    compares = contrast["contrast"].tolist()
    folder = [directory for directory in os.listdir(r"./result") if "Data_Summary" in directory]
    if len(folder) == 0:
        raise Exception("没有找到Data_Summary文件夹")
    file_path = rf"./result/{folder[0]}"
    xlsx = [fi for fi in os.listdir(file_path) if
            fi.endswith(".xlsx") and "SummaryStat" in fi and "$" not in fi]
    df = pd.read_excel(os.path.join(file_path, xlsx[0]))
    with open(r"./table/diff_analyse_result.txt", "w", encoding="utf-8") as f:
        f.write(f"Statistical methods\tFold change and student’s t test\tFold change > {FC}\tP value < {pvalue}\n")
        f.write("Comparisons\tUp Regulation\tDown Regulation\tAll\n")
        for compare in compares:
            up = df[df[f"{compare}.Regulation"] == "Up"].shape[0]
            down = df[df[f"{compare}.Regulation"] == "Down"].shape[0]
            f.write(
                f"{compare}\t{up}\t{down}\t{up + down}\n")
    with open(r"./table/diff_exp.txt", "w", encoding="utf-8", newline="") as f:
        f.write("Control_vs_Treat\tUp Regulation\tDown Regulation\tTotal\n")
        for compare in compares:
            up = df[df[f"{compare}.Regulation"] == "Up"].shape[0]
            down = df[df[f"{compare}.Regulation"] == "Down"].shape[0]
            f.write(f"{compare}\t{up}\t{down}\t{up + down}\n")

def create_variable():
    with open(r"./config/default.conf", "r", encoding="utf-8") as conf:
        for line in conf:
            if line.startswith("Fold_Change1"):
                FC = line.split(" ")[-1].strip()
            if line.startswith("pvalue"):
                pvalue = line.split(" ")[-1].strip()
                break
    fasta = ""
    directory = [folder for folder in os.listdir(r"./result") if "Attachments" in folder][0]
    for file in os.listdir(f"./result/{directory}"):
        if file.endswith(".fasta"):
            fasta = f"本次数据库为：{file}"
            break
    with open(r"./variable.txt", "w", encoding="utf-8") as f:
        f.write(f"FC\t{FC}\n")
        f.write(f"pvalue\t{pvalue}\n")
        f.write(f"fasta\t{fasta}\n")
        f.write(f"AA\t{'AA' if args.advanced_analyse else 'noAA'}\n")

def create_search_result():
    folder = [directory for directory in os.listdir(r"./result") if "Data_Summary" in directory]
    if len(folder) == 0:
        raise Exception("没有找到Data_Summary文件夹")
    file_path = rf"./result/{folder[0]}"
    Peptides_file = [fi for fi in os.listdir(file_path) if fi.endswith(".xlsx") and "Peptide" in fi and "$" not in fi][0]
    df = pd.read_excel(os.path.join(file_path, Peptides_file))
    Peptides = len(df)
    Identified_file = [fi for fi in os.listdir(file_path) if
                       fi.endswith(".xlsx") and "ProteinGroup" in fi and "$" not in fi][0]
    df = pd.read_excel(os.path.join(file_path, Identified_file))
    Identified = len(df)
    Quantified_file = [fi for fi in os.listdir(file_path) if
                       fi.endswith(".xlsx") and "SummaryStat" in fi and "$" not in fi][0]
    df = pd.read_excel(os.path.join(file_path, Quantified_file))
    Quantified = len(df)
    with open("./table/search_result.txt", "w", encoding="utf-8") as f:
        f.write("Result\tPeptides\tIdentified proteins\tQuantified proteins\n")
        f.write(f"Uniprot database\t{Peptides}\t{Identified}\t{Quantified}")

def create_annotations():
    folder = [directory for directory in os.listdir(r"./result") if "Data_Summary" in directory]
    if len(folder) == 0:
        raise Exception("没有找到Data_Summary文件夹")
    file_path = rf"./result/{folder[0]}"
    xlsx = [fi for fi in os.listdir(file_path) if
            fi.endswith(".xlsx") and "SummaryStat" in fi and "$" not in fi][0]
    df = pd.read_excel(os.path.join(file_path, xlsx))
    total = len(df)
    annotations = {}
    try:
        annotations["GO"] = df[df["Gene.Ontology.(GO)"].notnull()].shape[0]
    except:
        pass
    try:
        annotations["KEGG"] = df[df["pathway.id"].notnull()].shape[0]
    except:
        pass
    try:
        annotations["Reactome"] = df[df["Reactome"].notnull()].shape[0]
    except:
        pass
    try:
        annotations["Domain"] = df[df["InterPro.id"].notnull()].shape[0]
    except:
        pass
    try:
        annotations["TF"] = df[df["TF"].notnull()].shape[0]
    except:
        pass
    try:
        annotations["Kinase"] = df[df["kinase"]].shape[0]
    except:
        pass
    try:
        annotations["Subcellular location"] = df[df["subcellular_location"].notnull()].shape[0]
    except:
        pass
    try:
        annotations["eggNOG"] = df[df["eggNOG"].notnull()].shape[0]
    except:
        pass
    try:
        annotations["DO"] = df[df["DO"].notnull()].shape[0]
    except:
        pass
    try:
        annotations["Wikipathway"] = df[df["Wikipathway"].notnull()].shape[0]
    except:
        pass
    with open(r"./table/annotations.txt", "w", encoding="utf-8", newline="") as f:
        f.write("Database\tAnnotated\tPercent\n")
        for key, value in annotations.items():
            f.write(f"{key}\t{value}\t{round(value / total * 100, 2):.2f}%\n")

def style():
    def header_format(info):
        if info == 'base_info':
            bg_color = '#0586C0'
        elif info == 'sample_info':
            bg_color = '#F5780C'
        elif info == 'Regulation':
            bg_color = '#8064A2'
        elif info == 'annotation':
            bg_color = '#595959'
        else:
            bg_color = '#0586C0'
        return {
            'font_name': 'calibri',
            'bold': False,  # 设置字体加粗
            'font_size': 11,  # 设置字体大小
            'bg_color': bg_color,  # 设置背景颜色
            'color': 'white',  # 设置字体颜色
            'align': 'left',  # 设置文本居中对齐
            'valign': 'vcenter',  # 设置垂直居中对齐
            'border_color': '#CCCCCC',  # 设置边框颜色
            'border': 1,  # 边框宽度
            'top': 1,  # 上边框
            'left': 1,  # 左边框
            'right': 1,  # 右边框
            'bottom': 1  # 底边框
        }

    def cell_format(info):
        if info == 'Down':
            color = '#0070C0'
        elif info == 'Up':
            color = '#FF0000'
        else:
            color = 'black'
        return {
            'font_name': 'calibri',
            'bold': False,  # 设置字体加粗
            'font_size': 10,  # 设置字体大小
            'bg_color': '#FFFFFF',  # 设置背景颜色
            'color': color,  # 设置字体颜色
            'border_color': '#CCCCCC',  # 设置边框颜色
            'border': 1,  # 边框宽度
            'top': 1,  # 上边框
            'left': 1,  # 左边框
            'right': 1,  # 右边框
            'bottom': 1  # 底边框
        }
    for directory in os.listdir(r"./result"):
        if "Data_Summary" in directory:
            try:
                xlsx = [fi for fi in os.listdir(rf"./result/{directory}") if
                        fi.endswith(".xlsx") and "SummaryStat" in fi and "$" not in fi][0]
                xlsx_path = rf'./result/{directory}/{xlsx}'
                if os.path.exists(xlsx_path.replace('.xlsx', '_styled.xlsx')):
                    return
            except:
                raise Exception("没有找到SummaryStat文件")
    df = pd.read_excel(xlsx_path)
    writer = pd.ExcelWriter(xlsx_path.replace('.xlsx', '_styled.xlsx'), engine='xlsxwriter')
    workbook = writer.book
    ws = workbook.add_worksheet('Sheet 1')
    columns = df.columns.tolist()
    info = 'base_info'
    for i in range(len(columns)):
        if i > 0 and (columns[i - 1] == 'Coverage' or columns[i - 1].endswith('.Regulation')):
            info = 'sample_info'
        if columns[i].endswith('.Regulation'):
            info = 'Regulation'
        if columns[i] == ('Gene.Ontology.(GO)'):
            info = 'annotation'
        ws.write(0, i, columns[i], workbook.add_format(header_format(info)))
    ws.set_row(0, height=20)
    df = df.fillna('')
    rows, cols = df.shape
    for row in range(rows):
        for col in range(cols):
            value = df.iat[row, col]
            if col == cols - 2:
                ws.write(row + 1, col, str(value), workbook.add_format(cell_format(value)))
            else:
                ws.write(row + 1, col, value, workbook.add_format(cell_format(value)))
    writer._save()

def get_content(file_path: str) -> str:
    """
    获取文件内容
    :param file_path: 文件路径
    :return: 文件内容
    """
    with open(file_path, "r", encoding="utf-8") as f:
        content = f.read()
    return content + "\n"

def get_file_path(path: str, keyword: str, Format: str) -> str:
    """
    获取文件路径
    :param path: 路径
    :param keyword: 关键字
    :param Format: 格式
    :return: 文件路径
    """
    for file in os.listdir(path):
        if not file.endswith(Format):
            continue
        if keyword in file:
            abs_path = os.path.join(path, file)
            return abs_path.replace(os.getcwd(), ".").replace("\\", "/")
    return ""

def parse(md: str, path: str) -> str:
    """
    解析md
    :param md: md文件
    :param path: 路径
    :return: 解析后的文本
    """
    texts = ""
    file_list = []
    with open(md, encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if line.startswith("![") or line.startswith("&["):
                text = re.findall("\((.*)\)", line, re.S)[0]
                name = re.findall("\[(.*)\]", line, re.S)[0]
                variable = {x.split(":")[0].strip(): x.split(":")[1].strip() for x in text.split(",")}
                variable["name"] = name
                variable["type"] = "single"
                file_list.append(variable)
                if "depth" in variable.keys():
                    file_path = ""
                    for folder in os.listdir(path):
                        if os.path.isdir(os.path.join(path, folder)):
                            file_path = get_file_path(os.path.join(path, folder), variable["keyword"], variable["format"])
                            if file_path:
                                break
                else:
                    file_path = get_file_path(path, variable["keyword"], variable["format"])
                if not file_path:
                    return ""
                texts += line.replace(text, file_path) + "\n"
            elif line.startswith("[link]"):
                text = re.findall("\((.*)\)", line, re.S)[0]
                name = re.findall("\[(.*?)\]", line, re.S)[1]
                if name in ["3D PCA 图", "蛋白网络互作分析动态图"]:
                    variable = {x.split(":")[0].strip(): x.split(":")[1].strip() for x in text.split(",")}
                    if "depth" in variable.keys():
                        for folder in os.listdir(path):
                            file_path = get_file_path(os.path.join(path, folder), variable["keyword"],
                                                      variable["format"])
                            if file_path:
                                break
                    else:
                        file_path = get_file_path(path, variable["keyword"], variable["format"])
                    texts += line.replace(text, file_path) + "\n"
                else:
                    texts += line + "\n"
            elif line.startswith("imagebox["):
                text = re.findall("\((.*)\)", line, re.S)[0]
                name = re.findall("\[(.*)\]", line, re.S)[0]
                variable = {x.split(":")[0].strip(): x.split(":")[1].strip() for x in text.split(",")}
                variable["name"] = name
                variable["type"] = "multiple"
                file_list.append(variable)
                if "depth" in variable.keys():
                    if not any([os.path.isdir(os.path.join(subfolder, f)) for folder in os.listdir(path)
                                if os.path.isdir(subfolder := os.path.join(path, folder))
                                for f in os.listdir(subfolder)]):
                        return ""
                elif not any([get_file_path(os.path.join(path, folder), variable["keyword"], variable["format"])
                            for folder in os.listdir(path) if os.path.isdir(os.path.join(path, folder))]):
                    return ""
                texts += line.replace(")", f",path:{path})") + "\n"
            elif line.startswith("[result_file]"):
                texts += "[result_file_start]\n"
                for d in file_list:
                    if "图" in d["name"]:
                        format_list = ["png", "pdf", "html"]
                    elif "表" in d["name"]:
                        format_list = ["xlsx"]
                    else:
                        format_list = ["html"]
                    if d["type"] == "single":
                        if "depth" in d.keys():
                            for folder in os.listdir(path):
                                if os.path.isdir(os.path.join(path, folder)):
                                    for file in os.listdir(os.path.join(path, folder)):
                                        if d["keyword"] in file and file.split(".")[-1] in format_list:
                                            file_path = os.path.join(path, folder, file).replace("\\", "/")
                                            texts += f"[{folder} {d['name']}（{file.split('.')[-1]}格式）]({file_path.replace(os.getcwd(), '.')})\n"
                        else:
                            for file in os.listdir(path):
                                if d["keyword"] in file and file.split(".")[-1] in format_list:
                                    file_path = os.path.join(path, file).replace("\\", "/")
                                    texts += f"[{d['name']}（{file.split('.')[-1]}格式）]({file_path.replace(os.getcwd(), '.')})\n"
                    elif "depth" in d.keys():  # type为multiple
                        for folder in os.listdir(path):
                            if os.path.isdir(os.path.join(path, folder)):
                                for f in os.listdir(os.path.join(path, folder)):
                                    if os.path.isdir(os.path.join(path, folder, f)):
                                        folder_path = os.path.join(path, folder, f).replace(os.getcwd(), ".").replace("\\", "/")
                                        texts += f"[{folder} {d['name']}]({folder_path})\n"
                    else:
                        for folder in os.listdir(path):
                            if os.path.isdir(os.path.join(path, folder)):
                                for file in os.listdir(os.path.join(path, folder)):
                                    if d["keyword"] in file and file.split(".")[-1] in format_list:
                                        file_path = os.path.join(path, folder, file).replace(os.getcwd(), ".").replace("\\", "/")
                                        texts += f"[{folder} {d['name']}（{file.split('.')[-1]}格式）]({file_path})\n"
                texts += "[result_file_end]\n"
            else:
                texts += line + "\n"
    return texts



def main(project_type):
    """
    主函数
    """
    # 创建路径
    if not os.path.exists("./table"):
        os.mkdir("./table")
    if not os.path.exists("./image"):
        os.mkdir("./image")
    # 创建参数文件
    create_variable()
    # 生成md
    with open(f"./prot_report_{project_type}.md", 'w', newline='', encoding='utf-8') as f:

        # DIA、TMT、LFQ
        if project_type in ["DIA", "TMT", "LFQ"]:

            # 生成表格
            if not os.path.exists("./table/project_info.txt"):
                create_pro_info()
            if not os.path.exists("./table/sample_info.txt"):
                create_sample_info()
            if not os.path.exists("./table/diff_analyse_result.txt"):
                create_diff_analyse_result()
            if not os.path.exists("./table/search_result.txt"):
                create_search_result()
            if not os.path.exists("./table/annotations.txt"):
                create_annotations()
            style()
            if project_type == "DIA":
                f.write("# DIA 蛋白组\n")
                f.write(get_content("./temp/1.项目整体流程概况_DIA.md") + "\n")
            elif project_type == "TMT":
                f.write("# TMT 蛋白组\n")
                f.write(get_content("./temp/1.项目整体流程概况_TMT.md") + "\n")
            elif project_type == "LFQ":
                f.write("# Labelfree 蛋白组\n")
                f.write(get_content("./temp/1.项目整体流程概况_LFQ.md") + "\n")
            f.write(get_content("./temp/2.项目结果概览.md") + "\n")

            # 3_Quality_Control
            title_status = 0
            for directory in os.listdir("./result"):
                if "Quality_Control" in directory:

                    # 3.1_Quality_Control
                    sub_title_status = 0
                    for directory_1 in os.listdir(f"./result/{directory}"):
                        if "Quality_Control" in directory_1:
                            for directory_2 in os.listdir(f"./result/{directory}/{directory_1}"):
                                if "peptides_length_distribution" in directory_2:
                                    content = parse("./temp/3.1.1.肽段长度分布.md", f"./result/{directory}/{directory_1}/{directory_2}")
                                    if content:
                                        if not title_status:
                                            f.write("## 实验结果\n")
                                            title_status = 1
                                        if not sub_title_status:
                                            f.write("### 质控\n")
                                            sub_title_status = 1
                                        f.write(content)
                                    break
                            for directory_2 in os.listdir(f"./result/{directory}/{directory_1}"):
                                if "protein_molecular_distribution" in directory_2:
                                    content = parse("./temp/3.1.2.蛋白分子量分布.md", f"./result/{directory}/{directory_1}/{directory_2}")
                                    if content:
                                        if not title_status:
                                            f.write("## 实验结果\n")
                                            title_status = 1
                                        if not sub_title_status:
                                            f.write("### 质控\n")
                                            sub_title_status = 1
                                        f.write(content)
                                    break
                            for directory_2 in os.listdir(f"./result/{directory}/{directory_1}"):
                                if "proteins_coverage_distribution" in directory_2:
                                    content = parse("./temp/3.1.3.蛋白覆盖度分布.md", f"./result/{directory}/{directory_1}/{directory_2}")
                                    if content:
                                        if not title_status:
                                            f.write("## 实验结果\n")
                                            title_status = 1
                                        if not sub_title_status:
                                            f.write("### 质控\n")
                                            sub_title_status = 1
                                        f.write(content)
                                    break
                            for directory_2 in os.listdir(f"./result/{directory}/{directory_1}"):
                                if "proteins_isoelectric_distribution" in directory_2:
                                    content = parse("./temp/3.1.4.等电点分布.md", f"./result/{directory}/{directory_1}/{directory_2}")
                                    if content:
                                        if not title_status:
                                            f.write("## 实验结果\n")
                                            title_status = 1
                                        if not sub_title_status:
                                            f.write("### 质控\n")
                                            sub_title_status = 1
                                        f.write(content)
                                    break
                            for directory_2 in os.listdir(f"./result/{directory}/{directory_1}"):
                                if "IRT" in directory_2:
                                    content = parse("./temp/3.1.5.IRT保留时间.md", f"./result/{directory}/{directory_1}/{directory_2}")
                                    if content:
                                        if not title_status:
                                            f.write("## 实验结果\n")
                                            title_status = 1
                                        if not sub_title_status:
                                            f.write("### 质控\n")
                                            sub_title_status = 1
                                        f.write(content)
                                    break
                            break

                    # 3.2_Quantitative_Statistics
                    sub_title_status = 0
                    for directory_1 in os.listdir(f"./result/{directory}"):
                        if "Quantitative_Statistics" in directory_1:
                            for directory_2 in os.listdir(f"./result/{directory}/{directory_1}"):
                                if "Expression_Statistics" in directory_2:
                                    content = parse("./temp/3.2.1.蛋白鉴定数量统计.md", f"./result/{directory}/{directory_1}/{directory_2}")
                                    if content:
                                        if not title_status:
                                            f.write("## 实验结果\n")
                                            title_status = 1
                                        if not sub_title_status:
                                            f.write("### 样本统计分析\n")
                                            sub_title_status = 1
                                        f.write(content)
                                    break
                            for directory_2 in os.listdir(f"./result/{directory}/{directory_1}"):
                                if "Expression_Repeat" in directory_2:
                                    content = parse("./temp/3.2.2.样本鉴定重复性分析.md", f"./result/{directory}/{directory_1}/{directory_2}")
                                    if content:
                                        if not title_status:
                                            f.write("## 实验结果\n")
                                            title_status = 1
                                        if not sub_title_status:
                                            f.write("### 样本统计分析\n")
                                            sub_title_status = 1
                                        f.write(content)
                                    break
                            for directory_2 in os.listdir(f"./result/{directory}/{directory_1}"):
                                if "Expression_Violin" in directory_2:
                                    content = parse("./temp/3.2.3.样本定量值分析.md", f"./result/{directory}/{directory_1}/{directory_2}")
                                    if content:
                                        if not title_status:
                                            f.write("## 实验结果\n")
                                            title_status = 1
                                        if not sub_title_status:
                                            f.write("### 样本统计分析\n")
                                            sub_title_status = 1
                                        f.write(content)
                                    break
                            for directory_2 in os.listdir(f"./result/{directory}/{directory_1}"):
                                if "RSD" in directory_2:
                                    content = parse("./temp/3.2.4.RSD分布.md", f"./result/{directory}/{directory_1}/{directory_2}")
                                    if content:
                                        if not title_status:
                                            f.write("## 实验结果\n")
                                            title_status = 1
                                        if not sub_title_status:
                                            f.write("### 样本统计分析\n")
                                            sub_title_status = 1
                                        f.write(content)
                                    break
                            for directory_2 in os.listdir(f"./result/{directory}/{directory_1}"):
                                if "Expression_Scatter" in directory_2:
                                    content = parse("./temp/3.2.5.丰度分布散点图.md", f"./result/{directory}/{directory_1}/{directory_2}")
                                    if content:
                                        if not title_status:
                                            f.write("## 实验结果\n")
                                            title_status = 1
                                        if not sub_title_status:
                                            f.write("### 样本统计分析\n")
                                            sub_title_status = 1
                                        f.write(content)
                                    break
                            for directory_2 in os.listdir(f"./result/{directory}/{directory_1}"):
                                if "Correlation_Analysis" in directory_2:
                                    content = parse("./temp/3.2.6.相关性分析.md", f"./result/{directory}/{directory_1}/{directory_2}")
                                    if content:
                                        if not title_status:
                                            f.write("## 实验结果\n")
                                            title_status = 1
                                        if not sub_title_status:
                                            f.write("### 样本统计分析\n")
                                            sub_title_status = 1
                                        f.write(content)
                                    break
                            for directory_2 in os.listdir(f"./result/{directory}/{directory_1}"):
                                if "Global_Heatmap" in directory_2:
                                    content = parse("./temp/3.2.7.全局热图分析.md", f"./result/{directory}/{directory_1}/{directory_2}")
                                    if content:
                                        if not title_status:
                                            f.write("## 实验结果\n")
                                            title_status = 1
                                        if not sub_title_status:
                                            f.write("### 样本统计分析\n")
                                            sub_title_status = 1
                                        f.write(content)
                                    break
                            for directory_2 in os.listdir(f"./result/{directory}/{directory_1}"):
                                if "PCA" in directory_2:
                                    content = parse("./temp/3.2.8.主成分分析.md", f"./result/{directory}/{directory_1}/{directory_2}")
                                    if content:
                                        if not title_status:
                                            f.write("## 实验结果\n")
                                            title_status = 1
                                        if not sub_title_status:
                                            f.write("### 样本统计分析\n")
                                            sub_title_status = 1
                                        f.write(content)
                                    break
                    break

            # 4_Diff_Expressed
            title_status = 0
            for directory in os.listdir("./result"):
                if "Diff_Expressed" in directory:
                    for directory_1 in os.listdir(f"./result/{directory}"):
                        if "DiffStats" in directory_1:
                            content = parse("./temp/4.1.差异表达结果统计.md", f"./result/{directory}/{directory_1}")
                            if content:
                                if not title_status:
                                    f.write("## 表达差异分析\n")
                                    title_status = 1
                            f.write(content)
                            break
                    for directory_1 in os.listdir(f"./result/{directory}"):
                        if "Diff_Repeat" in directory_1:
                            content = parse("./temp/4.2.组间差异蛋白比较分析.md", f"./result/{directory}/{directory_1}")
                            if content:
                                if not title_status:
                                    f.write("## 表达差异分析\n")
                                    title_status = 1
                            f.write(content)
                            break
                    for directory_1 in os.listdir(f"./result/{directory}"):
                        if "Volcano" in directory_1:
                            content = parse("./temp/4.3.差异蛋白火山图.md", f"./result/{directory}/{directory_1}")
                            if content:
                                if not title_status:
                                    f.write("## 表达差异分析\n")
                                    title_status = 1
                            f.write(content)
                            break
                    for directory_1 in os.listdir(f"./result/{directory}"):
                        if "Heatmap" in directory_1:
                            content = parse("./temp/4.4.差异蛋白聚类热图.md", f"./result/{directory}/{directory_1}")
                            if content:
                                if not title_status:
                                    f.write("## 表达差异分析\n")
                                    title_status = 1
                            f.write(content)
                            break
                    for directory_1 in os.listdir(f"./result/{directory}"):
                        if "Trend" in directory_1:
                            content = parse("./temp/4.5.差异蛋白趋势分析.md", f"./result/{directory}/{directory_1}")
                            if content:
                                if not title_status:
                                    f.write("## 表达差异分析\n")
                                    title_status = 1
                            f.write(content)
                            break
                    for directory_1 in os.listdir(f"./result/{directory}"):
                        if "Radar" in directory_1:
                            content = parse("./temp/4.6.差异蛋白雷达图.md", f"./result/{directory}/{directory_1}")
                            if content:
                                if not title_status:
                                    f.write("## 表达差异分析\n")
                                    title_status = 1
                            f.write(content)
                            break
                    break

            db_intro = []
            # 5_Enrichment
            for directory in os.listdir("./result"):
                if "Enrichment" in directory:
                    f.write("## 富集分析\n")

                    # 5.1_GO
                    for directory_1 in os.listdir(f"./result/{directory}"):
                        if "GO" in directory_1:
                            content = parse("./temp/5.1.GO富集分析.md", f"./result/{directory}/{directory_1}")
                            if content:
                                db_intro.append("GO")
                                f.write(content)
                                f.write(parse("./temp/5.1.1.GO注释统计柱状图.md", f"./result/{directory}/{directory_1}"))
                                f.write(parse("./temp/5.1.2.差异蛋白GO注释统计柱状图.md", f"./result/{directory}/{directory_1}"))
                                f.write(parse("./temp/5.1.3.差异蛋白GO富集柱状图.md", f"./result/{directory}/{directory_1}"))
                                f.write(parse("./temp/5.1.4.差异蛋白GO富集因子图.md", f"./result/{directory}/{directory_1}"))
                                f.write(parse("./temp/5.1.5.差异蛋白GO富集圈图.md", f"./result/{directory}/{directory_1}"))
                                f.write(parse("./temp/5.1.6.差异蛋白GO富集网络图.md", f"./result/{directory}/{directory_1}"))
                                f.write(parse("./temp/5.1.7.差异蛋白GO富集和弦图.md", f"./result/{directory}/{directory_1}"))
                            break

                    # 5.2_KEGG
                    for directory_1 in os.listdir(f"./result/{directory}"):
                        if "KEGG" in directory_1:
                            content = parse("./temp/5.2.KEGG富集分析.md", f"./result/{directory}/{directory_1}")
                            if content:
                                db_intro.append("KEGG")
                                f.write(content)
                                f.write(parse("./temp/5.2.1.KEGG注释统计柱状图.md", f"./result/{directory}/{directory_1}"))
                                f.write(parse("./temp/5.2.2.差异蛋白KEGG注释统计柱状图.md", f"./result/{directory}/{directory_1}"))
                                f.write(parse("./temp/5.2.3.差异蛋白KEGG富集柱状图.md", f"./result/{directory}/{directory_1}"))
                                f.write(parse("./temp/5.2.4.差异蛋白KEGG富集因子图.md", f"./result/{directory}/{directory_1}"))
                                f.write(parse("./temp/5.2.5.差异蛋白KEGG富集圈图.md", f"./result/{directory}/{directory_1}"))
                                f.write(parse("./temp/5.2.6.差异蛋白KEGG富集网络图.md", f"./result/{directory}/{directory_1}"))
                                f.write(parse("./temp/5.2.7.差异蛋白KEGG富集和弦图.md", f"./result/{directory}/{directory_1}"))
                                f.write(parse("./temp/5.2.8.差异蛋白KEGG通路图.md", f"./result/{directory}/{directory_1}"))
                            break

                    db_intro.append("UniProt")

                    # 5.3_domain
                    for directory_1 in os.listdir(f"./result/{directory}"):
                        if "domain" in directory_1:
                            content = parse("./temp/5.3.结构域分析.md", f"./result/{directory}/{directory_1}")
                            if content:
                                db_intro.append("Domain")
                                f.write(content)
                                f.write(parse("./temp/5.3.1.差异蛋白结构域注释统计柱状图.md", f"./result/{directory}/{directory_1}"))
                                f.write(parse("./temp/5.3.2.差异蛋白结构域富集柱状图.md", f"./result/{directory}/{directory_1}"))
                                f.write(parse("./temp/5.3.3.差异蛋白结构域富集因子图.md", f"./result/{directory}/{directory_1}"))
                                f.write(parse("./temp/5.3.4.差异蛋白结构域富集圈图.md", f"./result/{directory}/{directory_1}"))
                                f.write(parse("./temp/5.3.5.差异蛋白结构域富集网络图.md", f"./result/{directory}/{directory_1}"))
                                f.write(parse("./temp/5.3.6.差异蛋白结构域富集和弦图.md", f"./result/{directory}/{directory_1}"))
                            break

                    # 5.4_Reactome
                    for directory_1 in os.listdir(f"./result/{directory}"):
                        if "Reactome" in directory_1:
                            content = parse("./temp/5.4.Reactome分析.md", f"./result/{directory}/{directory_1}")
                            if content:
                                db_intro.append("Reactome")
                                f.write(content)
                                f.write(parse("./temp/5.4.1.差异蛋白Reactome注释统计柱状图.md", f"./result/{directory}/{directory_1}"))
                                f.write(parse("./temp/5.4.2.差异蛋白Reactome富集柱状图.md", f"./result/{directory}/{directory_1}"))
                                f.write(parse("./temp/5.4.3.差异蛋白Reactome富集因子图.md", f"./result/{directory}/{directory_1}"))
                                f.write(parse("./temp/5.4.4.差异蛋白Reactome富集圈图.md", f"./result/{directory}/{directory_1}"))
                                f.write(parse("./temp/5.4.5.差异蛋白Reactome富集网络图.md", f"./result/{directory}/{directory_1}"))
                                f.write(parse("./temp/5.4.6.差异蛋白Reactome富集和弦图.md", f"./result/{directory}/{directory_1}"))
                            break

                    # 5.5_DO
                    for directory_1 in os.listdir(f"./result/{directory}"):
                        if "DO" in directory_1:
                            content = parse("./temp/5.5.DO富集分析.md", f"./result/{directory}/{directory_1}")
                            if content:
                                db_intro.append("DO")
                                f.write(content)
                                f.write(parse("./temp/5.5.1.差异蛋白DO注释统计柱状图.md", f"./result/{directory}/{directory_1}"))
                                f.write(parse("./temp/5.5.2.差异蛋白DO富集柱状图.md", f"./result/{directory}/{directory_1}"))
                                f.write(parse("./temp/5.5.3.差异蛋白DO富集因子图.md", f"./result/{directory}/{directory_1}"))
                                f.write(parse("./temp/5.5.4.差异蛋白DO富集圈图.md", f"./result/{directory}/{directory_1}"))
                                f.write(parse("./temp/5.5.5.差异蛋白DO富集网络图.md", f"./result/{directory}/{directory_1}"))
                                f.write(parse("./temp/5.5.6.差异蛋白DO富集和弦图.md", f"./result/{directory}/{directory_1}"))
                            break

                    # 5.6_WikiPathway
                    for directory_1 in os.listdir(f"./result/{directory}"):
                        if "WikiPathway" in directory_1:
                            content = parse("./temp/5.6.Wikipathway分析.md", f"./result/{directory}/{directory_1}")
                            if content:
                                db_intro.append("Wikipathway")
                                f.write(content)
                                f.write(parse("./temp/5.6.1.差异蛋白Wikipathway注释统计柱状图.md", f"./result/{directory}/{directory_1}"))
                                f.write(parse("./temp/5.6.2.差异蛋白Wikipathway富集柱状图.md", f"./result/{directory}/{directory_1}"))
                                f.write(parse("./temp/5.6.3差异蛋白Wikipathway富集因子图.md", f"./result/{directory}/{directory_1}"))
                                f.write(parse("./temp/5.6.4.差异蛋白Wikipathway富集圈图.md", f"./result/{directory}/{directory_1}"))
                                f.write(parse("./temp/5.6.5.差异蛋白Wikipathway富集网络图.md", f"./result/{directory}/{directory_1}"))
                                f.write(parse("./temp/5.6.6.差异蛋白Wikipathway富集和弦图.md", f"./result/{directory}/{directory_1}"))
                            break
                    break

            # 6_GSEA
            for directory in os.listdir(f"./result"):
                if "GSEA" in directory:
                    f.write("## GSEA 分析\n")
                    for directory_1 in os.listdir(f"./result/{directory}"):
                        if "GO_GSEA" in directory_1:
                            f.write(parse("./temp/6.1.GSEA-GO分析.md", f"./result/{directory}/{directory_1}"))
                            break
                    for directory_1 in os.listdir(f"./result/{directory}"):
                        if "KEGG_GSEA" in directory_1:
                            f.write(parse("./temp/6.2.GSEA-KEGG.md", f"./result/{directory}/{directory_1}"))
                            break
                    for directory_1 in os.listdir(f"./result/{directory}"):
                        if "Reactome_GSEA" in directory_1:
                            f.write(parse("./temp/6.3.GSEA-Reactome.md", f"./result/{directory}/{directory_1}"))
                            break
                    for directory_1 in os.listdir(f"./result/{directory}"):
                        if "Domain_GSEA" in directory_1:
                            f.write(parse("./temp/6.4.GSEA-Domain.md", f"./result/{directory}/{directory_1}"))
                            break
                    break

            # 7_TF
            title_status = 0
            for directory in os.listdir(f"./result"):
                if "TF" in directory:
                    content = parse("./temp/7.1.转录因子家族分布.md", f"./result/{directory}")
                    if content:
                        if not title_status:
                            f.write("## 转录因子分析")
                            title_status = 1
                        f.write(content)
                    content = parse("./temp/7.2.差异表达转录因子柱状统计图.md", f"./result/{directory}")
                    if content:
                        if not title_status:
                            f.write("## 转录因子分析")
                            title_status = 1
                        f.write(content)
                    if title_status:
                        db_intro.append("TF")
                    break

            # 8_Subcellular-location
            for directory in os.listdir(f"./result"):
                if "Subcellular-location" in directory:
                    content = parse("./temp/8.亚细胞定位分析.md", f"./result/{directory}")
                    if content:
                        db_intro.append("Subcellular_Localization")
                        f.write(content)
                        f.write(parse("./temp/8.1.差异蛋白亚细胞定位分布饼图.md", f"./result/{directory}"))
                        f.write(parse("./temp/8.2.差异蛋白亚细胞定位柱状统计图.md", f"./result/{directory}"))
                    break

            # 9_PPI
            for directory in os.listdir(f"./result"):
                if "PPI" in directory:
                    content = parse("./temp/9.蛋白网络互作分析.md", f"./result/{directory}")
                    if content:
                        db_intro.append("PPI")
                        f.write(content)
                    break

            # 10_Kinase
            for directory in os.listdir(f"./result"):
                if "Kinase" in directory:
                    content = parse("./temp/10.激酶分析.md", f"./result/{directory}")
                    if content:
                        f.write(content)
                        db_intro.append("Kinase")
                    break

            # 11_Advanced_analyse
            if args.advanced_analyse:
                f.write(get_content("./temp/11.高级分析.md"))

            # 12_附录
            f.write("## 附录\n")
            f.write("### 数据库介绍\n")
            for db in db_intro:
                if os.path.exists(f"./temp/12.1.{db}.md"):
                    f.write(get_content(f"./temp/12.1.{db}.md"))

            # 13_材料方法
            f.write(get_content(f"./temp/13.材料方法_{args.type}.md"))

            # 14_参考文献
            f.write("## 参考文献\n")
            f.write("[1] Wu J, An Y, Pu H, Shan Y, Ren X, An M, Wang Q, Wei S, Ji J. Enrichment ofserum low-molecular-weight proteins using C18 absorbent under urea/dithiothreitoldenatured environment. Anal Biochem. 2010 Mar 1;398(1):34-44.\n")
            f.write("[2] Wyant GA, Abu-Remaileh M, Frenkel EM, et al. NUFIP1 is a ribosome receptor for starvation-induced ribophagy. Science. 2018;360(6390):751-758.\n")
            i = 3
            for db in db_intro:
                for file in os.listdir(f"./temp/"):
                    if f"14.{db}" in file:
                        f.write(f"[{i}] " + get_content(f"./temp/{file}"))
                        i += 1



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Create md")
    parser.add_argument("--advanced-analyse", help="enable advanced analyse", action="store_true")
    parser.add_argument("-t", "--type", choices=["DIA", "LFQ", "TMT", "PQA"], help="project type: [DIA, LFQ, TMT, PQA]", type=str)
    args = parser.parse_args()
    main(args.type)