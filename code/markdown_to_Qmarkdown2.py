import os,sys
import re
import optparse
from  collections import defaultdict
import time
import itertools
"""
超链接:[Treat/3.1.1_phylo_profile/](../01.fastqc_dir/samples.stat.xls){target="_blank"} 
图片格式：![]()  #这种可以直接放pdf
"""



def Time():
    return time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())

def make_dir(dir):
    if not os.path.exists(dir):
        dir=os.path.expanduser(dir)
        os.makedirs(dir)

header = f"""---
#title: "test"
#author: "lhy"
execute:
  echo: false
  warning: false
  message: false
editor: visual
page-layout: full
title: DIA蛋白结果分析报告
format:
  html:
    embed-resources: true
    smooth-scroll: true
    code-link: true
    theme:
      light: Cerulean
      dark: Cyborg
    toc: true
    toc-depth: 6
    toc-location: left # left-靠左，right-靠右
    toc-title: <img class="logo" src="./static/img/logo.png"/><span>DIA蛋白结果分析报告</span>
    html-math-method: mathjax
    include-in-header:
      - text: |
          <link href="./static/css/index.css" type="text/css" rel="stylesheet" />
          <link href="./static/css/raxus.css" type="text/css" rel="stylesheet" />
          <link href="./static/css/carousel.css" type="text/css" rel="stylesheet" />
          <link href="./static/css/custom.css" type="text/css" rel="stylesheet" />
          <script src="./static/js/jquery-1.11.3.min.js" type="text/javascript"></script>
          <script src="./static/js/raxus-slider.min.js" type="text/javascript"></script>
          <script src="./static/js/vertical-carousel.js" type="text/javascript"></script>
          <script src="./static/js/common.js" type="text/javascript"></script>
    include-before-body: ./static/html/cover.html
    include-after-body: ./static/html/footer.html
---

"""
header += """
```{css, echo=FALSE}  
h6,.h6,h5,.h5,h4,.h4,h3,.h3,h2,.h2,h1,.h1 {
    margin-top: 0;
    margin-bottom: .5rem;
    font-weight: 500;
    line-height: 1.2;
    color: #2fa4e7;
}
```
"""


def table_convert(table,n):
    #解析&[](./table/project_info.txt)
    table_name=re.search(r"&\[(.*?)\]", table).group(1)
    table_path=re.search(r"\]\((.*?)\)", table).group(1)
    if table_name in ["亚细胞定位富集表", "结构域分析表", "亚细胞定位表", "GO 富集分析表", "KEGG 富集分析表", "结构域富集表",
                      "Reactome 富集表", "DO 富集表", "Wikipathway 分析表", "亚细胞定位富集表"]:
        align = 'l'
    else:
        align = 'c'
    if len(table_name) == 0:
            table_info=f"""
```{{r}}
#| label: tbl-table{n}
"""
    else:
        table_info=f"""
```{{r}}
#| label: tbl-table{n}
#| tbl-cap: "{table_name}"
"""
    table_info+=f"""
library(knitr)
#install.packages("kableExtra")
library(kableExtra)
library(migrate)
library(data.table)
if (tools::file_ext("{table_path}") == "xlsx"){{
library(readxl)
general <- as.data.table(read_xlsx("{table_path}", n_max = 500))
}}else{{general <- fread("{table_path}",header = T,sep="\\t",encoding = "UTF-8",nrows=500)}}

if (nrow(general) > 20) {{
            kable(general, align = "{align}",caption = "{table_name}") %>% kable_styling() %>%
                scroll_box(width = "100%", height = "400px") 
                }} else{{
            kable(general, align = "{align}")}}
```
"""
    return table_info


# def imagebox_convert(imagebox,all_info,n):
#     #解析imagebox[ QC 样本总离子流图重叠谱图](./image/1_QC/1.1_TIC/{model}/TIC.png)
#     imagebox_name=re.search(r'imagebox\[(.*?)\]', imagebox).group(1)
#     imagebox_path=re.search(r'\]\((.*?)\)', imagebox).group(1)
#     variable_info=set(re.findall("\{(.*?)\}",imagebox,re.S))
#     all_variable=[all_info[i] for i in variable_info if i in all_info]
#     imagebox_info=f"""<!--  多图模块  -->

# ::: {{.raxus-slider data-thumbnail="bottom" data-keypress="true" data-autoplay="3000" data-arrows="show"}}
# ::: slider-relative"""
#     all_imagebox_path=[]
#     imagebox_variable_path=imagebox_path
#     for i in variable_info:
#         if i in all_info:
#             for j in all_info[i]:
#                 pattern = f'{{{i}}}'
#                 imagebox_variable_path=re.sub(pattern, f"{j}", imagebox_variable_path)
#                 all_imagebox_path.append(imagebox_variable_path)
#                 print(imagebox_variable_path)
#     print(all_imagebox_path)
#     all_imagebox_path_uniq=set([ i for i in all_imagebox_path if not re.search("\{.*?\}",i )])
#     print(all_imagebox_path_uniq)
#     for i in  all_imagebox_path_uniq:
#         imagebox_info+=f"""
# ::: slide
# [![]({imagebox_variable_path})]({imagebox_variable_path})
# :::

# """
#     imagebox_info+=f"""
# :::
# :::

# ::: figcaption
# Fig {n}: {imagebox_name}
# :::
# """
#     return imagebox_info


def imagebox_convert(imagebox,all_info,n):
    #解析imagebox[ QC 样本总离子流图重叠谱图](./image/1_QC/1.1_TIC/{model}/TIC.png)
    text = re.findall("\((.*)\)", imagebox, re.S)[0]
    imagebox_name = re.findall("\[(.*)\]", imagebox, re.S)[0]
    variable = {x.split(":")[0].strip(): x.split(":")[1].strip() for x in text.split(",")}
    pngs = []
    imagebox_info = f"""<!--  多图模块  -->

::: {{.raxus-slider data-thumbnail="bottom" data-keypress="true" data-autoplay="3000" data-arrows="show"}}
::: {{.slider-relative}}"""
    for folder in os.listdir(variable["path"]):
        if os.path.isdir(os.path.join(variable["path"], folder)):
            if "depth" in variable.keys():
                for folder1 in os.listdir(os.path.join(variable["path"], folder)):
                    if os.path.isdir(os.path.join(variable["path"], folder, folder1)):
                        for file in os.listdir(os.path.join(variable["path"], folder, folder1)):
                            if file.endswith(variable["format"]):
                                pngs.append((folder, os.path.join(variable["path"], folder, folder1, file).replace("\\", "/")))
                                break
            else:
                for file in os.listdir(os.path.join(variable["path"], folder)):
                    if variable["keyword"] in file and file.endswith(variable["format"]):
                        pngs.append((folder, os.path.join(variable["path"], folder, file).replace("\\", "/")))
                        break
    for png_name, png_path in pngs:
        imagebox_info += f"""
::: {{.slide}}
[![{png_name}]({png_path})]({png_path})
:::

"""
    imagebox_info += f"""
:::
:::

::: {{.figcaption}}
Fig {n}: {imagebox_name}
:::
"""
    return imagebox_info

#     imagebox_name=re.search(r'imagebox\[(.*?)\]', imagebox).group(1)
#     imagebox_path=re.search(r'\]\((.*?)\)', imagebox).group(1)
#     variable_info=set(re.findall("\{(.*?)\}",imagebox,re.S))
#     all_variable=[all_info[i] for i in variable_info if i in all_info]
#     imagebox_info=f"""<!--  多图模块  -->
#
# ::: {{.raxus-slider data-thumbnail="bottom" data-keypress="true" data-autoplay="3000" data-arrows="show"}}
# ::: {{.slider-relative}}"""
#     imagebox_variable_path=imagebox_path
#     keys = variable_info
#     all_element=[]
#     for i in variable_info:
#         all_element.append(all_info[i])
#     # values = [list(all_info[i].values())[0] for i in all_variable_fixed ]
#     combinations = itertools.product(*all_element)
#     new_filepath=imagebox_path
#     # new_result_annotation=result_annotation
#     for combination in combinations:
#         png_name=""
#         filepath=new_filepath
#         # result_annotation=new_result_annotation
#         for i, key in enumerate(keys):
#             filepath = filepath.replace("{" + key + "}", combination[i])
#             png_name+=combination[i]+"."
#         # for i, key in enumerate(keys):
#         #     result_annotation = result_annotation.replace("{" + key + "}", combination[i])
#         pattern = r'{(.*?)}'
#         while re.search(pattern, filepath):
#             filepath = re.sub(pattern, r'\1', filepath)
#         result_path=os.path.dirname(filepath)
#         png_name=png_name[:-1]
#         print(png_name)
#         # while re.search(pattern, result_annotation):
#         #     result_annotation = re.sub(pattern, r'\1', result_annotation)
#         imagebox_info+=f"""
# ::: {{.slide}}
# [![{png_name}]({filepath})]({filepath})
# :::
#
# """
#     imagebox_info+=f"""
# :::
# :::
#
# ::: {{.figcaption}}
# Fig {n}: {imagebox_name}
# :::
# """
#     return imagebox_info
    

def imagebox_convert2(imagebox,all_info,n):
    #解析imagebox[ QC 样本总离子流图重叠谱图](./image/1_QC/1.1_TIC/{model}/TIC.png)
    imagebox_name=re.search(r'imagebox\[(.*?)\]', imagebox).group(1)
    imagebox_path=re.search(r'\]\((.*?)\)', imagebox).group(1)
    variable_info=set(re.findall("\{(.*?)\}",imagebox,re.S))
    all_variable=[all_info[i] for i in variable_info if i in all_info]
    imagebox_info_name=""
    imagebox_info_path=""
    imagebox_info_head=f"""<!--  多图模块  -->

::: {{.vertical-carousel}}
::: {{.vertical-carousel-left}}"""
    imagebox_variable_path=imagebox_path
    keys = variable_info
    print(keys)
    all_element=[]
    for i in variable_info:
        all_element.append(all_info[i])
    # values = [list(all_info[i].values())[0] for i in all_variable_fixed ]
    print(f"Ceshi{all_element}")
    print(variable_info)
    combinations = itertools.product(*all_element)
    new_filepath=imagebox_path
    # new_result_annotation=result_annotation
    for combination in combinations:
        png_name=""
        filepath=new_filepath
        # result_annotation=new_result_annotation
        for i, key in enumerate(keys):
            filepath = filepath.replace("{" + key + "}", combination[i])
            png_name+=combination[i]+"."
        # for i, key in enumerate(keys):
        #     result_annotation = result_annotation.replace("{" + key + "}", combination[i])
        pattern = r'{(.*?)}'
        while re.search(pattern, filepath):
            filepath = re.sub(pattern, r'\1', filepath)
        result_path=os.path.dirname(filepath)
        png_name=png_name[:-1]
        # while re.search(pattern, result_annotation):
        #     result_annotation = re.sub(pattern, r'\1', result_annotation)
        imagebox_info_name+=f"""
::: {{.vertical-carousel-label}}
{png_name}
:::"""
        imagebox_info_path+=f"""
::: {{.vertical-carousel-img}}
[![]({filepath}){{.carousel-img-percent-90}}]({filepath}){{target="_blank"}}
:::"""
        print(f'{imagebox_info_name}')
        print(f'{imagebox_info_path}')
    imagebox_info=f"""{imagebox_info_head}{imagebox_info_name}\n:::\n::: {{.vertical-carousel-right}}{imagebox_info_path}\n:::\n:::\n"""
    imagebox_info+=f"""
::: {{.figcaption}}
Fig {n}: {imagebox_name}
:::
"""
    return imagebox_info
    
    
    
#     for i in variable_info:
#         if i in all_info:
#             for j in all_info[i]:
#                 pattern = f'{{{i}}}'
#                 imagebox_variable_path=re.sub(pattern, f"{j}", imagebox_variable_path)
#                 all_imagebox_path.append(imagebox_variable_path)
#                 print(imagebox_variable_path)
#     print(all_imagebox_path)
#     all_imagebox_path_uniq=set([ i for i in all_imagebox_path if not re.search("\{.*?\}",i )])
#     print(all_imagebox_path_uniq)
#     for i in  all_imagebox_path_uniq:
#         imagebox_info+=f"""
# ::: slide
# [![]({imagebox_variable_path})]({imagebox_variable_path})
# :::

# """
#     imagebox_info+=f"""
# :::
# :::

# ::: figcaption
# Fig {n}: {imagebox_name}
# :::
# """
#     return imagebox_info



def image_convert(image,all_info,n):
    #解析![Fig 9: 基于MEGAN的多样本分类等级树图](./template/aa.pdf){width="100%" height="600"}
    image_name=re.search(r'!\[(.*?)\]', image).group(1)
    image_path=re.search(r'\]\((.*?)\)', image).group(1)
    if re.search(r'\]\(.*?\)\{(.*?)\}',image):
        image_add=re.search(r'\]\(.*?\)(.*?)',image).group(1)
    else:
        image_add=""
    image_info=f"![Fig{n} {image_name} ]({image_path}){image_add}\n\n"
    return image_info

def result_file(line,all_info):
    txt=""
    result_annotation=re.findall('\[(.*?)\]',line,re.S)[0]
    filepath = re.findall('\((.*?)\)',line,re.S)[0]
    result_path=os.path.dirname(filepath)
    
    if re.findall("\{.*?\}",line,re.S):
        new_result_annotation=result_annotation
        new_filepath=filepath
        variable=set(re.findall("\{(.*?)\}",line,re.S))
        keys = variable
        all_element=[]
        for i in variable:
            all_element.append(all_info[i])
        # values = [list(all_info[i].values())[0] for i in all_variable_fixed ]
        combinations = itertools.product(*all_element)
        # new_result_annotation=result_annotation
        for combination in combinations:
            result_annotation=new_result_annotation
            filepath=new_filepath
            for i, key in enumerate(keys):
                filepath = filepath.replace("{" + key + "}", combination[i])
            for i, key in enumerate(keys):
                result_annotation = result_annotation.replace("{" + key + "}", combination[i])
            pattern = r'{(.*?)}'
            while re.search(pattern, filepath):
                filepath = re.sub(pattern, r'\1', filepath)
            result_path=os.path.dirname(filepath)
            while re.search(pattern, result_annotation):
                result_annotation = re.sub(pattern, r'\1', result_annotation)
            # print(f"""{result_annotation}：[{filepath}]({filepath}){{.unindent target="_blank"}}""")
            txt+=f"""{result_annotation}：[{filepath}]({filepath}){{.unindent target="_blank"}}  \n"""

        # all_result_file=[]
        # variable=set(re.findall("\{(.*?)\}",line,re.S))
        # for i in variable:
        #     if i in all_info:
        #         for j in all_info[i]:
        #             pattern = f'{{{i}}}'
        #             result_annotation=re.sub(pattern, f"{j}", result_annotation)
        #             filepath=re.sub(pattern, f"{j}", filepath)
        #             result_path=re.sub(pattern, f"{j}", result_path)
        #             # print(f"""{result_annotation1}：[{result_path1}]({filepath1}){{target="_blank"}}""")
        #             all_result_file.append(f"""{result_annotation}：[{result_path}]({filepath}){{target="_blank"}}\n\n""")
        # all_result_file_uniq=set([ i for i in all_result_file if not re.search("\{.*?\}",i )])
        # for i in all_result_file_uniq:
        #     txt+=i
    else:
        txt+=f"""{result_annotation}：[{filepath}]({filepath}){{.unindent target="_blank"}}  \n"""
    return txt

def note(line):
    if line.startswith("[table_note]"):
        table_anntion_text=re.findall("\[table_note\](.*)",line,re.S)[0]
        text=f"""
::: {{.explain .custom-p-number-3}}
{table_anntion_text}
:::
"""
    elif line.startswith("[image_note]"):
        image_anntion_text=re.findall("\[image_note\](.*)",line,re.S)[0]
        text=f"""
::: {{.explain .custom-p-number-3}}
{image_anntion_text}
:::
"""
    else:
        text=""
    text+="\n\n"
    return text


def get_variable(line,all_info):
    keys = set(re.findall("\{(.*?)\}", line, re.S))
    for key in keys:
        line = line.replace("{" + key + "}", all_info[key][0] if all_info[key] else "")
    return line

    
def parse(md_f,all_info):
    twotitle_number=0
    threetitle_number=0
    fourtitle_number=0
    image_num=1
    table_num=1
    text = header
    result_file_button=0
    text_header_change=0
    with open(md_f, encoding='utf-8') as md:
        for line in md:
            line = line.strip()
            if line.startswith("#"):
                if  re.match(r'^(#+)\s+(.*)', line):
                    match = re.match(r'^(#+)\s+(.*)', line)
                    level_text = match.group(1)
                    level_info = match.group(2)
                    level = len(match.group(1))
                # print(str(level)+"\t"+line)
                    if level == 1:
                        if text_header_change==0:
                            text = header.replace("DIA蛋白结果分析报告", f"{level_info}")
                            text_header_change=1
                        continue
                    if level == 2:
                        twotitle_number += 1
                        threetitle_number=0
                        fourtitle_number=0
                    elif level == 3 :
                        threetitle_number += 1
                        fourtitle_number=0
                    elif level == 4:
                        fourtitle_number += 1
                    chapter_info=".".join([str(element) for element in [twotitle_number,threetitle_number,fourtitle_number]if element != 0])+" "+match.group(2)
                    text+=f"#{level_text} {chapter_info}\n\n"
                else:
                    continue
            elif line.startswith("!"):
                if set(re.findall("\{(.*?)\}", line, re.S)):
                    line = get_variable(line, all_info)
                text+=f"{line}\n\n"
                image_convert(line,all_info,image_num)
                image_num+=1
            elif line.startswith("&"):
                text+=table_convert(line,table_num)
                table_num+=1
            elif line.startswith("[table_note]") or line.startswith("[image_note]") :
                text+=note(line)
            elif line.startswith("[link]"):
                link_content=re.search(r"\[link\](.*)",line).group(1)+'{target="_blank"}\n\n'
                text+=link_content
            elif line.startswith("imagebox["):
                text+=imagebox_convert(line,all_info,image_num)
                image_num+=1
            elif line.startswith("[span]"):
                text+="<p></p>"

            elif line.startswith("[result_file_start]"):
                result_file_button=1
                text+=f"""::: {{.result-files}}
::: {{.result-button .custom-p-number-5}}
结果文件
:::
"""
            elif line.startswith("[result_file_end]"):
                result_file_button=0
                text+=""":::

"""
            else:
                if  result_file_button ==1 :
                    text+=result_file(line,all_info)
                    # if re.findall("\{.*?\}",line,re.S):
                    #     all_variable=set(re.findall("\{(.*?)\}",line,re.S))
                    #     text+=result_file(line,all_info,variable=all_variable)
                    # else:
                    #     text+=result_file(line,all_info,variable=[])
                elif set(re.findall("\{(.*?)\}", line, re.S)):
                    text+=get_variable(line, all_info)
                else:
                    text+=f"{line}\n\n"
    return text

    
if __name__ == "__main__":
    parser=optparse.OptionParser(usage='"usage:%prog [options] arg1,arg2"',version="%prog 1.2")  
    parser.add_option('-m','--markdown',  
                    action='store',dest='markdown',  
                    help='markdown file')
    parser.add_option('-p','--prefix',  
                    action='store',dest='prefix',default = 'report',  
                    help='report file prefix')
    parser.add_option('-o','--output',  
                    action='store',dest='output',default = './',  
                    help='report directory [default:%default]')
    parser.add_option('-i','--info',  
                    action='store',dest='info', 
                    help='All variable information files [default:%default]')
    
    options,args=parser.parse_args()
    if not options.markdown:
        os.system("python3 "+sys.argv[0]+" -h")
        sys.exit(1)
    outdir=os.path.realpath(options.output)
    prefix=options.prefix
    make_dir(outdir)
    
    all_varialble_info=defaultdict(dict)
    if os.path.isfile(options.info):
        with  open(options.info,"r",encoding="utf-8") as fi:
            for line in fi.readlines():
                lines=line.strip().split("\t")
                variable_name = lines[0]
                variable_info = lines[1:]
                all_varialble_info[variable_name]=variable_info
    with open(f'{outdir}/{prefix}.qmd',"w",encoding="utf-8") as fw:
        fw.write(parse(options.markdown,all_varialble_info))
    