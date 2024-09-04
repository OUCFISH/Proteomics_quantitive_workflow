#### 在result/1_Sampleinfo/中须有信息单的pdf文件
#! /usr/bin/bash

#第一个参数为项目类型
project_type=$1

#第二个参数为高级分析
option=$2

case $project_type in
    DIA|TMT|LFQ|PQA)
        echo "Project type: $project_type"
        ;;
    *)
        echo "Invalid project type. Valid types are: DIA, TMT, LFQ, PQA"
        exit 1
        ;;
esac

if [ "$option" = "--advanced-analyse" ]; then
    echo "Advanced analytics enabled."
else
    echo "Advanced analytics not enabled"
fi

cp -r ./config/static ./
cp ./code/markdown_to_Qmarkdown2.py ./markdown_to_Qmarkdown2.py
cp ./code/create_md.py ./create_md.py
#### 默认无高级分析，需要高级分析添加--advanced-analyse
if [ "$option" = "--advanced-analyse" ]; then
    /home/admin/miniconda3/envs/prot_report/bin/python ./create_md.py -t $1 $2
else
    /home/admin/miniconda3/envs/prot_report/bin/python ./create_md.py -t $1
fi
/home/admin/miniconda3/envs/prot_report/bin/python ./markdown_to_Qmarkdown2.py -m ./prot_report_$project_type.md -p report_$project_type -o ./ -i ./variable.txt
/PERSONALBIO1/software/quarto/quarto-1.5.45/bin/quarto render ./report_$project_type.qmd --to html
