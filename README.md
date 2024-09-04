#此流程为定量蛋白质组学分析流程，兼容DIA、TMT、LFQ的分析项目

#此代码仓库不含有结果文件和示例数据，只包含分析所用的code、config、shell脚本和temp文件夹下面的md文件

#如需测试请自行将测试数据放到指定data路径中

#不同分析项目的数据表格不同

##DIA
蛋白表：report.pg_matrix.tsv
肽段表：report.pr_matrix.tsv

##LFQ
蛋白表：Proteins.xlsx
肽段表：Peptide_Groups.xlsx
样本对应关系表：Input Files.xlsx (表头是File ID和File Name)

##TMT
蛋白表：Proteins.xlsx
肽段表：Peptide_Groups.xlsx
TMT标记样本对应表：标记对照表.xlsx (报告内样本名称	标记)

#运行项目首先请新建data文件夹存放原始数据表和分组表格文件 并传入对应搜库物种的fasta文件

#随后创建result空文件夹

#运行项目请注意修改default_config文件（主要修改数据表的路径、物种信息、卡值信息）和11.3的config文件（GO和KEGG背景文件）

#DIA项目原始矩阵表头格式应该为类似 /PERSONALBIO1/prote/project/SP240402174_0725/C1_DIA.mzML，代码中取最后一个/和_之间的元素

#如果带IRT则需要传入IRT.xlsx 和 report.tsv到data文件夹下 前者为标准品保留时间

##分析步骤
bash 1_Data_process.sh -DIA -TMT -LFQ (根据不同分析类型选择相应的参数)
bash 2_Quality_control.sh 
bash 3_Functional_Annotation.sh --advanced-analysis --no-advanced-analysis (分为高级分析和不含高级分析)
bash 4_Report.sh DIA, TMT, LFQ, PQA (根据不同分析类型选择相应的参数)
