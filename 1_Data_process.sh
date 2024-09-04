#!/usr/bin/bash
# Data process

# 定义颜色和格式代码
GREEN='\033[0;32m'
RED='\033[0;31m'
BOLD='\033[1m'
NC='\033[0m' # No Color

# 定义变量记录脚本运行状态
error_occurred=false
failed_analyses=()

# 创建1_Sampleinfo文件夹和12_Attachments
sampleinfo_dir="./result/1_Sampleinfo"
Attachments_dir="./result/12_Attachments"
data_dir="./data"

if [ ! -d "$sampleinfo_dir" ]; then
  mkdir -p "$sampleinfo_dir"
fi

if [ ! -d "$Attachments_dir" ]; then
  mkdir -p "$Attachments_dir"
fi

for file in "$data_dir"/*分组*.xlsx; 
do
  if [ -e "$file" ]; then
    cp "$file" "$sampleinfo_dir"
  fi
done

for file in "$data_dir"/*.fasta; 
do
  if [ -e "$file" ]; then
    cp "$file" "$Attachments_dir"
  fi
done

# 定义一个函数来运行命令并输出提示信息
run_rscript() {
  local script=$1
  local success_msg=$2
  local failure_msg=$3

  if Rscript "$script" 2>/dev/null; then
    echo -e "${GREEN}${success_msg}${NC}"
  else
    echo -e "${RED}${failure_msg}${NC}"
    error_occurred=true
    failed_analyses+=("$failure_msg")
  fi
}

# 定义一个函数来递归将PDF文件转换为PNG文件
convert_pdfs_to_pngs_recursive() {
  local folder_path=$1

  # 使用find命令递归查找所有PDF文件
  find "$folder_path" -type f -name "*.pdf" | while read -r pdf_file; do
    # 获取文件所在文件夹路径
    file_dir=$(dirname "$pdf_file")

    # 获取文件名，不包含扩展名
    base_name=$(basename "$pdf_file" .pdf)

    # 转换PDF文件为PNG文件，设置像素密度为300 DPI
    if ! convert -flatten -density 300 "$pdf_file" "$file_dir/$base_name.png" 2>/dev/null; then
      echo -e "${RED}Failed to convert $pdf_file to PNG.${NC}"
      error_occurred=true
    fi
  done
}

# 显示帮助信息的函数
show_help() {
    echo "Usage: $0 <analysis_type>"
    echo ""
    echo "Parameters:"
    echo "  <analysis_type>  The type of analysis to perform. Choose from:"
    echo "                   'DIA' - Data-Independent Acquisition"
    echo "                   'LFQ' - Label-Free Quantification"
    echo "                   'TMT' - Tandem Mass Tag"
    echo ""
    echo "Example:"
    echo "  $0 -DIA"
    echo "  $0 -LFQ"
    echo "  $0 -TMT"
}
# 检查是否传入了帮助命令或参数
if [ "$1" == "help" ] || [ "$1" == "-h" ] || [ "$1" == "--help" ]; then
    show_help
    exit 0
fi

# 检查是否传入了分析类型参数
if [ -z "$1" ]; then
    echo "Error: No analysis type provided."
    echo "Use 'help' or '-h' or '--help' for usage information."
    exit 1
fi

# 获取分析类型参数
analysis_type="$1"

# 根据分组表格生成group和contrast表
run_rscript "./code/1_sampleinfo.R" "Sample info generation was completed" "Sample info generation failed"

# 根据分析类型选择并运行相应的数据预处理命令
if [ "$analysis_type" == "-DIA" ]; then
    run_rscript "./code/2_data-process_DIA.R" "Data processing was completed" "Data processing failed"
elif [ "$analysis_type" == "-LFQ" ]; then
    run_rscript "./code/2_data-process_LFQ.R" "Data processing was completed" "Data processing failed"
elif [ "$analysis_type" == "-TMT" ]; then
    run_rscript "./code/2_data-process_TMT.R" "Data processing was completed" "Data processing failed"
else
    echo "Unknown analysis type: $analysis_type"
    exit 1
fi

# 3.2.7_Global_Heatmap
run_rscript "./code/3.2.7_Global_Heatmap.R" "Global Heatmap Analysis was completed" "Global Heatmap Analysis failed"
convert_pdfs_to_pngs_recursive "./result/3_Quality_Control/3.2_Quantitative_Statistics/3.2.7_Global_Heatmap/"

# 3.2.8_PCA
run_rscript "./code/3.2.8_quality-control-PCA.R" "PCA Analysis was completed" "PCA Analysis failed"
convert_pdfs_to_pngs_recursive "./result/3_Quality_Control/3.2_Quantitative_Statistics/3.2.8_PCA/"

# 4_Diff_Expressed
# 4.1_DiffStats
run_rscript "./code/4.1.diffexp.R" "Differential Expression Analysis was completed" "Differential Expression Analysis failed"
convert_pdfs_to_pngs_recursive "./result/4_Diff_Expressed/4.1_DiffStats/"

# 如果全部分析均运行成功，输出醒目信息
if [ "$error_occurred" = false ]; then
  echo -e "${GREEN}${BOLD}=========================================${NC}"
  echo -e "${GREEN}${BOLD}   DATA PROCESS ANALYSIS IS COMPLETE ${NC}"
  echo -e "${GREEN}${BOLD}             IT'S PRETTY COOL            ${NC}"
  echo -e "${GREEN}${BOLD}=========================================${NC}"
else
  echo -e "${RED}${BOLD}=========================================${NC}"
  echo -e "${RED}${BOLD}      SOME ANALYSES FAILED TO RUN        ${NC}"
  echo -e "${RED}${BOLD}=========================================${NC}"
  echo -e "${RED}${BOLD} Failed analyses: ${NC}"
  for failed in "${failed_analyses[@]}"; do
    echo -e "${RED}${BOLD}  - $failed ${NC}"
  done
  echo -e "${RED}${BOLD}=========================================${NC}"
fi
