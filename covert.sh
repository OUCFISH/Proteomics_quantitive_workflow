#!/bin/bash

# 定义颜色和格式代码
GREEN='\033[0;32m'
RED='\033[0;31m'
BOLD='\033[1m'
NC='\033[0m' # No Color

# 定义变量记录脚本运行状态
error_occurred=false
failed_analyses=()

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
      failed_analyses+=("$pdf_file") # 将失败的文件添加到列表
    else
      echo -e "${GREEN}Successfully converted $pdf_file to PNG.${NC}"
    fi
  done

  # 如果有错误发生，打印失败的文件列表
  if $error_occurred; then
    echo -e "${RED}The following files failed to convert:${NC}"
    for failed_file in "${failed_analyses[@]}"; do
      echo "$failed_file"
    done
  fi
}

# 脚本主入口点
# 检查是否提供了文件夹路径作为参数
if [ "$#" -eq 1 ]; then
  convert_pdfs_to_pngs_recursive "$1"
else
  echo -e "${BOLD}Usage: ${0} <folder_path>${NC}"
  exit 1
fi
