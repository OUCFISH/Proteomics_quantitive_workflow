#!/usr/bin/bash
# 功能分析部分

# 定义颜色和格式代码
GREEN='\033[0;32m'
RED='\033[0;31m'
BOLD='\033[1m'
NC='\033[0m' # No Color

# 定义变量记录脚本运行状态
error_occurred=false
failed_analyses=()

# 显示帮助信息的函数
show_help() {
    echo "Usage: $0 --advanced-analysis"
    echo ""
    echo "This script runs advanced analyses including:"
    echo "  - GO Enrichment Analysis"
    echo "  - KEGG Enrichment Analysis"
    echo "  - Other enrichment analyses"
    echo "  - GSEA analyses"
    echo "  - Transcription Factor Analysis"
    echo "  - Subcellular Location Analysis"
    echo "  - PPI Analysis"
    echo "  - Kinase Analysis"
    echo "  - Heatmap and Boxplot Analysis"
    echo ""
    echo "Options:"
    echo "  --advanced-analysis   Run all the advanced analysis scripts."
    echo "  --no-advanced-analysis   Skip all the advanced analysis scripts."
    echo "  -h, --help            Display this help message."
    echo ""
    echo "Example:"
    echo "  $0 --advanced-analysis"
    echo "  $0 --no-advanced-analysis"
}

if [ "$1" == "help" ] || [ "$1" == "-h" ] || [ "$1" == "--help" ] || 
   ([ "$1" != "--advanced-analysis" ] && [ "$1" != "--no-advanced-analysis" ]); then
    show_help
    exit 0
fi

# 定义一个函数来运行R脚本并输出提示信息
run_rscript() {
  local script=$1
  local config=$2
  local success_msg=$3
  local failure_msg=$4

  if [ -z "$config" ]; then
    if Rscript "$script" 2>/dev/null; then
      echo -e "${GREEN}${success_msg}${NC}"
    else
      echo -e "${RED}${failure_msg}${NC}"
      error_occurred=true
      failed_analyses+=("$failure_msg")
    fi
  else
    if Rscript "$script" "$config" 2>/dev/null; then
      echo -e "${GREEN}${success_msg}${NC}"
    else
      echo -e "${RED}${failure_msg}${NC}"
      error_occurred=true
      failed_analyses+=("$failure_msg")
    fi
  fi
}

# 定义一个函数来运行 Python 脚本并输出提示信息
run_pyscript() {
  local script=$1
  local arg=$2
  local success_msg=$3
  local failure_msg=$4

  if /home/admin/miniconda3/envs/prot/bin/python "$script" "$arg"; then
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
# 功能注释 - GO 富集分析
run_rscript "./code/5.1_GO_enrichment.R" "" "GO Enrichment Analysis was completed" "GO Enrichment Analysis failed"

# 功能注释 - KEGG 富集分析
run_rscript "./code/5.2_KEGG_enrichment.R" "" "KEGG Enrichment Analysis was completed" "KEGG Enrichment Analysis failed"
run_rscript "./code/5.2_KEGG_html.R" "" "KEGG pathway Analysis was completed" "KEGG pathway Analysis failed"

# 功能注释 - other 富集分析
run_rscript "./code/5.3_other_enrichment.R" "" "other enrichment Analysis was completed" "other enrichment Analysis failed"
convert_pdfs_to_pngs_recursive "./result/5_Enrichment/5.3_domain/"

# 功能注释 - plot
run_rscript "./code/5.plot.R" "" "Enrichment plot Analysis was completed" "Enrichment plot Analysis failed"

# 功能注释 - convert
convert_pdfs_to_pngs_recursive "./result/5_Enrichment/5.1_GO/"
convert_pdfs_to_pngs_recursive "./result/5_Enrichment/5.2_KEGG/"
convert_pdfs_to_pngs_recursive "./result/5_Enrichment/5.4_Reactome/"
convert_pdfs_to_pngs_recursive "./result/5_Enrichment/5.5_DO/"
convert_pdfs_to_pngs_recursive "./result/5_Enrichment/5.6_WikiPathway/"

# GSEA 分析 - GO
run_rscript "./code/6.1_GOGSEA_analysis.R" "" "GOGSEA Analysis was completed" "GOGSEA Analysis failed"
run_rscript "./code/6.2_GOGSEAPlot_analysis.R" "" "GOGSEA Plot Analysis was completed" "GOGSEA Plot Analysis failed"
convert_pdfs_to_pngs_recursive "./result/6_GSEA/6.1_GO_GSEA/"
# GSEA 分析 - KEGG
run_rscript "./code/6.3_KEGGGSEA_analysis.R" "" "KEGGGSEA Analysis was completed" "KEGGGSEA Analysis failed"
run_rscript "./code/6.4_KEGGGSEAPlot_analysis.R" "" "KEGGGSEA Plot Analysis was completed" "KEGGGSEA Plot Analysis failed"
convert_pdfs_to_pngs_recursive "./result/6_GSEA/6.2_KEGG_GSEA/"
# GSEA 分析 - Reactome
run_rscript "./code/6.5_ReactomeGSEA_analysis.R" "" "ReactomeGSEA Analysis was completed" "ReactomeGSEA Analysis failed"
run_rscript "./code/6.6_ReactomeGSEAPlot_analysis.R" "" "ReactomeGSEA Plot Analysis was completed" "ReactomeGSEA Plot Analysis failed"
convert_pdfs_to_pngs_recursive "./result/6_GSEA/6.3_Reactome_GSEA/"
# GSEA 分析 - Interpro
run_rscript "./code/6.7_InterproGSEA_analysis.R" "" "InterproGSEA Analysis was completed" "InterproGSEA Analysis failed"
run_rscript "./code/6.8_InterproGSEAPlot_analysis.R" "" "InterproGSEA Plot Analysis was completed" "InterproGSEA Plot Analysis failed"
convert_pdfs_to_pngs_recursive "./result/6_GSEA/6.4_Domain_GSEA/"
# 转录因子分析
run_rscript "./code/7_TF.R" "" "TF Analysis was completed" "TF Analysis failed"
convert_pdfs_to_pngs_recursive "./result/7_TF/"
# 亚细胞定位分析
run_rscript "./code/8_Subcellular-location-analysis.R" "" "Subcellular Location Analysis was completed" "Subcellular Location Analysis failed"
convert_pdfs_to_pngs_recursive "./result/8_Subcellular-location/"
# PPI 分析
run_rscript "./code/9_PPI.R" "" "PPI Analysis was completed" "PPI Analysis failed"
run_rscript "./code/9_PPI_plot.R" "" "PPI plot Analysis was completed" "PPI plot Analysis failed"
run_pyscript "./code/ppi_plot.py" "./result/9_PPI/" "PPI Plot Analysis was completed" "PPI Plot Analysis failed"

# 激酶分析
run_rscript "./code/10_Kinase.R" "" "Kinase Analysis was completed" "Kinase Analysis failed"
convert_pdfs_to_pngs_recursive "./result/10_Kinase/"
# Protein 匹配
run_pyscript "./code/match_protein.py" "" "Protein Matching was completed" "Protein Matching failed"

# 检查命令行参数
if [ "$1" == "--advanced-analysis" ]; then
  echo "Running advanced analyses..."

  # 高级分析 - 热图和箱线图
  run_rscript "./code/11.1_heatmap_boxplot_analysis.R" "./config/11.1_heatmap_boxplot_analysis.conf" "Heatmap Boxplot Analysis was completed" "Heatmap Boxplot Analysis failed"

  # 高级分析 - 热图和趋势图
  run_rscript "./code/11.2_heatmap_trendplot_analysis.R" "./config/11.2_heatmap_trendplot_analysis.conf" "Heatmap Trendplot Analysis was completed" "Heatmap Trendplot Analysis failed"

  # 高级分析 - 热图和富集分析
  run_rscript "./code/11.3_heatmap_Enrichment_analysis.R" "./config/11.3_heatmap_Enrichment_analysis.conf" "Heatmap Enrichment Analysis was completed" "Heatmap Enrichment Analysis failed"

  # 高级分析 - GSEA复合图
  run_rscript "./code/11.4_GOGSEA-multiPathwayPlot_analysis.R" "" "GOGSEA MultiPathway Plot Analysis was completed" "GOGSEA MultiPathway Plot Analysis failed"
  run_rscript "./code/11.5_KEGGGSEA-multiPathwayPlot-analysis.R" "" "KEGGGSEA MultiPathway Plot Analysis was completed" "KEGGGSEA MultiPathway Plot Analysis failed"
  convert_pdfs_to_pngs_recursive "./result/11_AdcancedPlot/11.4_GSEA-multiPathwayPlot/"

elif [ "$1" == "--no-advanced-analysis" ]; then
  echo "Skipping advanced analyses."

else
  echo "Invalid argument: $1"
  show_help
  exit 1
fi

# 如果全部分析均运行成功，输出醒目信息
if [ "$error_occurred" = false ]; then
  echo -e "${GREEN}${BOLD}=========================================${NC}"
  echo -e "${GREEN}${BOLD}   FUNCTIONAL MODULE ANALYSIS IS COMPLETE ${NC}"
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
