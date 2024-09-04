#!/usr/bin/bash
# Quality control

# 定义颜色和格式代码
GREEN='\033[0;32m'
RED='\033[0;31m'
BOLD='\033[1m'
NC='\033[0m' # No Color

# 定义变量记录脚本运行状态
error_occurred=false
failed_analyses=()

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


# 质控部分
run_rscript "./code/3.1.1_peptides_length_distribution.R" "Peptides Length Distribution Analysis was completed" "Peptides Length Distribution Analysis failed"
convert_pdfs_to_pngs_recursive "./result/3_Quality_Control/3.1_Quality_Control/3.1.1_peptides_length_distribution/"

run_rscript "./code/3.1.2_protein_molecular_distribution.R" "Protein Molecular Distribution Analysis was completed" "Protein Molecular Distribution Analysis failed"
convert_pdfs_to_pngs_recursive "./result/3_Quality_Control/3.1_Quality_Control/3.1.2_protein_molecular_distribution/"

run_rscript "./code/3.1.3_proteins_coverage_distribution.R" "Proteins Coverage Distribution Analysis was completed" "Proteins Coverage Distribution Analysis failed"
convert_pdfs_to_pngs_recursive "./result/3_Quality_Control/3.1_Quality_Control/3.1.3_proteins_coverage_distribution/"

run_rscript "./code/3.1.4_proteins_isoelectric_distribution.R" "Proteins Isoelectric Distribution Analysis was completed" "Proteins Isoelectric Distribution Analysis failed"
convert_pdfs_to_pngs_recursive "./result/3_Quality_Control/3.1_Quality_Control/3.1.4_proteins_isoelectric_distribution/"

run_rscript "./code/3.1.5_IRT.R" "IRT Analysis was completed" "IRT Analysis was failed"
convert_pdfs_to_pngs_recursive  " ./result/3_Quality_Control/3.1_Quality_Control/3.1.5_IRT/"

run_rscript "./code/3.2.1_proteins_number_barplot.R" "Proteins Number Barplot Analysis was completed" "Proteins Number Barplot Analysis failed"
convert_pdfs_to_pngs_recursive "./result/3_Quality_Control/3.2_Quantitative_Statistics/3.2.1_Expression_Statistics/"

run_rscript "./code/3.2.2_venn.R" "Venn Analysis was completed" "Venn Analysis failed"
run_rscript "./code/3.2.2_upset.R" "UpSet Analysis was completed" "UpSet Analysis failed"
convert_pdfs_to_pngs_recursive "./result/3_Quality_Control/3.2_Quantitative_Statistics/3.2.2_Expression_Repeat/"

run_rscript "./code/3.2.3_proteins_intensity_distribution.R" "Proteins Intensity Distribution Analysis was completed" "Proteins Intensity Distribution Analysis failed"
convert_pdfs_to_pngs_recursive " ./result/3_Quality_Control/3.2_Quantitative_Statistics/3.2.3_Expression_Violin/"

run_rscript "./code/3.2.4_quality-control-RSD.R" "Quality Control RSD Analysis was completed" "Quality Control RSD Analysis failed"
convert_pdfs_to_pngs_recursive "./result/3_Quality_Control/3.2_Quantitative_Statistics/3.2.4_RSD/"

run_rscript "./code/3.2.5_proteins_intensity_scatterplot.R" "Proteins Intensity Scatterplot Analysis was completed" "Proteins Intensity Scatterplot Analysis failed"
convert_pdfs_to_pngs_recursive "./result/3_Quality_Control/3.2_Quantitative_Statistics/3.2.5_Expression_Scatter/"

run_rscript "./code/3.2.6_process_correlation.R" "Process Correlation Analysis was completed" "Process Correlation Analysis failed"
run_rscript "./code/3.2.6_correlation_analysis.R" "Correlation Analysis was completed" "Correlation Analysis failed"
convert_pdfs_to_pngs_recursive "./result/3_Quality_Control/3.2_Quantitative_Statistics/3.2.6_Correlation_Analysis/"
#run_rscript "./code/3.2.7_Global_Heatmap.R" "Global Heatmap Analysis was completed" "Global Heatmap Analysis failed"
#run_rscript "./code/3.2.8_quality-control-PCA.R" "Quality Control PCA Analysis was completed" "Quality Control PCA Analysis failed"

# 差异表达部分
#run_rscript "./code/4.1.diffexp.R" "Diff Expression Analysis was completed" "Diff Expression Analysis failed"
run_rscript "./code/4.2_venn.R" "Differential Repeat Venn Analysis was completed" "Differential Repeat Venn Analysis failed"
run_rscript "./code/4.2_upset.R" "Differential Repeat UpSet Analysis was completed" "Differential Repeat UpSet Analysis failed"
convert_pdfs_to_pngs_recursive "./result/4_Diff_Expressed/4.2_Diff_Repeat/"

run_rscript "./code/4.3_Volcano.R" "Volcano Plot Analysis was completed" "Volcano Plot Analysis failed"
convert_pdfs_to_pngs_recursive "./result/4_Diff_Expressed/4.3_Volcano/"

run_rscript "./code/4.4_Heatmap.R" "Heatmap Analysis was completed" "Heatmap Analysis failed"
convert_pdfs_to_pngs_recursive "./result/4_Diff_Expressed/4.4_Heatmap/"

run_rscript "./code/4.5_Trend.R" "Trend Analysis was completed" "Trend Analysis failed"
convert_pdfs_to_pngs_recursive "./result/4_Diff_Expressed/4.5_Trend/"

run_rscript "./code/4.6_radar_chart.R" "Radar Analysis was completed" "Radar Analysis failed"
convert_pdfs_to_pngs_recursive "./result/4_Diff_Expressed/4.6_Radar/"

# 如果全部分析均运行成功，输出醒目信息
if [ "$error_occurred" = false ]; then
  echo -e "${GREEN}${BOLD}=========================================${NC}"
  echo -e "${GREEN}${BOLD}   QUALITY CONTROL ANALYSIS IS COMPLETE ${NC}"
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
