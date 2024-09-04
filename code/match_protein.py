import os

import pandas as pd


def match_protein(raw_file, match_file):
    raw_df = pd.read_csv(raw_file)
    match_df = pd.read_csv(match_file)
    result = raw_df.merge(match_df, how='inner', left_on='Symbol', right_on='Genes')
    result = result[['Protein','Symbol']]
    result = result.drop_duplicates()
    result_df = raw_df.merge(result, how='left', on='Symbol')
    result_df = result_df.dropna(subset=['Protein'])
    result_df.to_csv('./input/pathway2protein.tsv',sep='\t', index=False)


def format_annotation(annotation, form):
    if form == 'xlsx':
        df = pd.read_excel(annotation)
    else:
        df = pd.read_csv(annotation)
    df.to_csv(annotation.replace(f'.{form}', '.tsv'), sep='\t', index=False)


if __name__ == '__main__':
    # match_protein('./input/mmu_KO_pathway_description.csv', './input/EAPTG_vs_SPTG-difp.csv')
    # format_annotation('./input/Annotation.xlsx')
    for file in os.listdir('./temp/'):
        if file.endswith('.xlsx'):
            format_annotation('./temp/'+file,'xlsx')
        if file.endswith('.csv'):
            format_annotation('./temp/'+file,'csv')

    for file in os.listdir('./result/4_Diff_Expressed/4.1_DiffStats/'):
        if file.endswith('.xlsx'):
            format_annotation('./result/4_Diff_Expressed/4.1_DiffStats/'+file,'xlsx')
        if file.endswith('.csv'):
            format_annotation('./result/4_Diff_Expressed/4.1_DiffStats/'+file,'csv')