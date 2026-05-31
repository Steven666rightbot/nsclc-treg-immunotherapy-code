import pandas as pd

df = pd.read_csv('data/gse243013_geo_metadata.csv.gz')

print('=== MAJOR CELL TYPES ===')
print(df['major_cell_type'].value_counts())
print()

print('=== SUB CELL TYPES (all) ===')
sub_counts = df['sub_cell_type'].value_counts().sort_index()
for sub, count in sub_counts.items():
    major = df[df['sub_cell_type']==sub]['major_cell_type'].iloc[0]
    print(f'  {sub:35s} | {major:15s} | n={count:>8,}')

print(f'\nTotal cells: {len(df):,}')
print(f'Unique samples: {df["sampleID"].nunique()}')

# Check if all 51 subtypes from cell_proportions.csv are present
props = pd.read_csv('data/cell_proportions.csv')
print(f'\nSubtypes in cell_proportions.csv: {len(props.columns)}')
print(f'Subtypes in GEO metadata: {df["sub_cell_type"].nunique()}')

geo_subs = set(df['sub_cell_type'].unique())
prop_subs = set(props.columns)
print(f'\nMissing from GEO: {prop_subs - geo_subs}')
print(f'Extra in GEO: {geo_subs - prop_subs}')
