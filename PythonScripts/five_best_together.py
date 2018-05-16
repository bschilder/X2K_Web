import pandas as pd

df = pd.read_csv('process.csv').dropna()
ppis = [col for col in df.columns if col.startswith("enable_ppi.")]
df_results = df.groupby(ppis)['score'].agg({
    'total': 'count',
    'weighted_total': 'sum',
})
df_results['weighted_score'] = df_results['weighted_total'] / df_results['total']
df_results = df_results.sort_values('weighted_score', ascending=False)


def arr_to_names(names, vals):
    return [
        name
        for name, val in zip(names, vals)
        if val == 1
    ]

def print_top_n(df_results, n=5):
    for ind, (_, row) in enumerate(df_results.iterrows()):
        if ind > n:
            break
        print(
            arr_to_names(df_results.index.names, row.name),
            row.values
        )
