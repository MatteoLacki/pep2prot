import pandas as pd


def sum_top(X, top_no, group, top=True):
    cols = (X[[c]].sort_values(by=[group,c], ascending=top).groupby(group).head(top_no)
            for c in X.columns)
    return pd.concat(cols, axis=1).groupby(group).sum()
