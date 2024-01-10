def plot_tsne(features, valid, second_features=None, first_tag='NovelFrq',
              second_tag='TNF', labels=None):
    """
    Plot T-distributed Stochastic Neighbor Embedding for data.

    Parameters
    ----------
    data : numpy.ndarry, pandas.DataFrame of shape (n_dna, n_freq)
        features of data.
    valid : pandas.DataSeries or pandas.DataFrame of shape(n_dna, 1)
        valid cluster for data.
    Returns
    -------
    out : Nothing.
    """
    import numpy as np
    import pandas as pd
    from sklearn.manifold import TSNE

    if type(valid) is pd.core.frame.DataFrame:
        valid = valid[valid.columns[0]]
    if type(valid) is str:
        valid = pd.read_csv(valid, index_col=0, header=None)[1]

    Tsne = pd.DataFrame(columns=[first_tag + '_Tsne-one',
                                 first_tag + '_Tsne-two'],
                        index=features.index)
    Tsne[[first_tag + '_Tsne-one', first_tag + '_Tsne-two']] =\
        TSNE(n_components=2, init='random').fit_transform(features)

    if second_features is not None:
        se_tsne = pd.DataFrame(columns=[second_tag + '_Tsne-one',
                                        second_tag + '_Tsne-two'],
                               index=second_features.index)
        se_tsne[[second_tag + '_Tsne-one', second_tag + '_Tsne-two']] =\
            TSNE(n_components=2, init='random').fit_transform(second_features)
        Tsne[[second_tag + '_Tsne-one', second_tag + '_Tsne-two']] =\
            se_tsne[[second_tag + '_Tsne-one', second_tag + '_Tsne-two']]

    Tsne['label'] = valid

    if labels is None:
        labels = valid.unique()
    if type(labels) is int:
        from random import shuffle
        ll = valid.unique()
        shuffle(ll)
        labels = ll[:labels]

    if second_features is None:
        import plotly.express as px
        fig = px.scatter(Tsne[Tsne['label'].isin(labels)],
                         x=first_tag+"_Tsne-one",
                         y=first_tag+"_Tsne-two", color="label")
        fig.show()
    else:
        import matplotlib.pyplot as plt
        import seaborn as sns
        plt.figure(figsize=(16, 6))
        sns.scatterplot(
            x=first_tag+"_Tsne-one", y=first_tag+"_Tsne-two", hue="label",
            palette=sns.color_palette("hls", len(labels)),
            data=Tsne[Tsne['label'].isin(labels)],
            legend="full", alpha=0.3,
            ax=plt.subplot(1, 2, 1)
        )
        sns.scatterplot(
            x=second_tag+"_Tsne-one", y=second_tag+"_Tsne-two", hue="label",
            palette=sns.color_palette("hls", len(labels)),
            data=Tsne[Tsne['label'].isin(labels)],
            legend="full", alpha=0.3,
            ax=plt.subplot(1, 2, 2)
        )
        plt.show()

    return Tsne


def plot_dist(features, valid, kind='euq'):
    """
    """
    import numpy as np
    import pandas as pd
    if type(valid) is pd.core.frame.DataFrame:
        valid = valid[valid.columns[0]]
    if type(valid) is str:
        valid = pd.read_csv(valid, index_col=0, header=None)[1]

    def euq_dist(a, b):
        return np.abs(a - b[:, None]).sum(2) / np.sum(a[0])

    def square_dist(a, b):
        return np.power(a - b[:, None], 2).sum(2)

    if kind == 'square':
        dist_func = square_dist
    else:
        dist_func = euq_dist

    inter = []
    for i in range(len(valid.unique())):
        a = features.loc[valid == valid.unique()[i]].to_numpy()
        b = features.loc[valid.isin(valid.unique()[i+1:])].to_numpy()
        dists = dist_func(a, b)
        inter += list(dists.reshape(-1).astype('float32'))
    intra = []
    for i in range(len(valid.unique())):
        a = features.loc[valid == valid.unique()[i]].to_numpy()
        dists = dist_func(a, a)
        intra += [j for j in dists.reshape(-1).astype('float32') if j != 0]

    import matplotlib.pyplot as plt
    import seaborn as sns
    sns.distplot(inter, norm_hist=True, hist=False, kde=True, label='inter')
    sns.distplot(intra, norm_hist=True, hist=False, kde=True, label='intra')
    plt.legend(prop={'size': 16}, title='label')
    plt.show()
    return [np.mean(inter), np.std(intra)], [np.mean(intra), np.std(intra)]
