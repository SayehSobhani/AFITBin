def guess_clusters_num(sequences, init_K=4, end_situation=0.1):
    """
    """
    from sklearn.cluster import KMeans
    from sklearn.metrics import silhouette_score as score
    
    last_score = abs(score(sequences, KMeans(init_K).fit_predict(sequences)) -
                     end_situation)
    while init_K < len(sequences) // 250:
        init_K += 1
        new_score = abs(score(sequences,
                              KMeans(init_K).fit_predict(sequences)) -
                        end_situation)
        if new_score > last_score:
            break
        last_score = new_score
    return init_K

