from anndata import AnnData
from typing import Union, Optional, Literal
from numpy.random import RandomState
import scanpy as sc
from typing import Union, Optional, Any, Mapping, Callable
from types import MappingProxyType
import numpy as np
import warnings
from sklearn.decomposition import PCA

_Method = Literal["umap", "gauss", "rapids"]
_MetricFn = Callable[[np.ndarray, np.ndarray], float]
# from sklearn.metrics.pairwise_distances.__doc__:
_MetricSparseCapable = Literal[
    "cityblock", "cosine", "euclidean", "l1", "l2", "manhattan"
]
_MetricScipySpatial = Literal[
    "braycurtis",
    "canberra",
    "chebyshev",
    "correlation",
    "dice",
    "hamming",
    "jaccard",
    "kulsinski",
    "mahalanobis",
    "minkowski",
    "rogerstanimoto",
    "russellrao",
    "seuclidean",
    "sokalmichener",
    "sokalsneath",
    "sqeuclidean",
    "yule",
]
_Metric = Union[_MetricSparseCapable, _MetricScipySpatial]


def find_elbow(adata, reduction="pca", method="slope", ndims=30):
    """
    在降维结果中寻找肘点（最佳维度）

    参数:
        adata: AnnData对象（包含降维结果）
        reduction: 降维方法名称（默认为'pca'）
        method: 计算方法，'slope'（斜率变化）或'kerstin'（内在维度估计）
        ndims: 考虑的最大维度数（默认30）

    返回:
        optimal_dims: 推荐的最佳维度数
    """
    method = method.lower()
    if method not in ["slope", "kerstin"]:
        raise ValueError("method must be 'slope' or 'kerstin'")

    # 获取标准差数据
    if reduction not in adata.uns:
        raise ValueError(f"No information for reduction '{reduction}' found in .uns")

    variance_key = f"{reduction}_variance"
    if variance_key not in adata.uns:
        # 尝试从PCA结果中获取方差
        if "pca" in adata.uns and "variance" in adata.uns["pca"]:
            variance = adata.uns["pca"]["variance"]
        else:
            # 计算嵌入矩阵的标准差作为备选
            embedding_key = f"X_{reduction}"
            if embedding_key not in adata.obsm:
                raise ValueError(
                    f"Neither variance nor embeddings found for {reduction}"
                )
            warnings.warn(
                "Using standard deviation from embeddings (variance not found)"
            )
            embeddings = adata.obsm[embedding_key]
            stdev = np.std(embeddings, axis=0)
            variance = stdev**2
    else:
        variance = adata.uns[variance_key]

    stdev = np.sqrt(variance)

    # 检查ndims是否超出范围
    if ndims > len(stdev):
        warnings.warn(f"Object only has information for {len(stdev)} reductions")
        ndims = len(stdev)

    if method == "slope":
        return _slope_method(stdev, ndims)
    else:  # kerstin
        return _kerstin_method(adata, reduction, ndims)


def _slope_method(stdev, ndims):
    """斜率变化法寻找肘点"""
    if ndims <= 1:
        return 1

    # 取前ndims个标准差
    stdev = stdev[:ndims]
    n = len(stdev)

    # 计算相邻维度的斜率（应为负值）
    slopes = []
    for i in range(1, n):
        slope = stdev[i] - stdev[i - 1]
        slopes.append(slope)

    if not slopes:  # 只有一个维度
        return 1

    max_slope = max(slopes)  # 最小负值（最接近零的负数）
    min_stdev = min(stdev)
    max_stdev = max(stdev)

    # 遍历斜率寻找肘点
    for j in range(len(slopes)):
        current_slope = slopes[j]
        current_stdev = stdev[j + 1]  # 当前维度对应的标准差

        # 条件1：斜率变化显著（大于10倍历史最大斜率）
        cond1 = current_slope > (10 * max_slope)

        # 条件2：当前标准差接近最小值（小于最小值+5%的极差）
        cond2 = current_stdev < (min_stdev + 0.05 * (max_stdev - min_stdev))

        if cond1 and cond2:
            return j + 1  # 返回维度数（索引j对应第j+1个维度）

    return 1  # 未找到肘点时返回默认值


def _kerstin_method(adata, reduction, ndims):
    """内在维度估计法"""
    from skdim.id import MLE

    embedding_key = f"X_{reduction}"
    if embedding_key not in adata.obsm:
        raise ValueError(f"Embeddings not found for {reduction} in .obsm")

    data = adata.obsm[embedding_key][:, :ndims]

    try:
        # 使用全局最大似然估计内在维度
        mle = MLE(k1=20, k2=20)  # 使用20个最近邻
        dim_est = mle.fit(data).dimension_
        return round(dim_est)
    except ImportError:
        raise RuntimeError(
            "Please install scikit-dimension: pip install scikit-dimension"
        )


def neighbors(
    adata: AnnData,
    n_neighbors: int = 15,
    n_pcs: Optional[int] = None,
    use_rep: Optional[str] = None,
    knn: bool = True,
    random_state: Optional[Union[int, RandomState]] = 0,
    method: Optional[_Method] = "umap",
    metric: Union[_Metric, _MetricFn] = "euclidean",
    metric_kwds: Mapping[str, Any] = MappingProxyType({}),
    copy: bool = False,
    **kwargs,
) -> Optional[AnnData]:
    """\
    Compute a neighborhood graph of observations [McInnes18]_.
    The neighbor search efficiency of this heavily relies on UMAP [McInnes18]_,
    which also provides a method for estimating connectivities of data points -
    the connectivity of the manifold (`method=='umap'`). If `method=='gauss'`,
    connectivities are computed according to [Coifman05]_, in the adaption of
    [Haghverdi16]_.
    Parameters
    ----------
    adata
        Annotated data matrix.
    n_neighbors
        The size of local neighborhood (in terms of number of neighboring data
        points) used for manifold approximation. Larger values result in more
        global views of the manifold, while smaller values result in more local
        data being preserved. In general values should be in the range 2 to 100.
        If `knn` is `True`, number of nearest neighbors to be searched. If `knn`
        is `False`, a Gaussian kernel width is set to the distance of the
        `n_neighbors` neighbor.
    n_pcs
        Use this many PCs. If n_pcs==0 use .X if use_rep is None.
    use_rep
        Use the indicated representation. 'X' or any key for .obsm is valid.
        If None, the representation is chosen automatically: For .n_vars < 50, .X is used, otherwise ‘X_pca’ is used.
        If ‘X_pca’ is not present, it’s computed with default parameters.
    knn
        If `True`, use a hard threshold to restrict the number of neighbors to
        `n_neighbors`, that is, consider a knn graph. Otherwise, use a Gaussian
        Kernel to assign low weights to neighbors more distant than the
        `n_neighbors` nearest neighbor.
    random_state
        A numpy random seed.
    method
        Use ‘umap’ [McInnes18] or ‘gauss’ (Gauss kernel following [Coifman05] with adaptive width
        [Haghverdi16]) for computing connectivities. Use ‘rapids’ for the RAPIDS implementation of UMAP
        (experimental, GPU only).
    copy
        Return a copy instead of writing to adata.
    Returns
    -------
    Depending on `copy`, updates or returns `adata` with the following:
    **connectivities** : sparse matrix (`.uns['neighbors']`, dtype `float32`)
        Weighted adjacency matrix of the neighborhood graph of data
        points. Weights should be interpreted as connectivities.
    **distances** : sparse matrix (`.uns['neighbors']`, dtype `float32`)
        Instead of decaying weights, this stores distances for each pair of
        neighbors.
    """

    sc.pp.neighbors(
        adata,
        n_neighbors=n_neighbors,
        n_pcs=n_pcs,
        use_rep=use_rep,
        knn=knn,
        method=method,
        random_state=random_state,
        metric=metric,
        metric_kwds=metric_kwds,
        copy=copy,
        **kwargs,
    )

    print("Created k-Nearest-Neighbor graph in adata.uns['neighbors'] ")


def clustering(
    adata: AnnData,
    method: Union["louvain", "leiden"] = "leiden",
    resolution: float = 1.2,
    **kwargs,
):
    """
    Perform community detection on single-cell data using Louvain or Leiden algorithm.

    Parameters:
    - adata: AnnData
        An AnnData object containing single-cell data.
    - method: str, optional (default: 'louvain')
        The community detection algorithm to use. Options: 'louvain' or 'leiden'.
    - resolution: float, optional
        Resolution parameter for the community detection algorithm. If None, the default algorithm value is used.
    - louvain_key: str, optional (default: 'louvain')
        The key to store Louvain communities in adata.obs.
    - leiden_key: str, optional (default: 'leiden')
        The key to store Leiden communities in adata.obs.
    - **kwargs:
        Additional keyword arguments to be passed to the underlying community detection function.

    Returns:
    - None (in-place modification of adata)
    """
    if method == "louvain":
        sc.tl.louvain(adata, key_added=method, resolution=resolution, **kwargs)
    elif method == "leiden":
        sc.tl.leiden(adata, key_added=method, resolution=resolution, **kwargs)
    else:
        raise ValueError(f"Unsupported {method} method. Choose 'louvain' or 'leiden'.")
