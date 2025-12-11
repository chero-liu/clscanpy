'''
作者: xiufeng.yang xiufeng.yang@oebiotech.com
日期: 2025-02-19 13:49:01
最后编辑者: xiufeng.yang xiufeng.yang@oebiotech.com
最后编辑时间: 2025-02-27 10:26:55
文件路径: \oe-smt\run_celltypist.py
描述: Celltypist模型训练脚本
'''
import click
import scanpy as sc
import celltypist
import time
import numpy as np
import logging
import sys,os
from utils import setup_logger, convert_seurat_to_anndata
 
logger = setup_logger()
def show_config(config: dict):
    logger.info(f"⚙️ Configuration:")
    for key, value in config.items():
        logger.info(f"\t🛠️ {key}: {value}")
def load_and_preprocess_data(input_path):
    """加载并预处理输入数据
    
    参数:
        input_path: h5ad或者h5seurat文件路径
        
    返回:
        adata_ref: 预处理后的AnnData对象
    """
    logger.info(f"加载数据: {input_path}")
    if input_path.endswith('.h5seurat') or input_path.endswith('.rds') :
        adata_ref = convert_seurat_to_anndata(input_path,use_raw_counts=True)
    else:
        adata_ref = sc.read_h5ad(input_path)
 
    logger.info("检查数据是否已经进行了log转换")
    is_valid_range = (adata_ref.X[:1000].min() >= 0) and (adata_ref.X[:1000].max() <= 9.22)
    is_valid_sum = np.abs(np.expm1(adata_ref.X[0]).sum() - 1e4) <= 1
    
    logger.info(f"数据范围检查: {is_valid_range}")
    logger.info(f"标准化检查: {is_valid_sum}")
    
    if not (is_valid_range and is_valid_sum):   
        logger.info("数据未进行正确的log1p标准化，进行normalize_total和log1p处理...")
        if adata_ref.raw is not None:
            var_names = adata_ref.var_names
            adata_ref = sc.AnnData(
                X=adata_ref.raw.X, 
                obs=adata_ref.obs,
                var=adata_ref.raw.var,
                uns=adata_ref.uns,
                obsm=adata_ref.obsm
            )
            try:
                adata_ref.var_names = var_names
            except:
                pass 
            sc.pp.normalize_total(adata_ref, target_sum=1e4)
            sc.pp.log1p(adata_ref)
        else:
            logger.info("未找到raw数据,请核查数据是否正确")
            exit(0)
    else:
        logger.info("数据已经进行了正确的log1p标准化")

    return adata_ref

def train_model(adata_ref, 
                cell_type_key, 
                output_path, 
                prefix,
                downsample_cells=None, 
                use_sgd=False,
                feature_selection=False,
                n_features=300, 
                n_jobs=10, 
                max_iter=100):
    """训练celltypist模型
    
    参数:
        adata_ref: AnnData对象，包含训练数据
        cell_type_key: 细胞类型标签在obs中的列名
        output_path: 模型保存路径
        prefix: 模型保存文件名前缀
        downsample_cells: 降采样后的细胞数量，None表示不降采样
        use_sgd: 是否使用SGD优化
        feature_selection: 是否进行特征选择
        n_features: 特征选择保留的基因数量
        n_jobs: 并行作业数
        max_iter: 最大迭代次数
    """
    logger.info(f"检查adata_ref中各种{cell_type_key}细胞数目情况")
    
    # 添加以下代码，输出各细胞类型的数目情况
    if cell_type_key in adata_ref.obs.columns:
        cell_counts = adata_ref.obs[cell_type_key].value_counts()
        logger.info(f"细胞类型数量统计:")
        for cell_type, count in cell_counts.items():
            logger.info(f"  {cell_type}: {count}个细胞")
        logger.info(f"总计: {len(cell_counts)}种细胞类型, {adata_ref.n_obs}个细胞")
    else:
        logger.warning(f"警告: 在adata_ref.obs中未找到{cell_type_key}列")

    if downsample_cells:
        sampled_cell_index = celltypist.samples.downsample_adata(
            adata_ref, 
            mode='each', 
            n_cells=downsample_cells,
            by=cell_type_key, 
            return_index=True, 
            balance_cell_type=True
        )
        logger.info(f"降采样后用于训练的细胞数量: {len(sampled_cell_index)}")
        adata_ref = adata_ref[sampled_cell_index]
        
        # 添加以下代码，输出降采样后各细胞类型的数目情况
        if cell_type_key in adata_ref.obs.columns:
            cell_counts_after = adata_ref.obs[cell_type_key].value_counts()
            logger.info(f"降采样后细胞类型数量统计:")
            for cell_type, count in cell_counts_after.items():
                logger.info(f"  {cell_type}: {count}个细胞")
            logger.info(f"降采样后总计: {len(cell_counts_after)}种细胞类型, {adata_ref.n_obs}个细胞")
    #config settings
    config = { "X":adata_ref,
               "downsample_cells": downsample_cells,
               "labels": cell_type_key,
                "check_expression": True,
                ##LR param
                "C": 1.0,
                "solver": None, 
                "n_jobs": n_jobs, 
                "max_iter":max_iter, 
                ## SGD param
                "use_SGD": use_sgd, 
                "alpha": 0.0001,
                 ## GPU param
                "use_GPU": False,
                ## mini-batch 
                "mini_batch": True,
                "batch_number": 100,
                "batch_size": 1000,
                "epochs": 10,
                "balance_cell_type":True,
                ## feature selection
                "feature_selection":feature_selection, 
                "top_genes": n_features
    }
 
    #quiet or not
    show_config(config)

    t_start = time.time()
    model = celltypist.train(adata_ref, 
                            cell_type_key,
                            check_expression=False,
                            ##LR param
                            C=1.0,
                            solver=None, 
                            n_jobs=n_jobs, 
                            max_iter=max_iter, 
                            ## SGD param
                            use_SGD=use_sgd, 
                            alpha=0.0001,
                            ## GPU param
                            use_GPU=False,
                            ## mini-batch 
                            mini_batch=True,
                            batch_number=100,
                            batch_size=1000,
                            epochs=10,
                            balance_cell_type=True,
                            ## feature selection
                            feature_selection=feature_selection, 
                            top_genes=n_features)
 
    t_end = time.time()
    logger.info(f"模型训练耗时: {(t_end - t_start)/60:.2f} 分钟")

    # 保存模型
    os.makedirs(output_path,exist_ok=True)
    model.write(f"{output_path}/{prefix}.pkl")
    logger.info(f"模型已保存至: {output_path}/{prefix}.pkl")

@click.command()
@click.option('--input-path', required=True, 
              help='输入文件路径，支持h5ad,h5seurat格式，包含原始表达矩阵和细胞类型标注')
@click.option('--output-path', required=True, default='./model', 
              help='模型保存路径，将保存为pkl格式')
@click.option('--prefix', required=True, default='celltypist_model', 
              help='模型保存文件名前缀')
@click.option('--cell-type-key', 
              help='细胞类型标签在adata.obs中的列名， 需要提供')
@click.option('--downsample-cells', type=int, 
              help='每种细胞类型中降采样细胞数量，不设置则使用全部细胞')
@click.option('--use-sgd', type=bool, default=False, 
              help='训练方式，是否使用SGD (随机梯度下降) 逻辑回归，默认为False，使用标准逻辑回归')
@click.option('--feature-selection', type=bool, default=True,
              help='是否执行两阶段数据训练，其中第一阶段通过 SGD（随机梯度下降）学习筛选重要特征/基因。若设为 True，训练时间将延长')
@click.option('--n-features', default=300, type=int, 
              help='对于每一种细胞类型，在进行特征选择时候保留的基因数量，默认为300')
@click.option('--n-jobs', default=10, type=int, 
              help='并行计算的作业数，默认为10')
@click.option('--max-iter',  type=int, 
              help='最大迭代次数，不指定的话会自动设置， 200, 500, and 1000 for large (>500k cells), medium (50-500k), and small (<50k) datasets')
def main(input_path, output_path, prefix, cell_type_key, downsample_cells,
         use_sgd, feature_selection, n_features, n_jobs, max_iter):
    """Celltypist模型训练脚本

    此脚本用于训练Celltypist细胞类型注释模型。
     

    主要功能包括:
 
    1. 数据预处理：标准化和log转换
    2. 可选的细胞降采样
    3. 可选的特征选择
    4. 模型训练和保存
    """
    logger.info(f"step1:加载并预处理数据=======================================")
    adata_ref = load_and_preprocess_data(input_path)
    logger.info(f"加载数据形状: {adata_ref.shape}")
    
    logger.info(f"step2:训练并保存模型=========================================")
    train_model(
        adata_ref=adata_ref,
        cell_type_key=cell_type_key,
        output_path=output_path,
        prefix=prefix,
        downsample_cells=downsample_cells,
        use_sgd=use_sgd,
        feature_selection=feature_selection,
        n_features=n_features,
        n_jobs=n_jobs,
        max_iter=max_iter
    )

if __name__ == '__main__':
    main()