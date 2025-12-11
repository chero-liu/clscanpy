'''
作者: xiufeng.yang xiufeng.yang@oebiotech.com
日期: 2025-02-19 13:49:01
最后编辑者: xiufeng.yang xiufeng.yang@oebiotech.com
最后编辑时间: 2025-02-27 10:26:55
文件路径: \oe-smt\run_celltypist_annotate.py
描述: Celltypist细胞类型注释脚本
'''
import click
import os
import scanpy as sc
import celltypist
import time
import numpy as np
import pandas as pd
from pathlib import Path
from oescanpy.tools.read.utils import convert_seurat_to_anndata
from oescanpy.tools.utils import setup_logger
from utils import update_h5seurat_metadata
from color import select_colors   
import matplotlib.pyplot as plt
logger = setup_logger()

def show_config(config: dict):
    logger.info(f"⚙️ Configuration:")
    for key, value in config.items():
        logger.info(f"\t🛠️ {key}: {value}")
def show_help_and_exit(message: str):
    ctx = click.get_current_context()
    click.echo(click.style(message, fg="red"))
    click.echo()
    ctx.fail(ctx.get_help())
def load_and_preprocess_data(input_path,filter_gene=None):
    """加载并预处理输入数据
    
    参数:
        input_path: h5ad或者h5seurat文件路径
        
    返回:
        adata: 预处理后的AnnData对象
    """
    logger.info(f"加载数据: {input_path}")
    if input_path.endswith('.h5seurat') or input_path.endswith('.rds'):
        adata = convert_seurat_to_anndata(input_path,use_raw_counts=False)
    else:
        adata = sc.read_h5ad(input_path)
    # if filter_gene:
    #     logger.info("过滤表达细胞数少于10的基因...")
    #     adata= sc.pp.filter_genes(adata, min_cells=10)
    logger.info("检查数据是否已经进行了log转换")
    is_valid_range = (adata.X[:1000].min() >= 0) and (adata.X[:1000].max() <= 9.22)
    is_valid_sum = np.abs(np.expm1(adata.X[0]).sum() - 1e4) <= 1
    
    logger.info(f"数据范围检查: {is_valid_range}")
    logger.info(f"标准化检查: {is_valid_sum}")
    
    if not (is_valid_range and is_valid_sum):   
        logger.info("数据未进行正确的log1p标准化，进行normalize_total和log1p处理...")
        logger.info("检查是否存在raw数据")
        if adata.raw is not None:
            var_names = adata.var_names
            adata = sc.AnnData(
                X=adata.raw.X, 
                obs=adata.obs,
                var=adata.raw.var,
                uns=adata.uns,
                obsm=adata.obsm
            )
            try:
                adata.var_names = var_names
            except:
                pass 
            logger.info("进行normalize_total和log1p处理...")
            sc.pp.normalize_total(adata, target_sum=1e4)
            sc.pp.log1p(adata)
        else:
            logger.info("未找到raw数据,请核查数据是否正确")
            exit(0)
    else:
        logger.info("数据已经进行了正确的log1p标准化")  
    
    return adata

def prepare_clusters(adata, cluster_key=None, resolution=0.8):
    """准备聚类结果
    
    参数:
        adata: AnnData对象
        cluster_key: 已有的聚类结果列名，如果为None则重新聚类
        resolution: 聚类分辨率
    """
    if cluster_key is None or cluster_key not in adata.obs.columns:
        logger.info("未指定有效的聚类结果，进行新的聚类...")
        # 计算PCA
        if 'X_pca' not in adata.obsm:
            # 创建临时副本进行处理，保持原始.X不变
            temp_adata = adata.copy()
            
            # 过滤表达细胞数少于10的基因
            logger.info("过滤表达细胞数少于10的基因...")
            sc.pp.filter_genes(temp_adata, min_cells=10)
            logger.info(f"过滤后剩余基因数: {temp_adata.n_vars}")
            
            sc.pp.highly_variable_genes(temp_adata, n_top_genes=2000)
            sc.pp.scale(temp_adata)
            sc.tl.pca(temp_adata, svd_solver='arpack')
            # 将计算结果复制回原始对象
            adata.obsm['X_pca'] = temp_adata.obsm['X_pca']
            adata.var['highly_variable'] = pd.Series(False, index=adata.var_names)
            adata.var.loc[temp_adata.var_names, 'highly_variable'] = temp_adata.var['highly_variable']
            del temp_adata
        
        # 计算邻居图
        if 'neighbors' not in adata.uns:
            temp_adata = adata.copy()
            sc.pp.neighbors(temp_adata, n_neighbors=15, n_pcs=30)
            # 复制邻居图信息回原始对象
            adata.uns['neighbors'] = temp_adata.uns['neighbors']
            adata.obsp['distances'] = temp_adata.obsp['distances']
            adata.obsp['connectivities'] = temp_adata.obsp['connectivities']
            del temp_adata
        
        # 计算UMAP
        if 'X_umap' not in adata.obsm:
            temp_adata = adata.copy()
            sc.tl.umap(temp_adata)
            # 复制UMAP结果回原始对象
            adata.obsm['X_umap'] = temp_adata.obsm['X_umap']
            del temp_adata
        
        # 执行leiden聚类
        sc.tl.leiden(adata, resolution=resolution)
        cluster_key = 'leiden'
        logger.info(f"完成聚类，共{len(adata.obs[cluster_key].unique())}个簇")
    else:
        logger.info(f"使用已有的聚类结果: {cluster_key}")
        
        # 如果没有UMAP，也计算一下
        if 'X_umap' not in adata.obsm:
            # 确保有PCA和neighbors
            if 'X_pca' not in adata.obsm:
                temp_adata = adata.copy()
                
                # 过滤表达细胞数少于10的基因
                logger.info("过滤表达细胞数少于10的基因...")
                sc.pp.filter_genes(temp_adata, min_cells=10)
                logger.info(f"过滤后剩余基因数: {temp_adata.n_vars}")
                
                sc.pp.highly_variable_genes(temp_adata, n_top_genes=2000)
                sc.pp.scale(temp_adata)
                sc.tl.pca(temp_adata, svd_solver='arpack')
                adata.obsm['X_pca'] = temp_adata.obsm['X_pca']
                adata.var['highly_variable'] = pd.Series(False, index=adata.var_names)
                adata.var.loc[temp_adata.var_names, 'highly_variable'] = temp_adata.var['highly_variable']
                del temp_adata
            
            if 'neighbors' not in adata.uns:
                temp_adata = adata.copy()
                sc.pp.neighbors(temp_adata, n_neighbors=15, n_pcs=30)
                adata.uns['neighbors'] = temp_adata.uns['neighbors']
                adata.obsp['distances'] = temp_adata.obsp['distances']
                adata.obsp['connectivities'] = temp_adata.obsp['connectivities']
                del temp_adata
            
            # 计算UMAP
            temp_adata = adata.copy()
            sc.tl.umap(temp_adata)
            adata.obsm['X_umap'] = temp_adata.obsm['X_umap']
            del temp_adata
            logger.info("已计算UMAP降维")
    
    return cluster_key

def run_annotation(indata, 
                   model,
                   outdir,
                   update_models=False,
                   show_models=False,
                   mode="best_match", 
                   majority_voting=False, 
                   p_thres=0.5,
                   min_prop=0.0, 
                   prefix="", 
                   xlsx=False,
                   plot_results=False, 
                   use_gpu=False,
                   over_clustering='auto' ):
    """运行细胞类型注释"""
 
    if update_models:
        logger.info("更新Celltypist模型...")
        celltypist.models.download_models(force_update=True)
        exit(0)
 
    if show_models:
        logger.info("显示所有Celltypist模型...")
        md = celltypist.models.models_description()
        for _, row in md.iterrows():
            row = row.tolist()
            logger.info(row[0] + '   ' + row[1])
        exit(0)

    #validate model
    if model is None:
        model = celltypist.models.get_default_model()
        logger.info(f"🔖 未指定模型，使用默认模型: '{model}'")
    if model not in celltypist.models.get_all_models() and not os.path.exists(model):
        show_help_and_exit(f"🛑 无效模型名称: '{model}'. 可用模型: {', '.join(celltypist.models.get_all_models())}")
    
    # 检查模型特征数与输入数据基因数的兼容性
    logger.info("检查模型特征数与输入数据基因数的兼容性...")
    model_obj = celltypist.models.Model.load(model)
    model_genes = set(model_obj.features)
    data_genes = set(indata.var_names)
    common_genes = model_genes.intersection(data_genes)
    
    logger.info(f"模型特征数: {len(model_genes)}")
    logger.info(f"输入数据基因数: {len(data_genes)}")
    logger.info(f"共有基因数: {len(common_genes)}")
    
    overlap_percent = len(common_genes) / len(model_genes) * 100
    logger.info(f"基因重叠率: {overlap_percent:.2f}%")
    
    if overlap_percent < 50:
        logger.warning(f"⚠️ 警告: 基因重叠率低于50%，可能影响注释准确性")
    
    #output dir
    if outdir is None:
        outdir = os.getcwd()
        logger.warn(f"👀 输出目录未指定，使用当前目录: '{outdir}'")
    
    # 处理输出目录
    outdir = Path(str(outdir))
    outdir.mkdir(parents=True, exist_ok=True)  # 创建输出目录
    
    # 设置文件路径

    fig_dir = outdir / 'figures'  # 图片目录
    fig_dir.mkdir(exist_ok=True)
    
    #config settings
    config = {
        "indata": indata,
        "model": model,
        "mode": mode,
        "p-thres": p_thres,
        "majority-voting": majority_voting,
        "outdir": str(outdir),
        "prefix": prefix,
        "xlsx": xlsx,
        "plot-results": plot_results
    }
    if majority_voting:
        config["over-clustering"] = over_clustering if over_clustering != 'auto' else None
        config["use-gpu"] = use_gpu
        config["min-prop"] = min_prop

    #quiet or not
    show_config(config)
    #celltyping and majority voting
    logger.info("开始细胞类型预测...")
    t_start = time.time()

    # 处理over_clustering参数
    if over_clustering == 'auto':
        over_clustering = None

    result = celltypist.annotate(
                                filename=indata,
                                model=model,
                                mode=mode.replace("_", " "),
                                majority_voting=majority_voting,
                                over_clustering=over_clustering,  # 现在如果是'auto'会被设置为None
                                p_thres=p_thres,
                                min_prop=min_prop,
                                use_GPU=use_gpu
                                )
    t_end = time.time()
    logger.info(f"预测完成，耗时: {(t_end - t_start)/60:.2f} 分钟")
    
    # 将预测结果添加到adata中
    result.to_adata(indata)
    if majority_voting:
        indata.obs['celltypist_cell_type']=indata.obs['majority_voting'] 
    else:
        indata.obs['celltypist_cell_type']=indata.obs['predicted_labels'] 
    # 保存其他结果
    result.to_table(folder=str(outdir), prefix=prefix, xlsx=xlsx)
    
    # 生成可视化结果
    if plot_results:
        logger.info("生成可视化结果...")
        result.to_plots(
            folder=str(fig_dir),
            prefix=prefix,
            plot_probability=True,  # 只保留支持的参数
            format='pdf'  # 添加format参数
        )
        logger.info(f"可视化结果已保存到: {fig_dir}")
    
    return indata

def visualize_results(adata, 
                     reduction="X_umap",
                     cell_type_key="major_type",
                     cluster_key="clusters",
                     palette="customecol2",
                     outdir="results"):
    """可视化结果"""
    ## 输出目录检查
    os.makedirs(outdir, exist_ok=True)
    ## 绘制预测细胞类型图
    logger.info(f"绘制预测细胞类型({cell_type_key})图: {outdir}/{reduction}_{cell_type_key}.pdf/png")
    ## 对细胞类型进行排序并转换为category类型
    if cell_type_key in adata.obs.columns:
        # 获取排序后的唯一细胞类型
        unique_types = sorted(adata.obs[cell_type_key].unique())
        adata.obs[cell_type_key] = pd.Categorical(
            adata.obs[cell_type_key],
            categories=unique_types
        )
        cell_type_colors = select_colors(object=adata.obs, value=cell_type_key, palette=palette)
 
        # 添加颜色列到obs
        adata.obs[f"{cell_type_key}_col"] = adata.obs[cell_type_key].map(cell_type_colors )
        fig = sc.pl.umap(adata, color=cell_type_key, palette=cell_type_colors, show=False)
    else:
        fig = sc.pl.umap(adata, color=cell_type_key, show=False)
    
    plt.savefig(f'{outdir}/{reduction}_{cell_type_key}.pdf', bbox_inches='tight', dpi=300)
    plt.savefig(f'{outdir}/{reduction}_{cell_type_key}.png', bbox_inches='tight', dpi=300)
    plt.close()
 
    ## 保存聚类汇总
    if cluster_key:
        logger.info("\n各簇细胞类型组成:")
        # 创建聚类与细胞类型的交叉表
        cluster_celltype = pd.crosstab(
            adata.obs[cluster_key], 
            adata.obs['predicted_labels']
        )
        
        # 计算每个簇中各细胞类型的比例
        cluster_celltype_pct = cluster_celltype.div(cluster_celltype.sum(axis=1), axis=0) * 100
        
        # 合并数量和百分比
        detailed_summary = []
        for cluster in cluster_celltype.index:
            cluster_total = cluster_celltype.loc[cluster].sum()
            
            # 获取该簇中所有细胞类型及其数量和比例
            for cell_type in cluster_celltype.columns:
                count = cluster_celltype.loc[cluster, cell_type]
                percentage = cluster_celltype_pct.loc[cluster, cell_type]
                
                # 只添加数量大于0的细胞类型
                if count > 0:
                    detailed_summary.append({
                        'Cluster': cluster,
                        'Cell_Type': cell_type,
                        'Count': count,
                        'Percentage': percentage,
                        'Total_Cells': cluster_total
                    })
        
        # 创建详细汇总数据框
        detailed_df = pd.DataFrame(detailed_summary)
        
        # 按簇和百分比降序排序
        detailed_df = detailed_df.sort_values(['Cluster', 'Percentage'], ascending=[True, False])
        
        # 保存详细汇总
        detailed_path = f"{outdir}/cluster_celltype_composition.csv"
        detailed_df.to_csv(detailed_path, index=False)
        logger.info(f"簇细胞类型组成详细信息已保存到: {detailed_path}")
        
        # 同时保存原来的主要类型汇总（保持向后兼容）
        cluster_summary = []
        for cluster in cluster_celltype.index:
            major_type = cluster_celltype.loc[cluster].nlargest(1)
            pct = major_type.values[0] / cluster_celltype.loc[cluster].sum() * 100
            cluster_summary.append({
                'Cluster': cluster,
                'Major_type': major_type.index[0],
                'Percentage': pct
            })
        summary_df = pd.DataFrame(cluster_summary)
        summary_path = f"{outdir}/cluster_summary.csv"
        summary_df.to_csv(summary_path, index=False)
        logger.info(f"簇注释汇总已保存到: {summary_path}")

    return adata
@click.command()
@click.option('--input', required=True,
              help='输入文件路径，支持h5ad,h5seurat格式')
@click.option('--model-path', required=True, 
              help='Celltypist模型路径，pkl格式')
@click.option('--output-path', required=True, 
              help='结果保存路径，将保存为h5ad格式。如果未指定扩展名，将自动添加.h5ad')

@click.option('--cluster-key', default="clusters",
              show_default=True,
              help='已有的常规聚类结果列名，在最终汇总统计时候需要，注意该值并不用于过度聚类')
@click.option('--resolution', default=0.8, type=float,
              show_default=True,
              help='聚类分辨率，用于常规聚类，注意该值与过度聚类无关，主要用于cluster-key参数设置为空时，采用scanpy进行leiden聚类')

@click.option('--mode', default="best_match", 
              type=click.Choice(['best_match', 'prob_match']),
              show_default=True,
              help='执行细胞预测的方式：'
                   '默认模式best match会选择得分/概率最高的细胞类型作为最终预测结果。'
                   '设置为prob match将启用多标签分类，可以为每个查询细胞分配0个（即未分配）、1个或2个及以上的细胞类型标签')
@click.option('--p-thres', default=0.5, type=float,
              show_default=True,
              help='预测概率阈值，默认为0.5，仅在mode为prob match时使用')

@click.option("--majority-voting",
              is_flag=True, 
              default=False, 
              show_default=True,
              help="在后运行多数投票分类器以细化预测标签。")
@click.option("--over-clustering", 
              default='auto', 
              show_default=True,
              help="该参数和majority voting一起使用，默认是auto，celltypist将自动进行过度聚类，"
              "也可以指定obs中对应的过度聚类列名，但考虑到我们通常不会进行过度聚类，所以默认使用是auto，"
              "由软件自动进行过度聚类。",
              type=str)
@click.option('--use-gpu',  default=False, show_default=True,help='是否使用GPU，用于过度聚类scanpy加速，目前尚不支持，默认为False')
@click.option('--min-prop',default=0, type=float, show_default=True,help='majority voting时，簇内主要类型的最小占比阈值，默认0')
@click.option('--prefix', default="",show_default=True, type=str, help='输出文件的前缀,默认空')
@click.option('--xlsx', default=False,show_default=True, help='是否将结果合并为Excel文件,默认False')
@click.option('--plot-results',  default=False,show_default=True, help='是否生成可视化结果,默认False')
@click.option('--update-models',  default=False,show_default=True, help='是否更新Celltypist模型')
@click.option('--show-models',  default=False,show_default=True, help='是否显示所有Celltypist模型')
@click.option('--reduction', default="X_umap", show_default=True, help='UMAP降维方法，默认Xumap,用于可视化结果')

def main(input, model_path, output_path, cluster_key, resolution,
         mode, majority_voting, p_thres, min_prop, prefix, 
         xlsx, plot_results, update_models, show_models, use_gpu,
         over_clustering, reduction):
    """Celltypist细胞类型注释脚本

    此脚本用于使用Celltypist模型进行细胞类型注释。

    主要功能包括：
    1. 数据预处理：标准化和log转换
    2. 可选的聚类分析
    3. 使用指定模型进行细胞类型预测
    4. 保存注释结果和统计信息

    示例用法：
    python run_celltypist_annotate.py  
        --input input.h5seurat  
        --model-path model.pkl  
        --output-path output.h5ad  
        --mode best_match  
        --majority-voting   
    """

    logger.info("step1:加载并预处理数据=======================================")
    adata = load_and_preprocess_data(input)
    logger.info(f"数据形状: {adata.shape}") 
    print("输入数据基因名称示例:", list(adata.var_names[:5]))

    logger.info("step2:如果cluster_key参数为空，则采用scanpy进行聚类===========")
    if cluster_key is not None:
        cluster_key = prepare_clusters(adata, cluster_key, resolution)

    logger.info("step3:运行细胞类型注释=======================================")
    adata = run_annotation(
        indata=adata,
        model=model_path,
        outdir=output_path,
        mode=mode,
        majority_voting=majority_voting,
        p_thres=p_thres,
        min_prop=min_prop,
        prefix=prefix,
        xlsx=xlsx,
        plot_results=plot_results,
        update_models=update_models,
        show_models=show_models,
        use_gpu=use_gpu,
        over_clustering=over_clustering
    )

    logger.info("step4:可视化结果=============================================")
    adata=visualize_results(adata, 
                            reduction=reduction, 
                            cell_type_key="celltypist_cell_type",
                            cluster_key=cluster_key, 
                            outdir=output_path)

    logger.info("step5:保存对象==============================================")
    if input.endswith('.h5ad'):
        # 保存h5ad文件
        h5ad_file = f"{output_path}/adata.h5ad"  # h5ad文件
        logger.info(f"保存注释结果到: {h5ad_file}")
        # 'predicted_labels', 'over_clustering','majority_voting', 'conf_score'
        # del adata.obs['predicted_labels']
        # del adata.obs['over_clustering']
        # del adata.obs['majority_voting']
        # del adata.obs['conf_score']
        adata.write_h5ad(str(h5ad_file))
        metadata = adata.obs
        metadata.insert(0, "Barcode", metadata.index)
        metadata.to_csv(f"{output_path}/metadata.tsv", index=False, sep="\t")

    else:
        # 将adata.obs中特定列更新至seurat对象中
        h5seurat_file = f"{output_path}/adata.h5ad"
        update_h5seurat_metadata(adata=adata,
                                 h5seurat_path=input,
                                 output_h5seurat=h5seurat_file,
                                 col_names=['predicted_labels',
                                            'celltypist_cell_type',
                                            'celltypist_cell_type_col']  )
if __name__ == '__main__':
    main() 
