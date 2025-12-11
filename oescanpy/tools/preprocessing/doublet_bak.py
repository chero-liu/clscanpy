from typing import Optional, Union
import logging
import doubletdetection
import scanpy as sc
import numpy as np
import pandas as pd

# Set up logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def do_doubletdetection(
    adata: sc.AnnData, 
    pseudocount: float,
    n_iters: int = 10,
    clustering_algorithm: str = "louvain",
    standard_scaling: bool = True,
    n_jobs: int = -1
) -> sc.AnnData:
    """Detect doublets using doubletdetection method
    
    Args:
        adata: AnnData object containing single-cell data
        pseudocount: Expected doublet rate
        config: Configuration parameters for the classifier
        
    Returns:
        AnnData object with doublet predictions in adata.obs['doublet']
    """
    logger.info(f"Running doubletdetection with n_iters={n_iters}, clustering={clustering_algorithm}")
    
    clf = doubletdetection.BoostClassifier(
        n_iters=n_iters,
        clustering_algorithm=clustering_algorithm,
        standard_scaling=standard_scaling,
        pseudocount=pseudocount,
        n_jobs=n_jobs
    )
    
    logger.debug("Classifier initialized, starting prediction...")
    
    doublets = clf.fit(adata.X).predict()
    adata.obs["doublet"] = pd.Series(doublets, index=adata.obs.index).replace({1: True, 0: False})

    return adata

def do_scrublet(adata: sc.AnnData, expected_doublet_rate: float) -> sc.AnnData:
    """Detect doublets using scrublet method
    
    Args:
        adata: AnnData object containing single-cell data
        expected_doublet_rate: Expected doublet rate
        
    Returns:
        AnnData object with doublet predictions in adata.obs['doublet']
    """
    logger.info(f"Running scrublet with expected doublet rate: {expected_doublet_rate}")
    adata = sc.pp.scrublet(adata, expected_doublet_rate=expected_doublet_rate, copy=True)
    adata.obs = adata.obs.drop(columns = 'doublet_score')
    adata.obs = adata.obs.rename(columns = {'predicted_doublet' : 'doublet'})
    return adata

def calculate_doublet_rates(cell_count: int) -> float:
    """Calculate expected doublet rate based on cell count
    
    Args:
        cell_count: Number of cells in the sample
        
    Returns:
        Expected doublet rate
    """
    rate_ranges = [
        (0, 1000, 0.008),
        (1001, 2000, 0.015),
        (2001, 3000, 0.023),
        (3001, 4000, 0.030),
        (4001, 5000, 0.038),
        (5001, 6000, 0.046),
        (6001, 7000, 0.053),
        (7001, 8000, 0.061),
        (8001, 9000, 0.068),
        (9001, 10000, 0.080)
    ]
    
    for min_cells, max_cells, rate in rate_ranges:
        if min_cells <= cell_count <= max_cells:
            logger.debug(f"Cell count {cell_count} falls in range {min_cells}-{max_cells}, using rate {rate}")
            return rate
    logger.debug(f"Cell count {cell_count} exceeds maximum range, using default rate 0.1")
    return 0.1

def findDoublet(
    adata: sc.AnnData,
    method: str,
    batch_key: str = "sampleid",
    n_iters: int = 10,
    clustering_algorithm: str = "louvain",
    standard_scaling: bool = True,
    n_jobs: int = -1
) -> pd.Series:
    """Detect doublets in single-cell data
    
    Args:
        adata: AnnData object containing single-cell data
        method: Detection method ('doubletdetection' or 'scrublet')
        batch_key: Column in adata.obs containing batch information
        config: Configuration parameters for doublet detection
        
    Returns:
        Series containing doublet predictions
    """
    adatas = []
    samples = adata.obs[batch_key].unique()
    logger.info(f"Processing {len(samples)} samples with method: {method}")
    
    for sample in samples:
        logger.debug(f"Processing sample: {sample}")
        sample_mask = adata.obs[batch_key] == sample
        subadata = adata[sample_mask, :].copy()
        dbl_rates = calculate_doublet_rates(subadata.n_obs)
        
        try:
            if method == 'doubletdetection':
                subadata = do_doubletdetection(
                    subadata,
                    dbl_rates,
                    n_iters=n_iters,
                    clustering_algorithm=clustering_algorithm,
                    standard_scaling=standard_scaling,
                    n_jobs=n_jobs
                )
            elif method == 'scrublet':
                subadata = do_scrublet(subadata, expected_doublet_rate=dbl_rates)
            else:
                raise ValueError(f"Invalid method '{method}'. Choose 'doubletdetection' or 'scrublet'.")
            
            adatas.append(subadata)
        except Exception as e:
            logger.error(f"Error processing sample {sample}: {str(e)}")
            raise RuntimeError(f"Failed to process sample {sample}") from e
    
    adata = sc.concat(adatas, join="outer")

    return adata.obs['doublet']
