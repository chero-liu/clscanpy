from typing import Optional, Union
import logging
import doubletdetection
import scanpy as sc
import numpy as np
import pandas as pd
from concurrent.futures import ThreadPoolExecutor, as_completed
from math import ceil
import random

random.seed(2025)
# Set up logging
import logging

LOGGER = logging.getLogger(__name__)


def do_doubletdetection(
    adata: sc.AnnData,
    pseudocount: float,
    n_iters: int = 10,
    clustering_algorithm: str = "louvain",
    standard_scaling: bool = True,
    n_jobs: int = -1,
) -> sc.AnnData:
    """Detect doublets using doubletdetection method

    Args:
        adata: AnnData object containing single-cell data
        pseudocount: Expected doublet rate
        config: Configuration parameters for the classifier

    Returns:
        AnnData object with doublet predictions in adata.obs['doublet']
    """
    LOGGER.info(
        f"Running doubletdetection with n_iters={n_iters}, clustering={clustering_algorithm}"
    )

    clf = doubletdetection.BoostClassifier(
        n_iters=n_iters,
        clustering_algorithm=clustering_algorithm,
        standard_scaling=standard_scaling,
        pseudocount=pseudocount,
        n_jobs=n_jobs,
    )

    LOGGER.debug("Classifier initialized, starting prediction...")

    doublets = clf.fit(adata.X).predict()
    adata.obs["doublet"] = pd.Series(doublets, index=adata.obs.index).replace(
        {1: True, 0: False}
    )

    return adata


def do_scrublet(
    adata: sc.AnnData,
    expected_doublet_rate: float,
) -> sc.AnnData:
    """Detect doublets using scrublet method

    Args:
        adata: AnnData object containing single-cell data
        expected_doublet_rate: Expected doublet rate

    Returns:
        AnnData object with doublet predictions in adata.obs['doublet']
    """
    LOGGER.info(f"Running scrublet with expected doublet rate: {expected_doublet_rate}")
    adata = sc.pp.scrublet(
        adata, expected_doublet_rate=expected_doublet_rate, copy=True
    )
    adata.obs = adata.obs.drop(columns="doublet_score")
    adata.obs = adata.obs.rename(columns={"predicted_doublet": "doublet"})
    return adata


def calculate_doublet_rates(cell_count: int, version: str) -> float:
    """Calculate expected doublet rate based on cell count

    Args:
        cell_count: Number of cells in the sample

    Returns:
        Expected doublet rate
    """
    rate_rangesv3 = [
        (0, 1000, 0.008),
        (1001, 2000, 0.015),
        (2001, 3000, 0.023),
        (3001, 4000, 0.030),
        (4001, 5000, 0.038),
        (5001, 6000, 0.046),
        (6001, 7000, 0.053),
        (7001, 8000, 0.061),
        (8001, 9000, 0.068),
        (9001, 10000, 0.080),
    ]
    rate_rangesv4 = [
        (0, 1000, 0.004),
        (1001, 2000, 0.008),
        (2001, 3000, 0.012),
        (3001, 4000, 0.016),
        (4001, 5000, 0.020),
        (5001, 6000, 0.024),
        (6001, 7000, 0.028),
        (7001, 8000, 0.032),
        (8001, 9000, 0.036),
        (9001, 10000, 0.040),
        (10001, 11000, 0.044),
        (11001, 12000, 0.048),
        (12001, 13000, 0.052),
        (13001, 14000, 0.056),
        (14001, 15000, 0.060),
        (15001, 16000, 0.064),
        (16001, 17000, 0.068),
        (17001, 18000, 0.072),
        (18001, 19000, 0.076),
        (19001, 20000, 0.080),
    ]
    rate_rangesC4 = [(0, 6000, 0.0119), (6001, 12000, 0.0261), (12001, 18000, 0.0406)]
    if version == "v4":
        LOGGER.info("version is v4")
        for min_cells, max_cells, rate in rate_rangesv4:
            if min_cells <= cell_count <= max_cells:
                LOGGER.debug(
                    f"Cell count {cell_count} falls in range {min_cells}-{max_cells}, using rate {rate}"
                )
                return rate
        LOGGER.debug(
            f"Cell count {cell_count} exceeds maximum range, using default rate 0.1"
        )
        return 0.1
    elif version == "C4":
        LOGGER.info("version is C4")
        for min_cells, max_cells, rate in rate_rangesC4:
            if min_cells <= cell_count <= max_cells:
                LOGGER.debug(
                    f"Cell count {cell_count} falls in range {min_cells}-{max_cells}, using rate {rate}"
                )
                return rate
        LOGGER.debug(f"Cell count {cell_count} exceeds 18000, using rate 0.0637")
        return 0.0637
    else:
        LOGGER.info("version is v3")
        for min_cells, max_cells, rate in rate_rangesv3:
            if min_cells <= cell_count <= max_cells:
                LOGGER.debug(
                    f"Cell count {cell_count} falls in range {min_cells}-{max_cells}, using rate {rate}"
                )
                return rate
        LOGGER.debug(
            f"Cell count {cell_count} exceeds maximum range, using default rate 0.1"
        )
        return 0.1


def _process_sample(
    adata: sc.AnnData,
    sample: str,
    batch_key: str,
    method: str,
    dbl_rates: float,
    n_iters: int,
    clustering_algorithm: str,
    standard_scaling: bool,
    n_jobs: int,
) -> sc.AnnData:
    """Process a single sample for doublet detection"""
    LOGGER.debug(f"Processing sample: {sample}")
    sample_mask = adata.obs[batch_key] == sample
    subadata = adata[sample_mask, :].copy()

    if method == "doubletdetection":
        subadata = do_doubletdetection(
            subadata,
            dbl_rates,
            n_iters=n_iters,
            clustering_algorithm=clustering_algorithm,
            standard_scaling=standard_scaling,
            n_jobs=n_jobs,
        )
    elif method == "scrublet":
        subadata = do_scrublet(subadata, expected_doublet_rate=dbl_rates)
    else:
        raise ValueError(
            f"Invalid method '{method}'. Choose 'doubletdetection' or 'scrublet'."
        )

    return subadata


def find_doublet(
    adata: sc.AnnData,
    method: str,
    batch_key: str = "sampleid",
    n_iters: int = 10,
    clustering_algorithm: str = "louvain",
    standard_scaling: bool = True,
    n_jobs: int = -1,
    max_workers: int = 10,
    version: str = "v3",
) -> pd.Series:
    """Detect doublets in single-cell data with parallel processing

    Args:
        adata: AnnData object containing single-cell data
        method: Detection method ('doubletdetection' or 'scrublet')
        batch_key: Column in adata.obs containing batch information
        n_iters: Number of iterations for doubletdetection
        clustering_algorithm: Clustering algorithm to use
        standard_scaling: Whether to standard scale the data
        n_jobs: Number of jobs for parallel processing
        max_workers: Maximum number of parallel workers

    Returns:
        Series containing doublet predictions
    """
    samples = adata.obs[batch_key].unique()
    total_samples = len(samples)
    workers = max_workers
    batches = ceil(total_samples / workers)

    LOGGER.info(
        f"Processing {total_samples} samples with method: {method} in {batches} batches of {workers} workers"
    )

    adatas = []
    for batch_num in range(batches):
        start = batch_num * workers
        end = min((batch_num + 1) * workers, total_samples)
        batch_samples = samples[start:end]

        LOGGER.debug(
            f"Processing batch {batch_num + 1}/{batches} with samples: {batch_samples}"
        )

        with ThreadPoolExecutor(max_workers=workers) as executor:
            futures = {}
            for sample in batch_samples:
                sample_mask = adata.obs[batch_key] == sample
                subadata = adata[sample_mask, :].copy()
                dbl_rates = calculate_doublet_rates(subadata.n_obs, version)
                future = executor.submit(
                    _process_sample,
                    adata=adata,
                    sample=sample,
                    batch_key=batch_key,
                    method=method,
                    dbl_rates=dbl_rates,
                    n_iters=n_iters,
                    clustering_algorithm=clustering_algorithm,
                    standard_scaling=standard_scaling,
                    n_jobs=n_jobs,
                )
                futures[future] = sample

            for future in as_completed(futures):
                sample = futures[future]
                try:
                    subadata = future.result()
                    adatas.append(subadata)
                    LOGGER.info(
                        f"Sample {sample} processed: {subadata.obs['doublet'].value_counts().to_dict()}"
                    )
                except Exception as e:
                    LOGGER.error(f"Error processing sample {sample}: {str(e)}")
                    raise RuntimeError(f"Failed to process sample {sample}") from e

    adata = sc.concat(adatas, join="outer")
    return adata.obs["doublet"]
