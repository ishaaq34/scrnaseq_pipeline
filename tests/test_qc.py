"""
Unit tests for QC functions
"""

import pytest
import scanpy as sc
import numpy as np
import pandas as pd


def create_test_adata():
    """Create a test AnnData object"""
    np.random.seed(42)
    n_cells = 100
    n_genes = 50
    
    # Create random count matrix
    X = np.random.negative_binomial(5, 0.3, (n_cells, n_genes))
    
    # Create AnnData object
    adata = sc.AnnData(X)
    adata.var_names = [f"Gene_{i}" for i in range(n_genes)]
    adata.obs_names = [f"Cell_{i}" for i in range(n_cells)]
    
    # Add mitochondrial genes
    mt_genes = np.random.choice(n_genes, 5, replace=False)
    adata.var['mt'] = False
    adata.var.iloc[mt_genes, adata.var.columns.get_loc('mt')] = True
    
    return adata


def test_qc_metrics_calculation():
    """Test QC metrics calculation"""
    adata = create_test_adata()
    
    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=['mt'],
        percent_top=None,
        log1p=False,
        inplace=True
    )
    
    assert 'n_genes_by_counts' in adata.obs.columns
    assert 'total_counts' in adata.obs.columns
    assert 'pct_counts_mt' in adata.obs.columns


def test_cell_filtering():
    """Test cell filtering based on QC thresholds"""
    adata = create_test_adata()
    
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], inplace=True)
    
    n_cells_before = adata.n_obs
    
    # Apply filters
    sc.pp.filter_cells(adata, min_genes=10)
    adata = adata[adata.obs.pct_counts_mt < 50, :]
    
    n_cells_after = adata.n_obs
    
    assert n_cells_after <= n_cells_before


def test_gene_filtering():
    """Test gene filtering based on expression"""
    adata = create_test_adata()
    
    n_genes_before = adata.n_vars
    
    sc.pp.filter_genes(adata, min_cells=3)
    
    n_genes_after = adata.n_vars
    
    assert n_genes_after <= n_genes_before


def test_normalization():
    """Test normalization"""
    adata = create_test_adata()
    
    # Log normalization
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    
    # Check that data is normalized
    assert adata.X.max() > 0
    assert np.isfinite(adata.X).all()


def test_highly_variable_genes():
    """Test highly variable gene selection"""
    adata = create_test_adata()
    
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=20)
    
    assert 'highly_variable' in adata.var.columns
    assert adata.var['highly_variable'].sum() == 20


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
