"""
Perturbation embeddings for scPAGE

Creates perturbation-level transcriptional signatures and low-dimensional
embeddings from the normalized gene-by-cell matrix produced by
perturbation_preprocessing.py

A perturbation signature is defined as:

    mean(expression | perturbation, condition)
    -
    mean(expression | controls, condition)

Replicates are aggregated separately before being combined so that
replicates containing more cells do not dominate the signature.

Outputs
-------
perturbation_signatures.csv
    Rows are perturbations, columns are genes. Values are control-centered
    transcriptional responses.
perturbation_embeddings.csv
    PCA coordinates and, when umap-learn is available, UMAP coordinates

perturbation_similarity.csv
    Cosine similarity matrix between perturbation signatures

perturbation_neighbors.csv
    Nearest perturbations for every perturbation

Perturbation_embeddings.png
    Two-dimensional perturbation map
"""

import os
from typing import Iterable, Optional, Union

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from sklearn.decomposition import PCA
from sklearn.metrics.pairwise import cosine_similarity
from sklearn.preprocessing import StandardScaler


META_ROWS = ("condition", "replicate", "gene")


def _load_matrix(
    matrix: Union[str, os.PathLike, pd.DataFrame]
) -> pd.DataFrame:
    """
    Accept either the in-memory scPAGE perturbation matrix or a CSV path.
    """
    if isinstance(matrix, pd.DataFrame):
        return matrix.copy()

    if isinstance(matrix, (str, os.PathLike)):
        return pd.read_csv(matrix, index_col=0)

    raise TypeError(
        "MATRIX must be a pandas DataFrame or path to a CSV file."
    )


def _prepare_cell_table(
    matrix: pd.DataFrame,
    gene_subset: Optional[Iterable[str]] = None,
) -> tuple[pd.DataFrame, list[str]]:
    """
    Convert scPAGE's gene-by-cell matrix into a cell-by-feature table.

    Expected matrix shape
    ---------------------
                    cell1      cell2      cell3
    condition       day7       day7       day7
    replicate       rep1       rep1       rep2
    gene             TP53    Control       TP53
    ACTB             2.1        2.0        2.2
    GAPDH            4.3        4.1        4.4
    ...
    """
    missing = [row for row in META_ROWS if row not in matrix.index]

    if missing:
        raise ValueError(
            "Input matrix does not contain required metadata rows: "
            + ", ".join(missing)
        )

    expression_rows = [
        idx for idx in matrix.index
        if idx not in META_ROWS and idx not in ("cell", "grna")
    ]

    if gene_subset is not None:
        requested = set(gene_subset)
        expression_rows = [
            gene for gene in expression_rows
            if gene in requested
        ]

        if not expression_rows:
            raise ValueError(
                "None of the genes in gene_subset were found in MATRIX."
            )

    metadata = matrix.loc[list(META_ROWS)].T.copy()

    expression = (
        matrix.loc[expression_rows]
        .apply(pd.to_numeric, errors="coerce")
        .T
    )

    # Keep metadata and expression indexed identically.
    metadata.index = np.arange(len(metadata))
    expression.index = metadata.index

    cell_table = pd.concat([metadata, expression], axis=1)

    cell_table = cell_table.dropna(
        subset=["condition", "replicate", "gene"]
    )

    return cell_table, expression_rows


def build_perturbation_signatures(
    MATRIX: Union[str, os.PathLike, pd.DataFrame],
    gene_subset: Optional[Iterable[str]] = None,
    control_label: str = "Control",
    aggregation: str = "median",
    min_cells: int = 3,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Build control-centered perturbation signatures.

    The calculation occurs in two stages:

    1. aggregate cells within
       condition x replicate x perturbation
    2. average replicate-level profiles

    The matched control profile for each condition is then subtracted.

    Parameters
    ----------
    MATRIX
        Output of perturbation_preprocessing.combine_count_table(), or
        a CSV containing the same matrix.

    gene_subset
        Optional iterable of genes to use as embedding features.
        `va_genes` from perturbation.most_variable_genes() can be passed
        directly here.

    control_label
        Label assigned to control sgRNAs. scPAGE currently writes
        controls as "Control".

    aggregation
        "mean" or "median" cell-level aggregation.

    min_cells
        Minimum number of cells required for a perturbation in a replicate.

    Returns
    -------
    signatures
        DataFrame indexed by perturbation ID and containing one column
        per expression gene.

    metadata
        Perturbation-level metadata including cell and replicate counts.
    """
    matrix = _load_matrix(MATRIX)

    cell_table, expression_genes = _prepare_cell_table(
        matrix,
        gene_subset=gene_subset,
    )

    if aggregation not in {"mean", "median"}:
        raise ValueError(
            "aggregation must be either 'mean' or 'median'."
        )

    group_cols = ["condition", "replicate", "gene"]

    counts = (
        cell_table
        .groupby(group_cols, observed=True)
        .size()
        .rename("n_cells")
        .reset_index()
    )

    valid = counts[counts["n_cells"] >= min_cells][group_cols]

    cell_table = cell_table.merge(
        valid,
        on=group_cols,
        how="inner",
    )

    if cell_table.empty:
        raise ValueError(
            "No perturbations remain after applying min_cells="
            f"{min_cells}."
        )

    grouped = cell_table.groupby(
        group_cols,
        observed=True,
    )[expression_genes]

    if aggregation == "median":
        replicate_profiles = grouped.median()
    else:
        replicate_profiles = grouped.mean()

    replicate_profiles = replicate_profiles.reset_index()

    # Replicate profiles are now weighted equally.
    perturbation_profiles = (
        replicate_profiles
        .groupby(["condition", "gene"], observed=True)[expression_genes]
        .mean()
    )

    signatures = []
    metadata_records = []

    conditions = perturbation_profiles.index.get_level_values(
        "condition"
    ).unique()

    for condition in conditions:

        condition_profiles = perturbation_profiles.loc[condition]

        if control_label not in condition_profiles.index:
            raise ValueError(
                f"No control profile labelled '{control_label}' "
                f"was found for condition '{condition}'."
            )

        control_profile = condition_profiles.loc[control_label]

        perturbations = condition_profiles.drop(
            index=control_label,
            errors="ignore",
        )

        for perturbed_gene, profile in perturbations.iterrows():

            signature = profile - control_profile

            perturbation_id = (
                f"{condition}::{perturbed_gene}"
            )

            signature.name = perturbation_id
            signatures.append(signature)

            count_subset = counts[
                (counts["condition"] == condition)
                & (counts["gene"] == perturbed_gene)
                & (counts["n_cells"] >= min_cells)
            ]

            metadata_records.append(
                {
                    "perturbation_id": perturbation_id,
                    "condition": condition,
                    "gene": perturbed_gene,
                    "n_cells": int(count_subset["n_cells"].sum()),
                    "n_replicates": int(len(count_subset)),
                }
            )

    if not signatures:
        raise ValueError(
            "No non-control perturbations were available for embedding."
        )

    signatures = pd.DataFrame(signatures)
    signatures.index.name = "perturbation_id"

    metadata = (
        pd.DataFrame(metadata_records)
        .set_index("perturbation_id")
    )

    return signatures, metadata


def embed_signatures(
    signatures: pd.DataFrame,
    metadata: Optional[pd.DataFrame] = None,
    n_components: int = 10,
    scale_features: bool = True,
    use_umap: bool = True,
    random_state: int = 42,
) -> pd.DataFrame:
    """
    Calculate PCA and optional UMAP embeddings.

    PCA coordinates are always produced.

    UMAP is used when:
      * use_umap=True
      * umap-learn is available
      * at least 3 perturbations exist

    Otherwise UMAP1/UMAP2 fall back to PC1/PC2.
    """
    X = signatures.to_numpy(dtype=float)

    # Genes that never change across perturbations contribute nothing.
    variable = np.nanstd(X, axis=0) > 0

    if not np.any(variable):
        raise ValueError(
            "All perturbation signature features have zero variance."
        )

    X = X[:, variable]

    # Replace any residual NaNs/infs.
    X = np.nan_to_num(
        X,
        nan=0.0,
        posinf=0.0,
        neginf=0.0,
    )

    if scale_features:
        X_embedding = StandardScaler().fit_transform(X)
    else:
        X_embedding = X

    max_components = min(
        n_components,
        X_embedding.shape[0],
        X_embedding.shape[1],
    )

    if max_components < 1:
        raise ValueError(
            "Not enough perturbations/features to compute an embedding."
        )

    pca = PCA(
        n_components=max_components,
        random_state=random_state,
    )

    pcs = pca.fit_transform(X_embedding)

    embedding = pd.DataFrame(
        pcs,
        index=signatures.index,
        columns=[
            f"PC{i + 1}"
            for i in range(pcs.shape[1])
        ],
    )

    embedding["PCA_variance_explained"] = (
        pca.explained_variance_ratio_[0]
        if len(pca.explained_variance_ratio_)
        else np.nan
    )

    # PCA fallback gives a useful result without adding a hard dependency.
    if pcs.shape[1] >= 2:
        embedding["UMAP1"] = pcs[:, 0]
        embedding["UMAP2"] = pcs[:, 1]
        embedding["embedding_method"] = "PCA"
    else:
        embedding["UMAP1"] = pcs[:, 0]
        embedding["UMAP2"] = 0.0
        embedding["embedding_method"] = "PCA"

    if use_umap and X_embedding.shape[0] >= 3:

        try:
            import umap

            # UMAP requires n_neighbors < number of observations.
            n_neighbors = min(
                15,
                max(2, X_embedding.shape[0] - 1),
            )

            reducer = umap.UMAP(
                n_components=2,
                n_neighbors=n_neighbors,
                metric="cosine",
                min_dist=0.25,
                random_state=random_state,
            )

            umap_coordinates = reducer.fit_transform(X_embedding)

            embedding["UMAP1"] = umap_coordinates[:, 0]
            embedding["UMAP2"] = umap_coordinates[:, 1]
            embedding["embedding_method"] = "UMAP"

        except ImportError:
            # scanpy generally brings umap-learn with it, but keeping this
            # optional avoids changing scPAGE's dependency contract.
            pass

    if metadata is not None:
        embedding = metadata.join(embedding)

    return embedding


def perturbation_similarity(
    signatures: pd.DataFrame,
) -> pd.DataFrame:
    """
    Pairwise cosine similarity between perturbation signatures.

    +1 : highly similar transcriptional response
     0 : unrelated / orthogonal response
    -1 : opposing transcriptional response
    """
    X = np.nan_to_num(
        signatures.to_numpy(dtype=float),
        nan=0.0,
        posinf=0.0,
        neginf=0.0,
    )

    similarity = cosine_similarity(X)

    return pd.DataFrame(
        similarity,
        index=signatures.index,
        columns=signatures.index,
    )


def nearest_perturbations(
    similarity: pd.DataFrame,
    metadata: Optional[pd.DataFrame] = None,
    top_n: int = 5,
) -> pd.DataFrame:
    """
    Return nearest neighbours for each perturbation.
    """
    records = []

    for perturbation in similarity.index:

        neighbours = (
            similarity.loc[perturbation]
            .drop(index=perturbation)
            .sort_values(ascending=False)
            .head(top_n)
        )

        for rank, (neighbour, score) in enumerate(
            neighbours.items(),
            start=1,
        ):
            record = {
                "perturbation_id": perturbation,
                "neighbor_id": neighbour,
                "rank": rank,
                "cosine_similarity": float(score),
            }

            if metadata is not None:

                if perturbation in metadata.index:
                    record["condition"] = metadata.loc[
                        perturbation, "condition"
                    ]
                    record["gene"] = metadata.loc[
                        perturbation, "gene"
                    ]

                if neighbour in metadata.index:
                    record["neighbor_condition"] = metadata.loc[
                        neighbour, "condition"
                    ]
                    record["neighbor_gene"] = metadata.loc[
                        neighbour, "gene"
                    ]

            records.append(record)

    return pd.DataFrame(records)


def plot_perturbation_embedding(
    embedding: pd.DataFrame,
    out_file: str,
    label_genes: bool = True,
) -> None:
    """
    Plot perturbations in 2-D embedding space.
    """
    fig, ax = plt.subplots(figsize=(10, 8))

    conditions = embedding["condition"].astype(str).unique()

    for condition in conditions:

        subset = embedding[
            embedding["condition"].astype(str) == condition
        ]

        ax.scatter(
            subset["UMAP1"],
            subset["UMAP2"],
            s=55,
            alpha=0.8,
            label=condition,
        )

    if label_genes:
        for perturbation_id, row in embedding.iterrows():

            ax.annotate(
                str(row["gene"]),
                (row["UMAP1"], row["UMAP2"]),
                xytext=(4, 4),
                textcoords="offset points",
                fontsize=8,
                alpha=0.8,
            )

    method = (
        embedding["embedding_method"].iloc[0]
        if "embedding_method" in embedding.columns
        else "embedding"
    )

    ax.set_xlabel(f"{method} 1")
    ax.set_ylabel(f"{method} 2")
    ax.set_title("scPAGE perturbation embedding")

    if len(conditions) > 1:
        ax.legend(
            title="Condition",
            bbox_to_anchor=(1.02, 1),
            loc="upper left",
        )

    fig.tight_layout()
    fig.savefig(
        out_file,
        dpi=300,
        bbox_inches="tight",
    )
    plt.close(fig)


def perturbation_embeddings(
    MATRIX: Union[str, os.PathLike, pd.DataFrame],
    OUT_DIR: str,
    OUT_PLOT_DIR: Optional[str] = None,
    VA_GENES: Optional[Iterable[str]] = None,
    CONTROL_LABEL: str = "Control",
    AGGREGATION: str = "median",
    MIN_CELLS: int = 3,
    N_COMPONENTS: int = 10,
    TOP_NEIGHBORS: int = 5,
    SCALE_FEATURES: bool = True,
    USE_UMAP: bool = True,
    RANDOM_STATE: int = 42,
) -> dict:
    """
    scPAGE-compatible high-level perturbation embedding workflow.

    Parameters intentionally follow the uppercase naming style used by
    the existing perturbation functions.

    Returns
    -------
    dict containing:
        signatures
        metadata
        embedding
        similarity
        neighbors
    """
    os.makedirs(OUT_DIR, exist_ok=True)

    if OUT_PLOT_DIR is None:
        OUT_PLOT_DIR = OUT_DIR

    os.makedirs(OUT_PLOT_DIR, exist_ok=True)

    signatures, metadata = build_perturbation_signatures(
        MATRIX=MATRIX,
        gene_subset=VA_GENES,
        control_label=CONTROL_LABEL,
        aggregation=AGGREGATION,
        min_cells=MIN_CELLS,
    )

    embedding = embed_signatures(
        signatures=signatures,
        metadata=metadata,
        n_components=N_COMPONENTS,
        scale_features=SCALE_FEATURES,
        use_umap=USE_UMAP,
        random_state=RANDOM_STATE,
    )

    similarity = perturbation_similarity(signatures)

    neighbors = nearest_perturbations(
        similarity=similarity,
        metadata=metadata,
        top_n=TOP_NEIGHBORS,
    )

    signatures.to_csv(
        os.path.join(
            OUT_DIR,
            "perturbation_signatures.csv",
        )
    )

    embedding.to_csv(
        os.path.join(
            OUT_DIR,
            "perturbation_embeddings.csv",
        )
    )

    similarity.to_csv(
        os.path.join(
            OUT_DIR,
            "perturbation_similarity.csv",
        )
    )

    neighbors.to_csv(
        os.path.join(
            OUT_DIR,
            "perturbation_neighbors.csv",
        ),
        index=False,
    )

    plot_perturbation_embedding(
        embedding=embedding,
        out_file=os.path.join(
            OUT_PLOT_DIR,
            "Perturbation_embeddings.png",
        ),
    )

    return {
        "signatures": signatures,
        "metadata": metadata,
        "embedding": embedding,
        "similarity": similarity,
        "neighbors": neighbors,
    }