import os

S3_BUCKET = os.getenv("S3_BUCKET")

# Hugging Face dataset identifier
HF_DATASET_NAME = "rmanluo/RoG-webqsp"

# Base S3 path for storing intermediate files and results
S3_BASE_PATH = f"{S3_BUCKET}/{HF_DATASET_NAME.replace('/', '-')}"
S3_EXAMPLES_PATH = f"{S3_BASE_PATH}/examples" # Subdirectory for individual example data

EMBEDDINGS_PROVIDER = "google"

# Name of the embedding model (provider-specific)
# text-embedding-005
# text-embedding-large-exp-03-07
EMBEDDINGS_NAME = "text-embedding-005"

# Desired output dimensionality for embeddings
# For text-embedding-005, possible values are between 1 and 768
# For text-embedding-large-exp-03-07, possible values are between 1 and 3072
EMBEDDINGS_DIM = 256

# Base cost for edges in the PCST algorithm. This influences the trade-off
# between selecting high-prize nodes/edges and the size of the resulting tree.
PCST_EDGE_BASE_COST = 0.5
# PCST pruning approach
# gw = Goemans-Williamson pruning ; connectfinal = Dijkstra
PCST_PRUNING_APPROACH = "connectfinal"

import functools
import logging
import math
import sys
import time
from multiprocessing import Pool, cpu_count
from concurrent.futures import ThreadPoolExecutor, as_completed

import numpy as np
import pandas as pd
import torch
from datasets import concatenate_datasets, load_dataset
#from datasets.exceptions import DatasetNotFoundError
from s3path import S3Path
from s3torchconnector import S3Checkpoint
from torch_geometric.data.data import Data, DataEdgeAttr, DataTensorAttr
from torch_geometric.data.storage import GlobalStorage
from torch_geometric import utils as tgeo
from tqdm.auto import tqdm


# Create a CosineSimilarity instance to accelerate the compute
cosine_similarity = torch.nn.CosineSimilarity(dim=-1)

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s', force=True)
logger = logging.getLogger(__name__)

# Import PCST library
try:
    from pcst_fast import pcst_fast
    logger.info("Successfully imported pcst_fast.")
except ImportError as e:
    logger.error(f"Failed to import pcst_fast: {e}. Ensure compilation and installation were successful.")

try:
    aws_region = os.getenv("AWS_DEFAULT_REGION", "eu-west-1")
    s3_checkpoint = S3Checkpoint(region=aws_region)
    logger.info(f"S3Checkpoint initialized for region: {aws_region}")
except Exception as e:
    logger.error(f"Failed to initialize S3Checkpoint: {e}")

# Prevents potential issues when saving/loading PyG data with torch.save/load
torch.serialization.add_safe_globals([
    Data,
    DataEdgeAttr,
    DataTensorAttr,
    GlobalStorage
])

if EMBEDDINGS_PROVIDER == "google":
    try:
        import vertexai
        from vertexai.language_models import TextEmbeddingInput, TextEmbeddingModel
        from google.api_core.exceptions import ResourceExhausted, InternalServerError, GoogleAPICallError

        if not GCP_PROJECT or not GCP_LOCATION:
            logger.warning("GCP_PROJECT or GCP_LOCATION are not set. Vertex AI initialization might fail.")
            vertexai.init()
        else:
          vertexai.init(project=GCP_PROJECT, location=GCP_LOCATION)
          logger.info(f"Vertex AI initialized")
    except ImportError:
        logger.error("google-cloud-aiplatform or vertexai not installed. Cannot use Google provider.")
        EMBEDDINGS_PROVIDER = None
    except Exception as e:
        logger.error(f"Failed to initialize Vertex AI: {e}")
        EMBEDDINGS_PROVIDER = None

def load_and_preprocess_dataset(hf_dataset_name):
    """
    Loads the dataset, preprocesses it, and returns the concatenated dataset.

    Args:
        hf_dataset_name (str): The name or path of the dataset on the Hugging
                               Face Hub (e.g., "username/my_dataset").

    Returns:
        Dataset | None: A single `datasets.Dataset` object containing
                        the concatenated data from all valid splits
                        after loading and potential preprocessing.
                        Returns `None` if the original dataset cannot
                        be found, if an unexpected error occurs during
                        loading/processing, or if no valid data splits
                        are available after processing.
    """
    dataset_splits = []
    dataset_id_parts = hf_dataset_name.split("/")
    lettria_dataset_name = f"Lettria/{dataset_id_parts[-1]}" if len(dataset_id_parts) > 1 else f"Lettria/{hf_dataset_name}"
    processed_version_exists = False

    # Try loading the potentially preprocessed version first
    try:
        logger.info(f"Attempting to load potentially preprocessed dataset: {lettria_dataset_name}")
        dataset = load_dataset(lettria_dataset_name)
        logger.info(f"Successfully loaded dataset from {lettria_dataset_name}.")
        # Basic check if it looks processed (assuming preprocessing doesn't remove splits)
        if dataset:
            processed_version_exists = True
            for split_name in dataset.keys():
                 dataset_splits.append(dataset[split_name])

    except (Exception, ValueError) as e: # ValueError for cases like auth issues
        logger.warning(f"Could not load from {lettria_dataset_name}: {e}. Loading original: {hf_dataset_name}")
        processed_version_exists = False
    except Exception as e:
        logger.error(f"An unexpected error occurred while loading {lettria_dataset_name}: {e}")
        processed_version_exists = False

    # If processed version not loaded, load original and process it
    if not processed_version_exists:
        try:
            dataset = load_dataset(hf_dataset_name)
            logger.info(f"Loaded original dataset {hf_dataset_name}. Starting preprocessing...")

            mandatory_keys = {"id", "question", "graph", "answer"}
            num_procs = max(1, cpu_count() // 2) # Use half CPUs to avoid resource exhaustion
            logger.info(f"Using {num_procs} processes for map/filter operations.")

            # --- Preprocessing functions ---
            def _is_valid_example(ex):
                # Basic type and structure checks
                return (
                    isinstance(ex.get("question"), str) and
                    isinstance(ex.get("graph"), list) and
                    isinstance(ex.get("answer"), list) and
                    len(ex["graph"]) > 0 and
                    all(isinstance(tri, list) and len(tri) == 3 for tri in ex["graph"])
                )

            def _clean_graph(ex):
                # Remove triples with invalid entities
                removal_entities = {"", "n/a", None, "None"}
                # Create a new list instead of modifying in place during iteration
                ex["graph"] = [tri for tri in ex["graph"] if tri[0] not in removal_entities and tri[2] not in removal_entities]
                return ex

            def _update_answer(ex):
                # Ensure answers are nodes present in the *cleaned* graph
                nodes_in_graph = set()
                for s, _, o in ex["graph"]:
                    nodes_in_graph.add(s)
                    nodes_in_graph.add(o)

                # Filter answers to keep only those present in the graph nodes
                ex["answer"] = list(set(ex["answer"]) & nodes_in_graph)
                return ex
            # --- End Preprocessing functions ---

            processed_splits = {}
            for split_name in dataset.keys():
                logger.info(f"Processing split: {split_name}")
                split_data = dataset[split_name]

                # Check for mandatory columns
                missing_keys = mandatory_keys - set(split_data.column_names)
                if missing_keys:
                    logger.warning(f"Skipping split '{split_name}'. Missing keys: {missing_keys}")
                    continue

                # 1. Filter invalid examples first
                logger.info(f"  Filtering {split_name} for valid examples...")
                valid_split = split_data.filter(_is_valid_example, num_proc=num_procs)
                logger.info(f"  {len(valid_split)} / {len(split_data)} examples remaining after validation filter.")

                # 2. Clean graph structure
                logger.info(f"  Cleaning graph triples in {split_name}...")
                cleaned_split = valid_split.map(_clean_graph, num_proc=num_procs)

                 # 3. Remove examples with empty graphs after cleaning
                logger.info(f"  Filtering {split_name} for non-empty graphs after cleaning...")
                non_empty_graph_split = cleaned_split.filter(lambda ex: len(ex["graph"]) > 0, num_proc=num_procs)
                logger.info(f"  {len(non_empty_graph_split)} / {len(cleaned_split)} examples remaining after empty graph filter.")

                # 4. Update answers based on the cleaned graph
                logger.info(f"  Updating answers in {split_name}...")
                processed_split = non_empty_graph_split.map(_update_answer, num_proc=num_procs)

                dataset_splits.append(processed_split)
                processed_splits[split_name] = processed_split

            # Push the processed dataset to the Hub
            if dataset_splits:
                logger.info(f"Preprocessing complete. Pushing processed dataset to {lettria_dataset_name}...")
                try:
                    # Need to reconstruct a DatasetDict to push
                    from datasets import DatasetDict

                    processed_dataset_dict = DatasetDict(processed_splits)
                    processed_dataset_dict.push_to_hub(lettria_dataset_name, private=True, max_shard_size="90MB")
                    logger.info(f"Successfully pushed processed dataset to {lettria_dataset_name}.")
                except Exception as e:
                    logger.error(f"Failed to push processed dataset to Hugging Face Hub: {e}")
            else:
                logger.warning("No valid data splits found after preprocessing.")

        except Exception:
            logger.error(f"Original dataset '{hf_dataset_name}' not found on Hugging Face Hub.")
            return None
        except Exception as e:
            logger.error(f"An unexpected error occurred while loading/processing {hf_dataset_name}: {e}")
            return None

    # Concatenate all valid splits
    if dataset_splits:
        logger.info(f"Concatenating {len(dataset_splits)} dataset split(s).")
        concatenated_dataset = concatenate_datasets(dataset_splits)
        logger.info(f"Final dataset size: {len(concatenated_dataset)} examples.")
        return concatenated_dataset
    else:
        logger.error("No dataset splits available to concatenate.")
        return None

concatenated_dataset = load_and_preprocess_dataset(HF_DATASET_NAME)

def calculate_node_prizes(query_embedding, input_graph, top_nodes):
    """
    Calculate node prizes based on similarity to query embedding (Optimized).

    Nodes get prizes based on their rank among the top_nodes most similar nodes.
    Rank 1 (most similar) gets prize top_nodes, Rank top_nodes gets prize 1.
    All other nodes get prize 0.

    Args:
        query_embedding: Query embedding.
        input_graph: Input graph.
        top_nodes: Number of top nodes to assign prizes.

    Returns:
        Calculated node prizes tensor.
    """
    num_nodes = input_graph.num_nodes
    # Ensure prizes are created on the same device as the embeddings/graph features
    device = query_embedding.device if query_embedding is not None else (input_graph.x.device if input_graph.x is not None else 'cpu')
    node_prizes = torch.zeros(num_nodes, dtype=torch.float, device=device)

    if top_nodes > 0 and num_nodes > 0 and input_graph.x is not None:
        # Ensure top_nodes doesn't exceed the number of nodes
        k = min(top_nodes, num_nodes)

        # Calculate cosine similarity
        similarities = cosine_similarity(query_embedding, input_graph.x).float()

        # Get the indices of the top k similar nodes
        # We only need the indices, not the similarity values themselves for prize assignment
        _, top_indices = torch.topk(similarities, k, largest=True)

        # Assign prizes based on rank (top_nodes down to 1)
        # Create prize values tensor: [k, k-1, ..., 1]
        prize_values = torch.arange(k, 0, -1, dtype=torch.float, device=device)
        node_prizes[top_indices] = prize_values
    # Handle cases where input_graph.x might be None
    elif input_graph.x is None:
        print("Warning: input_graph.x is None in calculate_node_prizes. Returning zero prizes.")
        # node_prizes is already zeros

    return node_prizes

def calculate_edge_prizes(query_embedding, input_graph, top_edges):
    """
    Calculate edge prizes based on similarity to query embedding (Optimized).

    Prizes are assigned based on ranked unique similarity values among the top_edges.
    Handles ties by distributing prize budget. Ensures strictly decreasing prize values.

    Args:
        query_embedding: Query embedding.
        input_graph: Input graph.
        top_edges: Number of top ranked unique similarity values to consider.

    Returns:
        Calculated edge prizes tensor.
    """
    num_edges = input_graph.num_edges
    # Ensure prizes are created on the same device as the embeddings/graph features
    device = query_embedding.device if query_embedding is not None else (input_graph.edge_attr.device if input_graph.edge_attr is not None else 'cpu')
    edge_prizes = torch.zeros(num_edges, dtype=torch.float, device=device)

    # Use a small epsilon for strict decrease logic
    EPSILON = 1e-5

    if top_edges > 0 and num_edges > 0 and input_graph.edge_attr is not None:
        # Calculate cosine similarity
        similarities = cosine_similarity(query_embedding, input_graph.edge_attr).float()

        # Get unique similarity values and sort them descending
        unique_similarities = torch.unique(similarities)
        # Check if unique_similarities is empty before sorting
        if unique_similarities.numel() == 0:
            return edge_prizes # Return zeros if no unique similarities

        sorted_unique_similarities = torch.sort(unique_similarities, descending=True).values

        # Determine the actual number of top unique values to use
        k = min(top_edges, len(sorted_unique_similarities))

        if k > 0:
            # Get the top k unique similarity values
            top_k_unique_values = sorted_unique_similarities[:k]
            # Value below which prizes are zero (threshold)
            cutoff_value = top_k_unique_values[-1]

            # Initialize prize tracking
            last_top_edge_value = float('inf') # Initialize high

            # Iterate through the top k unique values to assign prizes
            for rank, value in enumerate(top_k_unique_values):
                # Find all edges matching the current unique similarity value
                matching_indices_mask = (similarities == value)
                count = matching_indices_mask.sum().item() # Ensure count is a scalar

                if count > 0:
                    # Calculate prize value based on rank and count
                    # Prize budget for this rank is (top_edges - rank)
                    # The original code used `top_edges` for budget calculation, preserving that.
                    calculated_prize = (top_edges - rank) / count

                    # Ensure prize is strictly less than the previous level's prize
                    # The min comparison ensures it doesn't exceed the previous prize
                    # The EPSILON subtraction ensures strict decrease for the *next* level
                    prize_value = min(calculated_prize, last_top_edge_value * (1 - EPSILON))

                    # Assign the calculated prize value
                    edge_prizes[matching_indices_mask] = prize_value

                    # Update the ceiling for the next iteration's prize calculation
                    last_top_edge_value = prize_value

            # Explicitly zero out prizes for edges below the cutoff similarity
            # This handles edges that might have had similarities equal to the cutoff
            # but weren't part of the top k unique values processed IF k < top_edges.
            # It also zeros out everything else correctly.
            edge_prizes[similarities < cutoff_value] = 0.0
        # else: k is 0, edge_prizes remains zeros

    # Handle case where edge_attr might be None or num_edges is 0
    elif input_graph.edge_attr is None:
         print("Warning: input_graph.edge_attr is None in calculate_edge_prizes. Returning zero prizes.")
         # edge_prizes is already zeros
    elif num_edges == 0:
         # edge_prizes is already zeros
         pass


    return edge_prizes

def prepare_pcst_data(input_graph, edge_prizes, edge_base_cost):
    """
    Prepare data structures for PCST algorithm (Optimized).

    Uses vectorized operations to separate low-prize and high-prize edges
    and construct the inputs for the PCST solver.

    Args:
        input_graph: Input graph.
        edge_prizes: Calculated edge prizes.
        edge_base_cost: Base cost for edges.

    Returns:
        tuple: (edge_costs_np, pcst_edges_np, virtual_node_prizes_np, virtual_edges_np, virtual_edge_costs_np, node_id_mapping, edge_id_mapping, original_pcst_edge_count)
               Returns numpy arrays suitable for pcst_fast, mappings, and edge count.
               Note: The first element returned in the original function was edge_costs (list), here it's edge_costs_np (numpy array).
               The second element was pcst_edges (list of tuples), here it's pcst_edges_np (numpy array).
               The structure of the returned tuple is slightly different but contains equivalent information.
               The numpy arrays replace the original lists for efficiency.
    """
    num_nodes = input_graph.num_nodes
    num_edges = input_graph.num_edges
    device = edge_prizes.device

    if num_edges == 0:
        # Handle empty graph case
        node_id_mapping = {}
        edge_id_mapping = {}
        edge_costs_np = np.array([], dtype=np.float64)
        pcst_edges_np = np.array([], dtype=np.int64).reshape(0, 2)
        virtual_node_prizes_np = np.array([], dtype=np.float64)
        virtual_edges_np = np.array([], dtype=np.int64).reshape(0, 2)
        virtual_edge_costs_np = np.array([], dtype=np.float64)
        original_pcst_edge_count = 0
        # Match original return signature structure as closely as possible with NumPy arrays
        return (edge_costs_np, pcst_edges_np, virtual_node_prizes_np,
                virtual_edges_np, virtual_edge_costs_np, node_id_mapping,
                edge_id_mapping, original_pcst_edge_count)

    # Identify high-prize edges (create virtual nodes) and low-prize edges (direct PCST edges)
    high_prize_mask = edge_prizes > edge_base_cost
    low_prize_mask = ~high_prize_mask
    pcst_edge_indices_orig = torch.where(low_prize_mask)[0]
    pcst_edges = input_graph.edge_index[:, low_prize_mask].t() # Transpose to get pairs
    edge_costs = torch.full_like(edge_prizes[low_prize_mask], edge_base_cost) - edge_prizes[low_prize_mask]
    original_pcst_edge_count = len(pcst_edges)
    # Create mapping from the new index in pcst_edges to original edge index
    edge_id_mapping = {new_idx: old_idx.item() for new_idx, old_idx in enumerate(pcst_edge_indices_orig)}
    virtual_edge_indices_orig = torch.where(high_prize_mask)[0]
    num_virtual_nodes = len(virtual_edge_indices_orig)

    if num_virtual_nodes > 0:
        # Calculate prizes for virtual nodes
        virtual_node_prizes = edge_prizes[high_prize_mask] - edge_base_cost

        # Assign IDs to virtual nodes (start after original nodes)
        virtual_node_ids = torch.arange(num_nodes, num_nodes + num_virtual_nodes, device=device)
        # Create mapping from virtual node ID to the original edge index it represents
        node_id_mapping = {vn_id.item(): orig_idx.item() for vn_id, orig_idx in zip(virtual_node_ids, virtual_edge_indices_orig)}

        # Get source and destination nodes of the original high-prize edges
        src_nodes_high = input_graph.edge_index[0, high_prize_mask]
        dst_nodes_high = input_graph.edge_index[1, high_prize_mask]

        # Create virtual edges: (src -> virtual_node) and (virtual_node -> dst)
        virtual_edges_part1 = torch.stack([src_nodes_high, virtual_node_ids], dim=0)
        virtual_edges_part2 = torch.stack([virtual_node_ids, dst_nodes_high], dim=0)
        virtual_edges = torch.cat([virtual_edges_part1, virtual_edges_part2], dim=1).t() # Transpose

        # Virtual edges have zero cost
        virtual_edge_costs = torch.zeros(2 * num_virtual_nodes, dtype=torch.float, device=device)

    else:
        # No high-prize edges
        virtual_node_prizes = torch.empty(0, dtype=torch.float, device=device)
        virtual_edges = torch.empty((0, 2), dtype=torch.long, device=device)
        virtual_edge_costs = torch.empty(0, dtype=torch.float, device=device)
        node_id_mapping = {}

    # Ensure costs are non-negative (important for PCST)
    edge_costs = torch.clamp(edge_costs, min=0.0)
    virtual_edge_costs = torch.clamp(virtual_edge_costs, min=0.0) # Should be 0 anyway

    edge_costs_np = edge_costs.cpu().numpy().astype(np.float64)
    pcst_edges_np = pcst_edges.cpu().numpy().astype(np.int64)
    virtual_node_prizes_np = virtual_node_prizes.cpu().numpy().astype(np.float64)
    virtual_edges_np = virtual_edges.cpu().numpy().astype(np.int64)
    virtual_edge_costs_np = virtual_edge_costs.cpu().numpy().astype(np.float64)

    return (edge_costs_np, pcst_edges_np, virtual_node_prizes_np,
            virtual_edges_np, virtual_edge_costs_np, node_id_mapping,
            edge_id_mapping, original_pcst_edge_count)

def run_pcst_algorithm(all_edges, all_prizes, all_edge_costs):
    """
    Run the PCST algorithm (Wrapper).

    Args:
        all_edges: Array of all edges (N, 2).
        all_prizes: Array of all node prizes (including virtual).
        all_edge_costs: Array of all edge costs (including virtual).

    Returns:
        tuple: (selected_vertices, selected_edges indices relative to all_edges)
    """
    # Parameters for pcst_fast (can be made arguments if needed)
    root_node = -1  # No root node specified
    num_clusters = 1  # We want a single connected component
    pruning_method = PCST_PRUNING_APPROACH  # gw = Goemans-Williamson pruning ; connectfinal = Dijkstra
    verbosity_level = 5  # No verbose output

    # Ensure costs are non-negative
    all_edge_costs = np.maximum(all_edge_costs, 0)

    # Check for empty inputs which might cause errors in pcst_fast
    if all_edges.shape[0] == 0:
        # If no edges, PCST result depends only on prizes.
        # Return nodes with positive prizes, and no edges.
        selected_vertices = np.where(all_prizes > 0)[0]
        selected_edges_indices = np.array([], dtype=np.int64)
        return selected_vertices, selected_edges_indices

    if all_prizes.shape[0] == 0:
         # If no nodes (shouldn't really happen if edges exist, but defensively)
         return np.array([], dtype=np.int64), np.array([], dtype=np.int64)


    # Run the PCST algorithm
    start_pcst_time = time.time()
    selected_vertices, selected_edges_indices = pcst_fast(all_edges,
                                                        all_prizes,
                                                        all_edge_costs,
                                                        root_node,
                                                        num_clusters,
                                                        pruning_method,
                                                        verbosity_level)
    end_pcst_time = time.time() - start_pcst_time

    return np.array(selected_vertices), np.array(selected_edges_indices), end_pcst_time

def create_subgraph(input_graph, selected_vertices, selected_pcst_edge_indices, edge_id_mapping, node_id_mapping, original_edge_count, pcst_edges, virtual_edges):
    """
    Create a subgraph from the PCST algorithm results (Optimized).

    Args:
        input_graph: Original input graph.
        selected_vertices: Selected vertices from PCST (including virtual).
        selected_pcst_edge_indices: Indices of selected edges *relative to the concatenated list* fed to PCST.
        edge_id_mapping: Mapping from pcst_edges index to original graph edge index.
        node_id_mapping: Mapping from virtual node ID to original graph edge index.
        original_edge_count: The number of non-virtual edges in the list fed to PCST.
        pcst_edges: The non-virtual edges passed to pcst_fast (Numpy array). Required by optimized logic.
        virtual_edges: The virtual edges passed to pcst_fast (Numpy array). Required by optimized logic.


    Returns:
        tuple: (subgraph, selected_node_original_indices, selected_edge_original_indices)
               The original returned selected_edges was ambiguous; this returns original indices.
    """
    num_original_nodes = input_graph.num_nodes

    # 1. Edges selected directly from the low-prize pool
    direct_mask = selected_pcst_edge_indices < original_edge_count
    selected_direct_indices = selected_pcst_edge_indices[direct_mask]
    # Map these indices (which refer to the `pcst_edges` array) back to original graph indices
    original_direct_edge_indices = [edge_id_mapping[idx] for idx in selected_direct_indices if idx in edge_id_mapping] # Added check if idx exists

    # 2. Edges selected implicitly via virtual nodes
    # A virtual node being selected means the original high-prize edge it represents should be included.
    selected_virtual_nodes = selected_vertices[selected_vertices >= num_original_nodes]
    # Map these virtual node IDs back to the original graph edge indices they represent
    original_virtual_edge_indices = [node_id_mapping[v_node] for v_node in selected_virtual_nodes if v_node in node_id_mapping] # Added check if v_node exists

    # Combine all selected original edge indices, ensuring uniqueness
    selected_edge_original_indices = np.array(list(set(original_direct_edge_indices + original_virtual_edge_indices)), dtype=np.int64)

    # Determine the set of original nodes involved
    # Start with nodes directly selected by PCST (non-virtual ones)
    selected_node_original_indices_direct = selected_vertices[selected_vertices < num_original_nodes]

    # Check if any edges were selected
    if len(selected_edge_original_indices) > 0:
        # Get edge index corresponding to the selected original edges
        subgraph_edge_index_orig = input_graph.edge_index[:, selected_edge_original_indices]
        # Find unique nodes involved in these selected edges
        nodes_from_edges = torch.unique(subgraph_edge_index_orig.flatten()).cpu().numpy()
        # Combine directly selected nodes and nodes from selected edges
        all_selected_nodes_orig_np = np.unique(np.concatenate([selected_node_original_indices_direct, nodes_from_edges]))
    else:
        # No edges selected, only use directly selected nodes
        all_selected_nodes_orig_np = np.unique(selected_node_original_indices_direct)


    # Handle case with no nodes selected at all
    if len(all_selected_nodes_orig_np) == 0:
        # No nodes or edges selected
        subgraph = Data(x=input_graph.x.new_empty((0, input_graph.x.size(1))) if input_graph.x is not None else None,
                        edge_index=input_graph.edge_index.new_empty((2,0)),
                        edge_attr=input_graph.edge_attr.new_empty((0, input_graph.edge_attr.size(1))) if input_graph.edge_attr is not None else None,
                        num_nodes=0)
        return subgraph, np.array([]), np.array([])

    all_selected_nodes_tensor = torch.from_numpy(all_selected_nodes_orig_np).to(input_graph.edge_index.device)

    # Perform subgraph extraction. This returns edge_index relative to the new node set,
    # and edge_attr containing the *original indices* of the edges included in the subgraph.
    subgraph_edge_index, subgraph_edge_attr_indices = tgeo.subgraph(
        subset=all_selected_nodes_tensor,
        edge_index=input_graph.edge_index,
        # Pass original indices as edge_attr to filter later
        edge_attr=torch.arange(input_graph.num_edges, device=input_graph.edge_index.device),
        relabel_nodes=True,
        num_nodes=input_graph.num_nodes
    )

    # Filter the edges returned by `subgraph` to keep only those selected by PCST logic.
    # Create a mask of allowed original edge indices based on PCST results.
    allowed_edge_indices_set = set(selected_edge_original_indices.tolist())

    # Create a boolean mask for edges returned by `subgraph`
    edge_mask_in_subgraph = torch.tensor(
        [idx.item() in allowed_edge_indices_set for idx in subgraph_edge_attr_indices],
        dtype=torch.bool,
        device=subgraph_edge_index.device
    )

    # Apply the mask to get the final edge index and the original indices of these edges
    final_subgraph_edge_index = subgraph_edge_index[:, edge_mask_in_subgraph]
    final_selected_original_edge_indices = subgraph_edge_attr_indices[edge_mask_in_subgraph]

    # Select node features for the subgraph nodes
    subgraph_x = input_graph.x[all_selected_nodes_tensor] if input_graph.x is not None else None

    # Select edge features using the final original indices
    subgraph_edge_attr = input_graph.edge_attr[final_selected_original_edge_indices] if input_graph.edge_attr is not None else None

    # Create the final Data object
    subgraph = Data(x=subgraph_x,
                    edge_index=final_subgraph_edge_index,
                    edge_attr=subgraph_edge_attr,
                    num_nodes=len(all_selected_nodes_orig_np)) # Number of nodes in the subgraph

    return subgraph, all_selected_nodes_orig_np, final_selected_original_edge_indices.cpu().numpy()


def run_retrieval_via_pcst(input_graph, query_embedding, node_text_data, edge_text_data, top_nodes, top_edges, edge_base_cost):
    """
    Perform retrieval using Prize-Collecting Steiner Tree (PCST) algorithm on a graph (Optimized).

    This function processes a graph to select relevant nodes and edges based on a query embedding.
    It uses the PCST algorithm to find a subgraph that maximizes the total prize of selected nodes
    while minimizing the total cost of selected edges.

    Args:
        input_graph (torch_geometric.data.Data): Input graph.
        query_embedding (torch.Tensor): Query embedding.
        node_text_data (pandas.DataFrame): Dataframe containing textual information for nodes.
        edge_text_data (pandas.DataFrame): Dataframe containing textual information for edges.
        top_nodes (int): Number of top nodes to consider based on similarity to query.
        top_edges (int): Number of top edges to consider based on similarity to query.
        edge_base_cost (float): Base cost for edges.

    Returns:
        tuple: (subgraph, description_nodes, description_edges)
            - subgraph (torch_geometric.data.Data): Selected subgraph.
            - description_nodes (pandas.DataFrame): Filtered node text data corresponding to subgraph nodes.
            - description_edges (pandas.DataFrame): Filtered edge text data corresponding to subgraph edges.
            (Note: Original function returned description, which seemed to be node text. This version returns both.)
    """

    # Basic check for graph features needed
    if input_graph.x is None and top_nodes > 0:
         print("Warning: top_nodes > 0 but input_graph.x is None. Node prizes cannot be calculated.")
         # Consider returning early or proceeding without node prizes
    if input_graph.edge_attr is None and top_edges > 0:
         print("Warning: top_edges > 0 but input_graph.edge_attr is None. Edge prizes cannot be calculated.")
         # Consider returning early or proceeding without edge prizes


    # Check if there's any textual data for nodes or edges - Proceed even if text is missing, return None for descriptions.
    # Early exit if graph itself is fundamentally unusable (e.g., no nodes)
    if input_graph.num_nodes == 0:
        print("Warning: Input graph has 0 nodes. Returning empty result.")
        empty_subgraph = Data(x=input_graph.x, edge_index=input_graph.edge_index, edge_attr=input_graph.edge_attr, num_nodes=0)
        empty_nodes_df = pd.DataFrame(columns=node_text_data.columns) if node_text_data is not None else None
        empty_edges_df = pd.DataFrame(columns=edge_text_data.columns) if edge_text_data is not None else None
        return empty_subgraph, empty_nodes_df, empty_edges_df, 0.0

    # 1. Calculate Node and Edge Prizes using the optimized functions (with original names)
    node_prizes_tensor = calculate_node_prizes(query_embedding, input_graph, top_nodes)
    edge_prizes_tensor = calculate_edge_prizes(query_embedding, input_graph, top_edges)

    # 2. Adjust edge_base_cost
    # Ensure cost isn't so high that no edge/virtual node can be selected
    if edge_prizes_tensor.numel() > 0:
         max_edge_prize = edge_prizes_tensor.max().item()
          # Avoid division by zero or issues if max_edge_prize is non-positive
         if max_edge_prize > 1e-9: # Use a small threshold
             # Adjust cost slightly below the max prize to allow selection
             adjusted_edge_base_cost = min(edge_base_cost, max_edge_prize * 0.995)
         else:
             adjusted_edge_base_cost = edge_base_cost # Keep original if max prize is effectively zero or negative
    else:
         # No edges or prizes calculated, keep original cost
         adjusted_edge_base_cost = edge_base_cost


    # 3. Prepare Data for PCST
    # The function returns numpy arrays directly and the original_pcst_edge_count
    (edge_costs_np, pcst_edges_np, virtual_node_prizes_np,
     virtual_edges_np, virtual_edge_costs_np, node_id_mapping,
     edge_id_mapping, original_pcst_edge_count) = prepare_pcst_data(
         input_graph, edge_prizes_tensor, adjusted_edge_base_cost
     )

    # 4. Combine Prizes, Edges, and Costs for PCST
    # Node prizes (original nodes)
    node_prizes_np = node_prizes_tensor.cpu().numpy().astype(np.float64)

    # Combine original node prizes and virtual node prizes
    # Ensure node_prizes_np has size input_graph.num_nodes
    if len(node_prizes_np) != input_graph.num_nodes:
         # Handle potential mismatch, e.g., pad with zeros if necessary
         print(f"Warning: Node prize dimension mismatch ({len(node_prizes_np)} vs {input_graph.num_nodes}). Padding with zeros.")
         padded_prizes = np.zeros(input_graph.num_nodes, dtype=np.float64)
         limit = min(len(node_prizes_np), input_graph.num_nodes)
         padded_prizes[:limit] = node_prizes_np[:limit]
         node_prizes_np = padded_prizes

    # Ensure virtual_node_prizes_np is numpy array even if empty
    if not isinstance(virtual_node_prizes_np, np.ndarray):
        virtual_node_prizes_np = np.array(virtual_node_prizes_np, dtype=np.float64)


    all_prizes_np = np.concatenate([node_prizes_np, virtual_node_prizes_np])

    # Combine direct PCST edges and virtual edges
    # Ensure edge arrays are 2D even if empty
    if pcst_edges_np.ndim == 1: pcst_edges_np = pcst_edges_np.reshape(-1, 2)
    if virtual_edges_np.ndim == 1: virtual_edges_np = virtual_edges_np.reshape(-1, 2)

    if virtual_edges_np.shape[0] > 0:
        all_edges_np = np.vstack([pcst_edges_np, virtual_edges_np])
        all_edge_costs_np = np.concatenate([edge_costs_np, virtual_edge_costs_np])
    else:
        all_edges_np = pcst_edges_np
        all_edge_costs_np = edge_costs_np

    # Handle case where there are no edges at all for PCST
    if all_edges_np.shape[0] == 0 and all_prizes_np.sum() <= 0: # Check if sum is non-positive
        print("Warning: No edges and no positive prizes for PCST. Returning empty subgraph.")
        empty_subgraph = Data(x=input_graph.x.new_empty((0, input_graph.x.size(1))) if input_graph.x is not None else None,
                              edge_index=input_graph.edge_index.new_empty((2, 0)),
                              edge_attr=input_graph.edge_attr.new_empty((0, input_graph.edge_attr.size(1))) if input_graph.edge_attr is not None else None,
                              num_nodes=0)
        empty_nodes_df = pd.DataFrame(columns=node_text_data.columns) if node_text_data is not None else None
        empty_edges_df = pd.DataFrame(columns=edge_text_data.columns) if edge_text_data is not None else None
        return empty_subgraph, empty_nodes_df, empty_edges_df, 0.0


    # 5. Run PCST Algorithm using the wrapper
    selected_vertices, selected_pcst_edge_indices, end_pcst_time = run_pcst_algorithm(
        all_edges_np, all_prizes_np, all_edge_costs_np
    )

    # 6. Create Subgraph from Results using the optimized function (with original name)
    # Pass the necessary arguments, including the components fed to PCST
    subgraph, selected_node_idxs, selected_edge_idxs = create_subgraph(
        input_graph,
        selected_vertices,
        selected_pcst_edge_indices,
        edge_id_mapping,
        node_id_mapping,
        original_pcst_edge_count, # Use the count returned by prepare_pcst_data
        pcst_edges_np, # Pass the actual edges fed to PCST
        virtual_edges_np # Pass the actual virtual edges fed to PCST
    )

    # 7. Generate Textual Description
    # Filter pandas DataFrames based on the selected original indices
    description_nodes = node_text_data.iloc[selected_node_idxs] if node_text_data is not None and len(selected_node_idxs) > 0 else \
                       (pd.DataFrame(columns=node_text_data.columns) if node_text_data is not None else None) # Handle empty selection or missing df

    description_edges = edge_text_data.iloc[selected_edge_idxs] if edge_text_data is not None and len(selected_edge_idxs) > 0 else \
                       (pd.DataFrame(columns=edge_text_data.columns) if edge_text_data is not None else None) # Handle empty selection or missing df


    return subgraph, description_nodes, description_edges, end_pcst_time

def find_empty_answer_examples(dataset):
    """
    Returns indices of examples with empty answer lists and the count of non-empty ones.

    Args:
        dataset (Dataset): A Hugging Face Dataset.

    Returns:
        tuple: A tuple containing:
            - list[int]: A list of integer indices corresponding to the examples
                         in the input `dataset` that have an empty or missing 'answer'.
            - int: The total count of examples in the `dataset` that have a
                   non-empty 'answer'. Returns ([], 0) if the input dataset is None.
    """
    if dataset is None:
        return [], 0

    empty_answer_indices = [
        rank for rank, ex in enumerate(dataset)
        if not ex.get("answer") # Check if 'answer' key exists and is not empty
    ]
    total_examples = len(dataset)
    non_empty_answer_count = total_examples - len(empty_answer_indices)
    logger.info(f"Found {len(empty_answer_indices)} examples with empty answers.")
    logger.info(f"{non_empty_answer_count} examples have ground truth answers.")
    return empty_answer_indices, non_empty_answer_count

if concatenated_dataset:
    examples_with_empty_answers_indices, total_examples_for_recall = find_empty_answer_examples(concatenated_dataset)
else:
    examples_with_empty_answers_indices, total_examples_for_recall = [], 0

s3_examples_path = S3_EXAMPLES_PATH
example_id = "WebQTrn-141"
base_embedding_suffix = f'{EMBEDDINGS_NAME}_{EMBEDDINGS_DIM}'
graph_embedding_path = f'{s3_examples_path}/{example_id}/{base_embedding_suffix}_graph_embedding.pt'
question_embedding_path = f'{s3_examples_path}/{example_id}/{base_embedding_suffix}_question_embedding.pt'
nodes_csv_path = f"{s3_examples_path}/{example_id}/nodes.csv"
edges_csv_path = f"{s3_examples_path}/{example_id}/edges.csv"

try:
    with s3_checkpoint.reader(graph_embedding_path) as reader:
        graph = torch.load(reader, weights_only=True)
    with s3_checkpoint.reader(question_embedding_path) as reader:
        q_emb = torch.load(reader, weights_only=True)
    nodes_df_orig = pd.read_csv(nodes_csv_path)
    edges_df_orig = pd.read_csv(edges_csv_path)
except Exception as e:
    logger.error(f"Rank {rank} (ID: {example_id}): Failed to load input data: {e}")
    # Return structure indicating failure/no processing for aggregation

k_nodes = math.ceil((1 / 100.0) * len(nodes_df_orig)) if len(nodes_df_orig) > 100 else 1 * 10 # Use fixed count for small graphs
k_edges = math.ceil((1 / 100.0) * len(edges_df_orig)) if len(edges_df_orig) > 100 else 1 * 10 # Use fixed count for small graphs
subgraph, subgraph_nodes_df, subgraph_edges_df, pcst_time = run_retrieval_via_pcst(
    graph,
    q_emb,
    nodes_df_orig,
    edges_df_orig,
    k_nodes,
    k_edges,
    PCST_EDGE_BASE_COST
)