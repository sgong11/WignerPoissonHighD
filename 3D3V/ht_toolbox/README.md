# HT Toolbox for Hierarchical Tucker Decomposition (HTD)

This is a lightweight MATLAB toolbox for building and manipulating tensors in the Hierarchical Tucker (HT) format. The toolbox is specifically designed for clarity, extensibility, and efficiency—suitable for both research and teaching. It provides a functional (non-class-based) interface, customized for use in **Hierarchical Tucker Adaptive Cross Approximation (HTACA)** algorithms.

---

## 📌 Highlights

- ✅ **Inspired by the `htucker` toolbox** developed by *Daniel Kressner* and *Christine Tobler*. We gratefully acknowledge their foundational contribution.  
  Our implementation retains the same fundamental tree structure—namely `children` and `dim2ind`—but replaces the class-based design with a fully functional one. We have reimplemented all core routines from scratch, and extended the tree structure with additional attributes (e.g., nonleaf, depth, freedom, orders) to better suit our own algorithmic needs.

- 🔧 **Developed for research-grade prototyping of HTACA algorithms.**  
  The current MATLAB implementation focuses on algorithm development and reproducible experimentation for Hierarchical Tucker Adaptive Cross Approximation (HTACA).  
  It is intended as a reference platform for verifying mathematical ideas, exploring truncation strategies, and testing recursive compression pipelines.

  For production-level performance or GPU acceleration, future implementations in C++/CUDA are being planned.

- 📂 **Self-contained and class-free.**  
  All functions are standalone, with no external dependencies or object-oriented overhead—ideal for customization and embedding.

- 📚 **Designed for extensibility and clarity—ideal for research and teaching.**

---

## 🗂️ Directory Structure

The toolbox is organized as follows:

### 🔹 `+ht/`: Core HTD Functional Routines

This MATLAB package contains all core functions used for constructing, transforming, evaluating, and manipulating Hierarchical Tucker (HT) tensors.

| Function Group              | Representative Functions                                                                 | Description |
|----------------------------|--------------------------------------------------------------------------------------------|-------------|
| **Tree Construction**      | `make_tree`, `validate_tree`, `tree_depth`, `rank`                                        | Construct and validate HT trees with metadata such as depth, orders, and rank |
| **Initializers**           | `zeros`, `ones`, `rand`, `rand_complex`                                                   | Create HT tensors with preset values (zero, one, random real/complex) |
| **Basic Algebra**          | `add`, `add_scalar`, `multiply_scalar`                                                    | Perform HTD-level addition, scalar shifts, and scaling |
| **Norm & Inner Product**   | `norm`, `inner_product`                                                                   | Compute Frobenius norm and pairwise inner product |
| **Orthogonalization**      | `orthogonalize`                                                                           | Left-orthogonalize all nodes of an HTD |
| **Low-Rank Truncation**    | `truncate`, `gramians_orthog`                                                             | Compute Gramians and truncate using Frobenius norm threshold |
| **Tensor Transformations** | `mode_matrix_product`, `mode_apply_function`, `ttm`                                       | Matrix-tensor mode products and per-mode transformations |
| **Evaluation Operators**   | `evaluate_index`, `evaluate_slice`, `evaluate_fiber`, `evaluate_index_sum`, `evaluate_fiber_sum` | Pointwise, fiber, or slice evaluation and summation |
| **Rank Queries**           | `get_default_min_rank`, `get_default_max_rank`, `choose_rank_fro`                         | Automatic rank suggestion and tolerance-based rank selection |
| **Subtree Tools**          | `get_left_subtree`, `get_right_subtree`, `get_subtree_pool`                               | Retrieve subtree structure and subtree indices |
| **HTACA Construction**     | `HTACA`                                                                                   | Adaptive sampling-based HT construction from function handle input |
| **Utility Functions**      | `squeeze`, `determine_sample_size`                                                        | Remove singleton modes, determine sample count based on tree depth |

---

### 📁 `examples/`: Usage Demonstration Scripts

This folder contains lightweight example scripts for each functionality group:

| Script                             | Demonstration Purpose |
|------------------------------------|------------------------|
| `example_build_tree.m`            | Constructing trees and rank initialization |
| `example_construct_htd.m`         | Creating random HT tensors |
| `example_add_and_scalar.m`        | Addition and scalar operations |
| `example_norm_and_inner.m`        | Norm and inner product computation |
| `example_mode_transform.m`        | Applying mode-wise functions or projections |
| `example_truncate_and_gramian.m`  | Orthogonalization, Gramians, and truncation |
| `example_ttm.m`                   | Mode-n matrix multiplication using `ht.ttm` for general multidimensional arrays|
| `example_evaluate_samples.m`      | Pointwise, slice, and fiber evaluation |
| `example_squeeze.m`               | Reducing dimensions using `ht.squeeze` |
| `example_rank_query.m`            | Auto rank computation and query utilities |
| `example_get_subtree_pool.m`      | Extracting subtree pools for recursive compression |
| `example_htaca.m`                 | HTACA construction of compressible high-dimensional tensors |


---

## 🔧 Requirements

- **MATLAB R2023b** or newer is recommended.  
  The toolbox may be compatible with earlier versions, but it has been developed and tested only on MATLAB **2023b**.

---

## 📚 Suggested References

For users who are new to the **Hierarchical Tucker Decomposition (HTD)**, we recommend the following foundational papers:

- **Introductory**:  
  D. Kressner and C. Tobler,  
  *htucker—A MATLAB toolbox for tensors in hierarchical Tucker format*,  
  Mathicse, EPF Lausanne, 2012.  
  [Link: https://infoscience.epfl.ch/record/180138](https://infoscience.epfl.ch/record/180138)  

- **Theoretical Foundation**:  
  L. Grasedyck,  
  *Hierarchical singular value decomposition of tensors*,  
  SIAM J. Matrix Anal. Appl., 31(4), 2029–2054, 2010.  
  [Link: https://doi.org/10.1137/090752286](https://doi.org/10.1137/090752286)

These provide excellent context on the HT structure, SVD-based compression, and basic algorithms.

For readers interested in our **Hierarchical Tucker Adaptive Cross Approximation (HTACA)** framework, we encourage you to follow our upcoming article:

> **A Semi-Lagrangian Adaptive Rank (SLAR) Method II: High-dimensional Vlasov Dynamics**  
> *Zheng, Nanyi*, et al. (Preprint, 2025).  
> [arXiv / DOI link to be added upon release]

This paper introduces HTACA in the context of high-dimensional kinetic PDEs and details algorithmic insights behind this toolbox.

---

## 📜 License and Attribution
This toolbox is open for academic and non-commercial use.

If you find this toolbox useful in your work, we kindly ask you to acknowledge or cite it appropriately, depending on your use case:

🔬 For theoretical contributions (e.g., HTACA algorithm, adaptive rank selection, PDE solver design, etc.), please cite our research paper:

Zheng, Nanyi, et al.
A Semi-Lagrangian Adaptive Rank (SLAR) Method II: High-dimensional Vlasov Dynamics
(Preprint, 2025). [DOI / arXiv link to be added here]

🧰 For code reuse or adaptation (e.g., using ht.build_subtree, htaca, choose_rank_fro, etc.), please cite the software:

Zheng, Nanyi.
Lightweight Hierarchical Tucker Toolbox for HTACA.
GitHub, 2025. [Repository link / Zenodo DOI to be added here]

Please include both citations if your work involves both theoretical understanding and practical reuse of the implementation.
---

## ✏️ Contact and Contribution

For suggestions, bug reports, or collaborative extensions, please feel free to reach out or submit pull requests (if used in a repository).  
We welcome constructive feedback to make the toolbox more robust and usable.

---

