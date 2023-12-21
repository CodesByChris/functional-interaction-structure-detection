# Functions for interaction structure detection in teams.

library(assertthat)
library(ghypernet)
library(igraph)
library(magrittr)
library(pbmcapply)
library(purrr)
library(scales)
library(tibble)
library(tidygraph)



################################################################################
### Step 1: Interaction Structure Detection
################################################################################

#' Detect functional interaction structure as interaction counts.
#'
#' This approach is introduced in Section 4.1. An example is Figure 2.
#'
#' The functional interaction structure is essentially a fully-connected igraph
#' graph whose nodes are the roles and the edges are weighted by the number of
#' interactions between the individuals with the corresponding roles in the
#' original network.
#'
#' @param network igraph graph from which to deduce the interaction preferences
#'   between the roles. The roles can be supplied either as a vertex property
#'   "roles" or via the second argument roles.
#' @param roles (Optional) Vector of length `vcount(network)` whose entry i is a
#'   string specifying the role of vertex i.
#' @returns Functional interaction structure as an igraph graph whose nodes are
#'   the roles and edge weights are the interaction counts.
#' @export
role_network_counts <- function(network, roles = V(network)$role) {

  # Get roles
  role_names <- unique(roles)

  # Count-based functional networks
  network %<>%
    as_tbl_graph() %>%
    activate(edges) %>%
    mutate(from_role = .N()$role[from], to_role = .N()$role[to])

  # Count the number of interactions between each pair of roles
  n_roles <- length(role_names)
  interaction_counts <- matrix(NA, n_roles, n_roles) %>%
    set_rownames(role_names) %>%
    set_colnames(role_names)

  for (role_from in role_names)
    for (role_to in role_names) {
      interaction_counts[role_from, role_to] <- network %>%
        activate(edges) %>%
        tidygraph::filter(from_role == role_from & to_role == role_to) %>%
        ecount()
    }

  # Construct role network
  #  `directed = is_directed(network)` because the number of edges between two
  #  roles is symmetric (i.e. undirected) for undirected networks.
  mode <- if (is_directed(network)) "directed" else "undirected"
  return(graph_from_adjacency_matrix(interaction_counts, mode = mode,
                                     weighted = TRUE))
}


#' Fit BCCM ensemble with blocks corresponding to roles.
#'
#' This function fits a BCCM ensemble as explained in Section 4.4. The fitted
#' ensemble serves as input (i) for functional interaction structure detection
#' as explained in Section 4.5 (see `role_network_bccm` and
#' `role_network_bccm_norm`) and (ii) for interaction structure optimization as
#' explained in Section 7 (see `benchmark_scenario`).
#'
#' @param network igraph graph from which to deduce the role interaction
#'   preferences.
#' @param roles (Optional) Vector of role names. The entry `roles[i]` is the
#'   name of node i's role. If the argument is not given, the roles are taken
#'   from a vertex property `roles` in `network`. If the argument is given, it
#'   takes precedence over the vertex property.
#' @param ... Further args that are forwarded unchanged to `ghypernet::bccm`.
#'   Useful args might be `directed` and `selfloops`.
#' @returns BCCM ensemble fitted to the given network and roles.
#' @export
role_ensemble <- function(network, roles = V(network)$role, ...) {

  # Handle special cases
  if (length(unique(roles)) == 0) {
    warning("Empty network given, which can have no roles.")
    return(NA)
  }
  if (length(unique(roles)) == 1) {
    warning("All nodes have the same role: \"", unique(roles), "\"")
    return(NA)
  }

  # Fit BCCM
  network %>%
    as_adj(sparse = FALSE) %>%
    bccm(labels = roles, directedBlocks = is_directed(network), ...) %>%
    return()
}


#' Detect functional interaction structure as BCCM interaction propensities.
#'
#' This approach is explained in Section 4.4. An example is Figure 3(b).
#'
#' @param role_ensemble: BCCM ensemble from which to compute the functional
#'   interaction structure.
#' @returns Functional interaction structure as an igraph graph whose nodes are
#'   the roles and edge weights are the propensities between the roles.
#' @export
role_network_bccm <- function(role_ensemble) {
  #   mode = "directed" because blockOmegas are not symmetric in general even
  #   for undirected interaction networks. The propensity Omega_ij of role i
  #   depends on i's interactions with all other roles and not only with role j.
  return(graph_from_adjacency_matrix(role_ensemble$blockOmega,
                                     mode = "directed", weighted = TRUE))
}


#' Detect functional interaction structure as normalized BCCM propensities.
#'
#' This approach normalizes the edge weights from `role_network_bccm` across the
#' out-edges of each role to obtain interaction probabilities for roles. The
#' approach is explained in Section 4.5. Examples are:
#' 1. Figure 3(c) for category = NULL.
#' 2. Figure 5(b) for positive and negative preferences.
#'
#' @param role_ensemble: BCCM ensemble from which to compute the functional
#'   interaction structure.
#' @param category: If NULL, the normalized propensities are the edge weights.
#'   If "positive" or "negative", the positive or negative interaction
#'   preferences as explained in Section 5 are the edge weights.
#' @returns Functional interaction structure as an igraph graph whose nodes are
#'   the roles and edge weights are the normalized propensities between roles.
#' @export
role_network_bccm_norm <- function(role_ensemble, category = NULL) {

  # Normalize omega over out-edges
  normalize_row <- function(row) if (any(row != 0)) row / sum(row) else row
  normalized_omega <- role_ensemble$blockOmega %>%
    apply(MARGIN = 1, normalize_row) %>%
    t()  # transpose because apply stores results as columns

  roles <- rownames(role_ensemble$blockOmega)
  rownames(normalized_omega) <- roles
  colnames(normalized_omega) <- roles

  # Return network
  if (is.null(category))
    return(graph_from_adjacency_matrix(normalized_omega, mode = "directed",
                                       weighted = TRUE))

  omega_at_random <- 1 / length(roles)
  if (category == "positive")
    # e.g. 2 means twice as many as expected at random
    categorized_omega <- ifelse(normalized_omega > omega_at_random & normalized_omega != 0,
                                normalized_omega / omega_at_random, 0)
  else
    # e.g. 2 means half as many as expected at random
    categorized_omega <- ifelse(normalized_omega < omega_at_random & normalized_omega != 0,
                                omega_at_random / normalized_omega, 0)
  return(graph_from_adjacency_matrix(categorized_omega, mode = "directed",
                                     weighted = TRUE))
}


#' Plots a functional interaction structure given as a role network.
#'
#' It plots networks generated by `role_network_counts`, `role_network_bccm`,
#' and `role_network_bccm_norm`.
#'
#' @param role_network Weighted network representing the functional interaction
#'   structure to be plotted.
#' @param role_names (Optional) Vector containing the (unique) names of all
#'   roles. This argument can be useful, for example, when calling
#'   `plot_role_network` for a series of networks with different roles but whose
#'   plots should display all roles to preserve the layout across plots.
#' @param role_colors (Optional) Vector defining the colors of the roles.
#' @param edge_weights_range_from (Optional) Range of the smallest and largest
#'   possible edge weight. Setting this argument can, for example, ensure
#'   comparable edge widths across separate plots. The default is the minimum
#'   and maximum edge weight in `role_network`.
#' @param node_sizes (Optional) If `"Outdegree"`, a node's size is proportional
#'   to its outdegree. If a numeric vector with the same length as `role_names`,
#'   the i'th entry specifies the node size of the i'th node. If a number, all
#'   nodes have this size. The default is `"Outdegree"`.
#' @param node_sizes_range_from (Optional) Like `edge_weights_range_from` but
#'   for node sizes. The default is the minimum and maximum value in
#'   `node_sizes`.
#' @param node_sizes_range_to (Optional) Interval to which the node sizes are
#'   rescaled before plotting. The default is the interval $[1, 30]$.
#' @param threshold (Optional) If given, all edges with weights below this
#'   threshold are omitted in the plot. The default is to plot all edges.
#' @param show_weights (Optional) If TRUE, edges are labelled with their
#'   weights. The default is `FALSE`.
#' @param ... Further args that are forwarded unchanged to `igraph::plot`.
#' @export
plot_role_network <- function(role_network,
                              role_names = V(role_network)$name,
                              role_colors = NULL,
                              edge_weights_range_from =
                                range(E(role_network)$weight),
                              node_sizes = "Outdegree",
                              node_sizes_range_from = NULL,
                              node_sizes_range_to = c(1, 30),
                              threshold = 0,
                              show_weights = FALSE,
                              ...) {

  # Adjust args
  if (is.string(node_sizes) && node_sizes == "Outdegree")
    node_sizes <- degree(role_network, mode = "out")
  if (is.null(node_sizes_range_from))
    node_sizes_range_from <- range(min(node_sizes), max(node_sizes))
  role_network %<>% as_tbl_graph()

  # Add inactive nodes
  role_network %<>%
    activate(nodes) %>%
    mutate(active = TRUE)

  inactive_roles <- setdiff(role_names, V(role_network)$name) %>%
    tibble(name = ., active = FALSE)
  role_network %<>% bind_nodes(inactive_roles)

  # Set labels
  role_network %<>% mutate(label = name)
  if (show_weights)
    role_network %<>%
      activate(edges) %>%
      mutate(label = as.character(round(weight, 2)))

  # Set node sizes
  V(role_network)$size <- node_sizes %>%
    rescale(from = node_sizes_range_from, to = node_sizes_range_to)

  # Filter and set edge sizes
  role_network %<>%
    activate(edges) %>%
    filter(weight > threshold) %>%
    mutate(width = rescale(weight, from = edge_weights_range_from, to = c(1, 15)))

  # Colors
  if (!is.null(role_colors)) {
    # Node color: role's color
    role_network %<>%
      activate(nodes) %>%
      mutate(color = ifelse(active,
                            role_colors[match(name, names(role_colors))],
                            "gray85"))

    # Edge color: "from"-node's color
    role_network %<>%
      activate(edges) %>%
      mutate(color = role_colors[.N()[from, ]$name])
  }

  # Plot
  expanded_network <- make_full_graph(length(role_names))
  V(expanded_network)$name <- role_names
  layout <- layout_in_circle(expanded_network)

  plot(role_network, layout = layout, vertex.label = NA, edge.label.cex = 1.75,
       edge.label.color = "black", edge.curved = 0.15, ...)
}



################################################################################
### Step 2: Interaction Structure Optimization
################################################################################

#' Simplifies its argument.
#'
#' Simplifies a list into a vector if its entries are of an atomic type (i.e.
#' numeric, boolean, etc.).
#'
#' @param values Object to simplify.
#' @returns Simplified object.
.simplify <- function(values) {
  if (is.atomic(values[[1]]) && all(sapply(values, length) == 1))
    return(unlist(values))
  else
    return(values)
}


#' Bootstraps measure by sampling n_samples networks from the role_ensemble.
#'
#' @param role_ensemble ...
#' @param measure ...
#' @param n_samples ...
#' @param n_cores ...
#' @returns ...
#' @export
benchmark_scenario <- function(role_ensemble, measure, n_samples = 10000,
                               n_cores = 1) {

  # Define helper function run in parallel
  boot <- function(adj, measure, directed) {
    mode <- if (directed) "directed" else "undirected"
    return(measure(graph_from_adjacency_matrix(adj, mode = mode)))
  }

  # Bootstrap measurements
  adjs <- rghype(n_samples, role_ensemble)
  if (n_cores > 1)
    values <- pbmclapply(adjs, possibly(boot, NA), mc.cores = n_cores,
                         measure = measure, directed = role_ensemble$directed)
  else
    values <- lapply(adjs, possibly(boot, NA), measure = measure,
                     directed = role_ensemble$directed)
  return(.simplify(values))
}
