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

#' Detect functional interaction structure with interaction count approach.
#'
#' This approach is introduced in Section 4.1.
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


#' Computes the BCCM ensemble with blocks corresponding to roles.
#'
#' ...
#'
#' @param network ...
#' @param roles ...
#' @param ... ...
#' @returns ...
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


#' Detect functional interaction structure with BCCM interaction propensities.
#'
#' ...
#'
#' @param role_ensemble ...
#' @returns ...
#' @export
role_network_bccm <- function(role_ensemble) {
  #   mode = "directed" because blockOmegas are not symmetric in general even
  #   for undirected interaction networks. The propensity Omega_ij of role i
  #   depends on i's interactions with all other roles and not only role j.
  return(graph_from_adjacency_matrix(role_ensemble$blockOmega,
                                     mode = "directed", weighted = TRUE))
}


#' Detect functional interaction structure with normalized BCCM approach.
#'
#' ...
#'
#' @param role_ensemble ...
#' @param category ...
#' @returns ...
#' @export
role_network_bccm_norm <- function(role_ensemble, category = NULL) {

  # Normalize omega over out-degrees
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


#' Plots a role network.
#'
#' ...
#'
#' @param role_network ...
#' @param role_names ...
#' @param role_colors ...
#' @param edge_weights_range_from ...
#' @param node_sizes ...
#' @param node_sizes_range_from ...
#' @param node_sizes_range_to ...
#' @param threshold ...
#' @param show_weights ...
#' @param ... ...
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
#' ...
#'
#' @param measurements ...
#' @returns ...
.simplify <- function(measurements) {
  if (is.atomic(measurements[[1]]) && all(sapply(measurements, length) == 1))
    return(unlist(measurements))
  else
    return(measurements)
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
