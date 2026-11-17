library(dplyr)
library(plotly)
library(fields)
library(rmcfs)
library(devtools)
library(R.ROSETTA)
library(tidyverse)
library(csVisuNet)
library(RCy3)

ros <- readRDS("rosetta_results.rds")
rules <- readRDS("rules_project4.rds")

#Visunet visualization
cytoscapePing()
vis <- visunetcyto(rules[1:100, ])



make_topo_plot <- function(
    network,
    network_name,
    grid_size = 100,
    sigma = 0.2,
    edge_quantile = 0.80
) {

  # 1. Extract existing VisuNet node metrics

  node_points <- data.frame(
    x = network$nodes$NodeConnection,
    y = network$nodes$meanDecisionCoverage,
    z = network$nodes$meanSupp,
    node = network$nodes$label,
    node_id = network$nodes$id
  )


  # 2. Remove invalid rows

  node_points <- node_points[
    is.finite(node_points$x) &
    is.finite(node_points$y) &
    is.finite(node_points$z),
  ]

  if (nrow(node_points) < 3) {
    stop(paste("Not enough nodes in", network_name))
  }


  # 3. Identify duplicate X/Y coordinates

  duplicate_groups <- split(
    node_points,
    interaction(node_points$x, node_points$y, drop = TRUE)
  )

  duplicate_check <- do.call(
    rbind,
    lapply(duplicate_groups, function(group) {
      data.frame(
        x = group$x[1],
        y = group$y[1],
        number_of_nodes = nrow(group),
        number_of_different_z = length(unique(group$z)),
        z_values = paste(group$z, collapse = ", ")
      )
    })
  )

  cat("\nDuplicate coordinates in", network_name, ":\n")

  print(
    duplicate_check[
      duplicate_check$number_of_nodes > 1,
    ]
  )


  # 4. Create unique interpolation coordinates

  surface_points <- aggregate(
    z ~ x + y,
    data = node_points,
    FUN = mean
  )

  surface_points <- surface_points[
    is.finite(surface_points$x) &
    is.finite(surface_points$y) &
    is.finite(surface_points$z),
  ]

  if (nrow(surface_points) < 3) {
    stop(paste("Not enough unique coordinates in", network_name))
  }

  if (length(unique(surface_points$x)) < 2 ||
      length(unique(surface_points$y)) < 2) {
    stop(paste("Coordinates do not span two dimensions in", network_name))
  }


  # 5. Normalise X and Y

  surface_points$x_norm <- (
    surface_points$x - min(surface_points$x)
  ) / diff(range(surface_points$x))

  surface_points$y_norm <- (
    surface_points$y - min(surface_points$y)
  ) / diff(range(surface_points$y))

  node_coords <- as.matrix(
    surface_points[, c("x_norm", "y_norm")]
  )


  # 6. Calculate Gaussian RBF kernel

  distance_matrix <- as.matrix(
    dist(node_coords)
  )

  kernel_matrix <- exp(
    -(distance_matrix^2) / (2 * sigma^2)
  )


  # 7. Check matrix dimensions and conditioning

  cat("\nInterpolation points:", nrow(surface_points), "\n")
  cat("Kernel matrix dimensions:", dim(kernel_matrix), "\n")
  cat("Number of Z values:", length(surface_points$z), "\n")

  if (nrow(kernel_matrix) != length(surface_points$z)) {
    stop("RBF matrix and Z values have incompatible dimensions")
  }

  condition_number <- kappa(kernel_matrix)

  cat("Condition number:", condition_number, "\n")

  if (!is.finite(condition_number) ||
      condition_number > 1e12) {
    stop(paste(
      "RBF matrix is ill-conditioned in",
      network_name,
      "- condition number:",
      condition_number
    ))
  }


  # 8. Calculate RBF weights

  rbf_weights <- solve(
    kernel_matrix,
    matrix(surface_points$z, ncol = 1)
  )


  # 9. Create grid including original coordinates

  x_grid <- sort(unique(c(
    seq(
      min(surface_points$x),
      max(surface_points$x),
      length.out = grid_size
    ),
    surface_points$x
  )))

  y_grid <- sort(unique(c(
    seq(
      min(surface_points$y),
      max(surface_points$y),
      length.out = grid_size
    ),
    surface_points$y
  )))

  grid <- expand.grid(
    x = x_grid,
    y = y_grid
  )


  # 10. Normalise grid coordinates

  grid$x_norm <- (
    grid$x - min(surface_points$x)
  ) / diff(range(surface_points$x))

  grid$y_norm <- (
    grid$y - min(surface_points$y)
  ) / diff(range(surface_points$y))

  grid_coords <- as.matrix(
    grid[, c("x_norm", "y_norm")]
  )


  # 11. Calculate distances from grid to original nodes

  grid_distances <- fields::rdist(
    grid_coords,
    node_coords
  )


  # 12. Interpolate Z values

  grid_kernel <- exp(
    -(grid_distances^2) / (2 * sigma^2)
  )

  grid$z <- as.vector(
    grid_kernel %*% rbf_weights
  )


  # 13. Convert predictions into surface matrix

  z_matrix <- t(matrix(
    grid$z,
    nrow = length(x_grid),
    ncol = length(y_grid)
  ))


  # 14. Preserve exact original surface heights

  for (i in seq_len(nrow(surface_points))) {

    x_index <- match(surface_points$x[i], x_grid)
    y_index <- match(surface_points$y[i], y_grid)

    z_matrix[y_index, x_index] <- surface_points$z[i]
  }


  # 15. Filter network edges by connection strength

  network_edges <- network$edges

  if (is.null(network_edges) || nrow(network_edges) == 0) {
    stop(paste("No network edges found in", network_name))
  }

  if (!"conn" %in% names(network_edges)) {
    stop("Network edges do not contain a conn column")
  }

  edge_threshold <- quantile(
    network_edges$conn,
    edge_quantile,
    na.rm = TRUE
  )

  filtered_edges <- network_edges[
    !is.na(network_edges$conn) &
    network_edges$conn >= edge_threshold,
  ]

  cat("\nEdge filtering for", network_name, ":\n")
  cat("Original edges:", nrow(network_edges), "\n")
  cat("Selected edges:", nrow(filtered_edges), "\n")
  cat("Connection threshold:", edge_threshold, "\n")


  # 16. Match edge endpoints to original node IDs

  from_index <- match(
    filtered_edges$from,
    node_points$node_id
  )

  to_index <- match(
    filtered_edges$to,
    node_points$node_id
  )

  valid_edges <- !is.na(from_index) & !is.na(to_index)

  matched_edges <- filtered_edges[valid_edges, ]

  unmatched_edges <- filtered_edges[!valid_edges, ]

  from_index <- from_index[valid_edges]
  to_index <- to_index[valid_edges]

  cat("Matched edges:", nrow(matched_edges), "\n")
  cat("Unmatched edges:", nrow(unmatched_edges), "\n")


  # 17. Create 3D edge segments

  edge_segments <- data.frame(
    x = as.vector(t(cbind(
      node_points$x[from_index],
      node_points$x[to_index],
      NA_real_
    ))),

    y = as.vector(t(cbind(
      node_points$y[from_index],
      node_points$y[to_index],
      NA_real_
    ))),

    z = as.vector(t(cbind(
      node_points$z[from_index],
      node_points$z[to_index],
      NA_real_
    )))
  )


  # 18. Select five nodes with highest Z values

  top5_nodes <- node_points[
    order(node_points$z, decreasing = TRUE)[
      seq_len(min(5, nrow(node_points)))
    ],
  ]


  # 19. Create 3D plot

  p <- plot_ly() %>%

    add_surface(
      x = x_grid,
      y = y_grid,
      z = z_matrix,
      colorscale = "Viridis",
      opacity = 0.85,
      name = "RBF surface",
      colorbar = list(title = list(text = "Mean support"))
    ) %>%

    add_trace(
      x = edge_segments$x,
      y = edge_segments$y,
      z = edge_segments$z,
      type = "scatter3d",
      mode = "lines",
      line = list(
        color = "black",
        width = 2.5
      ),
      opacity = 0.5,
      hoverinfo = "none",
      showlegend = TRUE,
      name = "Network connections",
      connectgaps = FALSE
    ) %>%

    add_markers(
      data = node_points,
      x = ~x,
      y = ~y,
      z = ~z,
      type = "scatter3d",
      mode = "markers",
      marker = list(
        color = "red",
        size = 3
      ),
      text = ~paste0(
        "Node: ", node,
        "<br>ID: ", node_id,
        "<br>Node connection: ", round(x, 2),
        "<br>Mean decision coverage: ", round(y, 3),
        "<br>Mean support: ", round(z, 2)
      ),
      hoverinfo = "text",
      name = "Original nodes"
    ) %>%

    add_trace(
      data = top5_nodes,
      x = ~x,
      y = ~y,
      z = ~z,
      type = "scatter3d",
      mode = "text",
      text = ~node,
      textposition = "top center",
      textfont = list(
        size = 7,
        color = "red"
      ),
      showlegend = FALSE
    ) %>%

    layout(
      title = paste("Topographical map:", network_name),
      scene = list(
        xaxis = list(title = "Node connection"),
        yaxis = list(title = "Mean decision coverage"),
        zaxis = list(
          title = "Mean node support",
          range = c(0, max(node_points$z) + 2)
        ),
        aspectmode = "auto"
      )
    )


  # 20. Return plot and underlying data

  return(list(
    plot = p,
    node_points = node_points,
    surface_points = surface_points,
    edges = filtered_edges,
    matched_edges = matched_edges,
    unmatched_edges = unmatched_edges,
    edge_segments = edge_segments,
    x_grid = x_grid,
    y_grid = y_grid,
    z_matrix = z_matrix,
    surface = grid,
    top5 = top5_nodes,
    condition_number = condition_number,
    sigma = sigma,
    edge_threshold = edge_threshold
  ))
}


# Generate maps

asian_map <- make_topo_plot(
  vis$Asian,
  "Asian",
  sigma = 0.05
)

cau_map <- make_topo_plot(
  vis$Cau,
  "Caucasian",
  sigma = 0.05
)

afro_map <- make_topo_plot(
  vis$Afro,
  "Afro",
  sigma = 0.05
)


# Display maps

asian_map$plot
cau_map$plot
afro_map$plot

afro_map$plot$x$data
