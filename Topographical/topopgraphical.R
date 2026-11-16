library(dplyr)
library(plotly)
library(akima)
library(fields)
library(RANN)
library(tidyr)

ros <- readRDS("C:/Users/sigri/Documents/GitHub/3D-visualization/Topographical/rosetta_results.rds")

# Check the number of unique values and coordinate pairs
length(unique(ros$main$accuracyRHS))
length(unique(ros$main$supportLHS))
length(unique(ros$main$supportRHS))

nrow(unique(ros$main[, c("accuracyRHS", "supportLHS")]))
sum(duplicated(ros$main[, c("accuracyRHS", "supportLHS")]))

# 1. Prep the data

# Original rules from the data with their feature and decision info
rule_points <- ros$main %>%
  transmute(x = accuracyRHS, y = supportLHS, z = coverageRHS, features, decision)

# Average coverage for rules with identical X/Y coordinates
points <- aggregate(z ~ x + y, data = rule_points, FUN = mean)

head(points)
summary(points)

# 2. Explore Akima interpolation

# Default Akima interpolation
surface <- interp(
  x = points$x,
  y = points$y,
  z = points$z
)


# Plot the interpolated surface
plot_ly() %>%
  add_surface(
    x = surface$x,
    y = surface$y,
    z = surface$z
  )

# Compare with linear interpolation
surface_linear <- interp(
  x = points$x,
  y = points$y,
  z = points$z,
  linear = TRUE
)

# Check how different the surfaces are
max(surface$z, na.rm = TRUE)
max(surface_linear$z, na.rm = TRUE)
max(abs(surface$z - surface_linear$z), na.rm = TRUE)


# 3. Explore thin-plate spline interpolation (Tps)

# Create a regular prediction grid
grid_size <- 50

prediction_grid <- expand.grid(
  x = seq(min(points$x), max(points$x), length.out = grid_size),
  y = seq(min(points$y), max(points$y), length.out = grid_size)
)

# Fit the default thin-plate spline model
rbf_model <- Tps(
  x = as.matrix(points[, c("x", "y")]),
  Y = points$z
)

# Predict coverage across the grid
prediction_grid$z <- predict(
  rbf_model,
  as.matrix(prediction_grid[, c("x", "y")])
)


# Convert predictions into a matrix for plotly
z_matrix <- matrix(
  prediction_grid$z,
  nrow = grid_size,
  ncol = grid_size
)

# Plot the default Tps surface
plot_ly(
  x = sort(unique(prediction_grid$x)),
  y = sort(unique(prediction_grid$y)),
  z = z_matrix
) %>%
  add_surface()


# Less smoothing: lambda = 0.001
rbf_model_less_smooth <- Tps(
  x = as.matrix(points[, c("x", "y")]),
  Y = points$z,
  lambda = 0.001
)

prediction_grid$z_001 <- predict(
  rbf_model_less_smooth,
  as.matrix(prediction_grid[, c("x", "y")])
)

z_matrix_001 <- matrix(
  prediction_grid$z_001,
  nrow = grid_size,
  ncol = grid_size
)

plot_ly(
  x = sort(unique(prediction_grid$x)),
  y = sort(unique(prediction_grid$y)),
  z = z_matrix_001
) %>%
  add_surface()


# Even less smoothing: lambda = 0.00001
rbf_model_less_smooth2 <- Tps(
  x = as.matrix(points[, c("x", "y")]),
  Y = points$z,
  lambda = 0.00001
)

prediction_grid$z <- predict(
  rbf_model_less_smooth2,
  as.matrix(prediction_grid[, c("x", "y")])
)

z_matrix_less_smooth2 <- matrix(
  prediction_grid$z,
  nrow = grid_size,
  ncol = grid_size
)

plot_ly(
  x = sort(unique(prediction_grid$x)),
  y = sort(unique(prediction_grid$y)),
  z = z_matrix_less_smooth2
) %>%
  add_surface()


# 4. Show some info upon hovering about the closest rule

# Find the nearest original rule for each grid point
# Extract grid and original rule coordinates
grid_coords <- as.matrix(prediction_grid[, c("x", "y")])
rule_coords <- as.matrix(rule_points[, c("x", "y")])

# Normalise both coords using the same original coordinate ranges
x_range <- range(points$x)
y_range <- range(points$y)

grid_scaled <- scale(
  grid_coords,
  center = c(x_range[1], y_range[1]),
  scale = c(diff(x_range), diff(y_range))
)

rules_scaled <- scale(
  rule_coords,
  center = c(x_range[1], y_range[1]),
  scale = c(diff(x_range), diff(y_range))
)

# Find the nearest original rule with euclidean distance
nearest_rules <- nn2(
  data = rules_scaled,
  query = grid_scaled,
  k = 1
)$nn.idx[, 1]

head(nearest_rules)
length(nearest_rules)

# Attach information from the nearest original rule
prediction_grid$features <- rule_points$features[nearest_rules]
prediction_grid$decision <- rule_points$decision[nearest_rules]

prediction_grid$rule_accuracy <- rule_points$x[nearest_rules]
prediction_grid$rule_support <- rule_points$y[nearest_rules]
prediction_grid$rule_coverage <- rule_points$z[nearest_rules]

# 5. Add rule information to the Tps surface hover


hover_matrix <- matrix(
  paste0(
    "Accuracy: ", round(prediction_grid$x, 3),
    "<br>Support: ", prediction_grid$y,
    "<br>Estimated coverage: ", round(prediction_grid$z, 3),
    "<br><br>Nearest rule:",
    "<br>Features: ", prediction_grid$features,
    "<br>Decision: ", prediction_grid$decision,
    "<br>Original accuracy: ", round(prediction_grid$rule_accuracy, 3),
    "<br>Original support: ", prediction_grid$rule_support,
    "<br>Original coverage: ", round(prediction_grid$rule_coverage, 3)
  ),
  nrow = grid_size,
  ncol = grid_size
)

# Plot the Tps surface with hover info
plot_ly(
  x = sort(unique(prediction_grid$x)),
  y = sort(unique(prediction_grid$y)),
  z = z_matrix_less_smooth2
) %>%
  add_surface(
    text = hover_matrix,
    hoverinfo = "text",
    name = "Tps surface"
  )


# Plot the Akima interpolation with hover info

# Create a grid matching the Akima surface
akima_grid <- expand.grid(
  x = surface$x,
  y = surface$y
)

akima_coords <- as.matrix(akima_grid[, c("x", "y")])

# Normalise the Akima grid
akima_grid_scaled <- scale(
  akima_coords,
  center = c(x_range[1], y_range[1]),
  scale = c(diff(x_range), diff(y_range))
)

# Find the nearest original rule for each Akima grid point
akima_nearest_rules <- nn2(
  data = rules_scaled,
  query = akima_grid_scaled,
  k = 1
)$nn.idx[, 1]

# Attach nearest rule info
akima_grid$features <- rule_points$features[akima_nearest_rules]
akima_grid$decision <- rule_points$decision[akima_nearest_rules]

akima_grid$rule_accuracy <- rule_points$x[akima_nearest_rules]
akima_grid$rule_support <- rule_points$y[akima_nearest_rules]
akima_grid$rule_coverage <- rule_points$z[akima_nearest_rules]

# Create hover text for the Akima surface
akima_hover_matrix <- matrix(
  paste0(
    "Accuracy: ", round(akima_grid$x, 3),
    "<br>Support: ", akima_grid$y,
    "<br>Estimated coverage: ", round(as.vector(surface$z), 3),
    "<br><br>Nearest rule:",
    "<br>Features: ", akima_grid$features,
    "<br>Decision: ", akima_grid$decision,
    "<br>Original accuracy: ", round(akima_grid$rule_accuracy, 3),
    "<br>Original support: ", akima_grid$rule_support,
    "<br>Original coverage: ", round(akima_grid$rule_coverage, 3)
  ),
  nrow = length(surface$x),
  ncol = length(surface$y)
)

# Plot the Akima surface with hover info
plot_ly(
  x = surface$x,
  y = surface$y,
  z = surface$z
) %>%
  add_surface(
    text = akima_hover_matrix,
    hoverinfo = "text",
    name = "Akima surface"
  )

###
### ADD NODES AS FEATURES FOR THE PLOT, CODE BY cLAUDE
###

# Settings
rds_path    <- "C:/Users/sigri/Documents/GitHub/3D-visualization/Topographical/rosetta_results.rds"
grid_size   <- 50      # nr grid points per axis
lambda      <- NA      # Tps: NA = chosen automatically (GCV). 
# Smaller value (eg. 0.001)
# produces a surface that follows the points more closely, meaning a less even surface
show_points <- TRUE    # TRUE = shows the 20 features as black dots on the surface

ros <- readRDS(rds_path)


# Rules --> features
# features/levels/cuts are comma separated strings that are stored as pairs,
# eg. "UTS2_Th17_48,IFITM3_Th17_48". One row per feature for each rule.
feature_rows <- ros$main %>%
  mutate(rule_id = row_number()) %>%
  separate_rows(features, levels, cuts, sep = ",")

cat("Rows after separation:", nrow(feature_rows), "\n")
cat("Nr unique features: ", n_distinct(feature_rows$features), "\n")

# Hur många regler per beslutsklass som innehåller varje feature (för hovertexten)
# for hover text, how many rules per decision class that includes all features
decision_info <- feature_rows %>%
  distinct(rule_id, features, decision) %>%
  count(features, decision) %>%
  group_by(features) %>%
  summarise(decisions = paste0(decision, ": ", n, collapse = ", "),
            .groups = "drop")

# One point per feature: mean value of all rules that include the feature
feature_points <- feature_rows %>%
  group_by(features) %>%
  summarise(
    x       = mean(accuracyRHS),
    y       = mean(supportLHS),
    z       = mean(coverageRHS),
    n_rules = n_distinct(rule_id),
    .groups = "drop"
  ) %>%
  left_join(decision_info, by = "features")

# Control: no features can have the exact same (x, y)
cat("Dubles in (x, y):", sum(duplicated(feature_points[, c("x", "y")])), "\n")

feature_points$hover <- paste0(
  "<b>", feature_points$features, "</b>",
  "<br>Antal regler: ", feature_points$n_rules,
  "<br>Regler per beslut: ", feature_points$decisions,
  "<br>Medel-accuracy: ", round(feature_points$x, 3),
  "<br>Medel-support: ", round(feature_points$y, 1),
  "<br>Medel-coverage: ", round(feature_points$z, 3)
)


# Help Functions
# Find the closest feature for each grid point (x and y scaled to 0-1 first, 
# otherwise support would dominate the distance)
nearest_feature <- function(grid) {
  x_range <- range(feature_points$x)
  y_range <- range(feature_points$y)
  scale01 <- function(m) {
    cbind((m[, 1] - x_range[1]) / diff(x_range),
          (m[, 2] - y_range[1]) / diff(y_range))
  }
  nn2(
    data  = scale01(as.matrix(feature_points[, c("x", "y")])),
    query = scale01(as.matrix(grid[, c("x", "y")])),
    k = 1
  )$nn.idx[, 1]
}

# Hover text fro each point on the surdface: position + closest feature
make_hover <- function(grid, z_est) {
  idx <- nearest_feature(grid)
  paste0(
    "Position på ytan:",
    "<br>Accuracy: ", round(grid$x, 3),
    "<br>Support: ", round(grid$y, 1),
    "<br>Skattad coverage: ", round(z_est, 3),
    "<br><br>Närmaste feature:<br>", feature_points$hover[idx]
  )
}

# Plotly assumes matrices with y as rows and x as columns.
# for vectors from expand.grid(x,y), x has the "fastest variation (maybe steepest curve?)", 
# because of this, byrow = True
to_matrix <- function(v, x_vals, y_vals) {
  matrix(v, nrow = length(y_vals), ncol = length(x_vals), byrow = TRUE)
}

plot_feature_surface <- function(x_vals, y_vals, z_mat, hover_mat, title) {
  p <- plot_ly() %>%
    add_surface(
      x = x_vals, y = y_vals, z = z_mat,
      text = hover_mat, hoverinfo = "text",
      colorbar = list(title = "Coverage")
    ) %>%
    layout(
      title = title,
      scene = list(
        xaxis = list(title = "Mean-accuracy"),
        yaxis = list(title = "Mean-support"),
        zaxis = list(title = "Mean-coverage")
      )
    )
  
  if (show_points) {
    p <- p %>%
      add_trace(
        data = feature_points,
        x = ~x, y = ~y, z = ~z,
        type = "scatter3d", mode = "markers",
        text = ~hover, hoverinfo = "text",
        marker = list(size = 4, color = "black"),
        name = "Features"
      )
  }
  p
}


# Surface med thin-plate spline (Tps) 
x_vals <- seq(min(feature_points$x), max(feature_points$x), length.out = grid_size)
y_vals <- seq(min(feature_points$y), max(feature_points$y), length.out = grid_size)
grid_tps <- expand.grid(x = x_vals, y = y_vals)

tps_model <- Tps(
  x = as.matrix(feature_points[, c("x", "y")]),
  Y = feature_points$z,
  lambda = lambda,
  scale.type = "range"   # scales x and y to the same interval before the spline is adjusted
)

grid_tps$z <- as.vector(predict(tps_model, as.matrix(grid_tps[, c("x", "y")])))

plot_feature_surface(
  x_vals, y_vals,
  z_mat     = to_matrix(grid_tps$z, x_vals, y_vals),
  hover_mat = to_matrix(make_hover(grid_tps, grid_tps$z), x_vals, y_vals),
  title     = "Features: Tps-surface"
)


# Surface with Akima-interpolation
# Akima ger NA utanför de 20 punkternas konvexa hölje, så ytan kan få "hål" i kanterna
# Akima sets NA outside the convex surface of the 20 points, meaning the edge might have visual holes
ip <- interp(
  x = feature_points$x,
  y = feature_points$y,
  z = feature_points$z,
  nx = grid_size, ny = grid_size
)

grid_akima   <- expand.grid(x = ip$x, y = ip$y)
grid_akima$z <- as.vector(ip$z)   # ip$z has x on the rows; as.vector gives x fastest

plot_feature_surface(
  ip$x, ip$y,
  z_mat     = to_matrix(grid_akima$z, ip$x, ip$y),
  hover_mat = to_matrix(make_hover(grid_akima, grid_akima$z), ip$x, ip$y),
  title     = "Features: Akima-surface"
)
