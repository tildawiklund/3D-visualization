library(dplyr)
library(plotly)
library(akima)
library(fields)
library(RANN)

ros <- readRDS("rosetta_results.rds")


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

