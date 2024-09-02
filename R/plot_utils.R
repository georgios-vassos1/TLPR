#' Create a Gradient Color Palette
#'
#' This function generates a gradient color palette of `n` colors ranging from `low` to `high`.
#' If a vector of normalized values (`norm_`) is provided, the function maps these values
#' to colors in the gradient palette.
#'
#' @param n Integer. The number of colors to generate in the gradient palette.
#' @param norm_ Numeric vector. A vector of normalized values to map to colors in the gradient palette. Default is `NULL`.
#' @param low Character. The color representing the low end of the gradient. Default is `"#80FFFF"`.
#' @param high Character. The color representing the high end of the gradient. Default is `"#FF80FF"`.
#'
#' @return A character vector of hex color codes. If `norm_` is provided, returns a vector of colors corresponding to `norm_`.
#' If `norm_` is not provided, returns the full gradient palette.
#'
#' @importFrom grDevices colorRampPalette
#'
#' @examples
#' # Generate a gradient of 5 colors
#' create_gradient_palette(5)
#'
#' # Generate a gradient and map normalized values to colors
#' create_gradient_palette(5, norm_ = c(0.1, 0.5, 0.9))
#'
#' @export
create_gradient_palette <- function(n, norm_=NULL, low = "#80FFFF", high = "#FF80FF") {
  # Create a gradient color palette
  gradient_colors <- colorRampPalette(c(low, high))(n)
  if (is.null(norm_)) return(gradient_colors)

  # Map normalized z-values to colors in the gradient palette
  gradient_colors[cut(norm_, breaks = n, include.lowest = TRUE)]
}

#' Plot State Value Surface
#'
#' This function creates a 3D plot of the state value surface at a given time `t`.
#' The z-values are normalized and mapped to colors using a gradient palette.
#'
#' @param VT Numeric matrix. The value matrix to plot.
#' @param t Integer. The time index for which the value surface is plotted.
#' @param theta Numeric. The angle defining the viewing direction. Default is `30`.
#' @param phi Numeric. The angle defining the viewing direction. Default is `30`.
#'
#' @return A 3D plot of the state value surface.
#'
#' @importFrom grDevices colorRampPalette
#'
#' @examples
#' # Assuming 'VT' is a predefined value matrix and 't' is a valid time index
#' plotVsurface(VT, t = 5)
#'
#' @export
plotVsurface <- function(env, V.t, t, theta = 30, phi = 30) {
  # Normalize the z-values to the range [0, 1] for mapping to colors
  norm_V <- (V.t - min(V.t, na.rm = T)) / (max(V.t, na.rm = T) - min(V.t, na.rm = T))

  # Map normalized z-values to colors in the gradient palette
  vertex_colors <- create_gradient_palette(1000L, sort(norm_V))

  # Create a 3D plot using plot3D
  plot3D::persp3D(x = env$stateSupport, y = env$extendedStateSupport, z = t(V.t), col = vertex_colors, theta = theta, phi = phi,
                  xlab = "Entry State", ylab = "Exit State", zlab = "V", main = paste0("State Value Surface at t = ", t))

  # Add contours to the surface
  plot3D::contour3D(x = env$stateSupport, y = env$extendedStateSupport, z = t(V.t), colvar = t(V.t), add = TRUE, col = "darkblue", colkey = FALSE, lwd = 2)
}

#' Plot Action Value Surface for Fixed Entry State
#'
#' This function creates a 3D plot of the action value surface at a given time `t` for a fixed entry state.
#' The z-values are normalized and mapped to colors using a gradient palette.
#'
#' @param Qx Numeric matrix. The action value matrix to plot.
#' @param entry_state Integer. The fixed entry state to use in the plot.
#' @param t Integer. The time index for which the action value surface is plotted.
#' @param theta Numeric. The angle defining the viewing direction. Default is `66`.
#' @param phi Numeric. The angle defining the viewing direction. Default is `30`.
#'
#' @return A 3D plot of the action value surface for the fixed entry state.
#'
#' @importFrom grDevices colorRampPalette
#'
#' @examples
#' # Assuming 'Qx' is a predefined action value matrix and 'entry_state' is defined
#' plotQsurfaceForFixedEntry(Qx, entry_state = 2, t = 5)
#'
#' @export
plotQsurfaceForFixedEntry <- function(env, Qx, entry_state, t, theta = 66, phi = 30, zlab = 'Action Value') {
  ## Plot the action value surface at time t
  norm_Q <- (Qx - min(Qx, na.rm = T)) / (max(Qx, na.rm = T) - min(Qx, na.rm = T))

  vertex_colors <- create_gradient_palette(1000L, sort(norm_Q))

  plot3D::persp3D(x = env$extendedStateSupport, y = env$actionSupport, z = t(Qx), col = vertex_colors, theta = theta, phi = phi,
          xlab = "Exit State", ylab = "Action", zlab = zlab, main = paste0("Action Value (Entry State: ", entry_state, ") at t = ", t))

  plot3D::contour3D(x = env$extendedStateSupport, y = env$actionSupport, z = t(Qx), colvar = t(Qx), add = TRUE, col = "darkblue", colkey = FALSE, lwd = 2)
}

#' Plot Action Value Surface for Fixed Exit State
#'
#' This function creates a 3D plot of the action value surface at a given time `t` for a fixed exit state.
#' The z-values are normalized and mapped to colors using a gradient palette.
#'
#' @param Qx Numeric matrix. The action value matrix to plot.
#' @param exit_state Integer. The fixed exit state to use in the plot.
#' @param t Integer. The time index for which the action value surface is plotted.
#' @param theta Numeric. The angle defining the viewing direction. Default is `66`.
#' @param phi Numeric. The angle defining the viewing direction. Default is `30`.
#'
#' @return A 3D plot of the action value surface for the fixed exit state.
#'
#' @importFrom grDevices colorRampPalette
#'
#' @examples
#' # Assuming 'Qx' is a predefined action value matrix and 'exit_state' is defined
#' plotQsurfaceForFixedExit(Qx, exit_state = 3, t = 5)
#'
#' @export
plotQsurfaceForFixedExit <- function(env, Qx, exit_state, t, theta = 66, phi = 30, zlab = 'Action Value') {
  ## Plot the action value surface at time t
  norm_Q <- (Qx - min(Qx, na.rm = T)) / (max(Qx, na.rm = T) - min(Qx, na.rm = T))

  vertex_colors <- create_gradient_palette(1000L, sort(norm_Q))

  plot3D::persp3D(x = env$stateSupport, y = env$actionSupport, z = t(Qx), col = vertex_colors, theta = theta, phi = phi,
          xlab = "Entry State", ylab = "Action", zlab = zlab, main = paste0("Action Value (Exit State: ", exit_state, ") at t = ", t))

  plot3D::contour3D(x = env$stateSupport, y = env$actionSupport, z = t(Qx), colvar = t(Qx), add = TRUE, col = "darkblue", colkey = FALSE, lwd = 2)
}

#' Plot 2x2 Instance
#'
#' This function plots a series of graphs representing the state of the environment at different times.
#'
#' @param env A list containing the environment variables including `tau`, `nI`, and `nJ`.
#' @param graph.dt A data.table containing the graph data with columns `t`, `origin`, `destination`, and `assignment`.
#' @param S.I A vector representing the source nodes.
#' @param S.J A vector representing the destination nodes.
#' @param Q A vector representing some quantities related to the nodes.
#' @param D A vector representing some other quantities related to the nodes.
#' 
#' @import data.table igraph
#' @examples
#' # Assuming `env`, `graph.dt`, `S.I`, `S.J`, `Q`, and `D` are defined
#' plot2x2instance(env, graph.dt, S.I, S.J, Q, D)
#'
#' @export
plot2x2instance <- function(env, graph.dt, S.I, S.J, Q, D) {
  n <- 2L
  reflex <- diag(n)[, seq(n, 1L)]
  
  # Set up the plotting area to display multiple plots
  par(mfrow = c(3L, 4L))
  for (t_ in seq(env$tau)) {
    # Create an igraph object for the current time step
    g <- graph_from_data_frame(graph.dt[t == t_, -1L], directed = TRUE, vertices = seq(env$nI + env$nJ))
    
    # Set vertex attributes
    V(g)$color <- "lightblue"
    V(g)$size <- 50L  # Size of the vertices
    
    # Set edge attributes
    E(g)$arrow.size <- 0.25  # Size of the arrows
    E(g)$weight <- graph.dt$assignment[(t_ - 1L) * env$nL + seq(env$nL)]
    
    # Plot the graph with edge labels
    plot(g, edge.label = E(g)$weight, layout = igraph::layout_on_grid(g) %*% reflex, main = paste0("t = ", t_))
    
    # Add text labels to the plot
    text(x = c(-1.34, -1.34, 1.34, 1.34, -1.00, -1.00, 1.00, 1.00), 
         y = c(-1.00, 1.00, -1.00, 1.00, -0.60, 0.60, -0.60, 0.60), 
         labels = c(
           Q[(t_ - 1L) * env$nI + seq(env$nI)],
           D[(t_ - 1L) * env$nI + seq(env$nI)],
           S.I[(t_ - 1L) * env$nI + env$I_, 1L],
           S.J[(t_ - 1L) * env$nJ + env$J_, 1L]
         ),
         adj = c(.5, .5))
  }
  # Reset plotting area to single plot
  par(mfrow = c(1L, 1L))
}
