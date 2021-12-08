generate_overview_graph <- function() {
  sealKey()
  scCCnet <- readRDS("./data/scCCnet.rds") # S. cerevisiae Cell Cycle network
  scCCnet_ns <- scCCnet[-3] # Remove scores
  graph <- igraph::graph_from_edgelist(as.matrix(scCCnet_ns), directed = FALSE)
  degrees <- igraph::degree(graph, normalized = TRUE)
  graphxy <- igraph::layout_with_graphopt(graph, charge=0.30, mass=30)

  # Generate color palette:
  colors = heat.colors(
    max(igraph::degree(graph)+1)
  )[igraph::degree(graph)+1]

  plot(
    graph,
    layout=graphxy,
    rescale=TRUE, # Default is true, explicit as reminder.
    vertex.size = (125 * degrees),
    vertex.color = colors,
    vertex.label=NA
  )

  legend(
    "topright",
    legend=c(max(igraph::degree(graph)), 0),
    fill=c(colors[which.max(igraph::degree(graph))], colors[1]),
    title = "Degrees"
  )
  sealKey()
}

generate_log_freq_vs_log_rank <- function() {
  sealKey()
  scCCnet <- readRDS("./data/scCCnet.rds") # S. cerevisiae Cell Cycle network
  scCCnet_ns <- scCCnet[-3] # Remove scores
  graph <- igraph::graph_from_edgelist(as.matrix(scCCnet_ns), directed = FALSE)
  degrees <- igraph::degree(graph, normalized = FALSE)
  freqs <- table(degrees)
  scatter.smooth(log10(as.numeric(names(freqs)) + 1),
                 log10(as.numeric(freqs)), type = "b",
                 bg = "#A5F5CC",
                 xlab = "log(Rank)", ylab = "log(frequency)",
                 main = "Log Frequency against Log-Rank")
  sealKey()
}


get_top_degrees <- function(n = 5) {
  sealKey()
  scCCnet <- readRDS("./data/scCCnet.rds") # S. cerevisiae Cell Cycle network
  scCCnet_ns <- scCCnet[-3] # Remove scores
  graph <- igraph::graph_from_edgelist(as.matrix(scCCnet_ns), directed = FALSE)
  degrees <- igraph::degree(graph, normalized = FALSE)

  sorted_degrees <- degrees[order(degrees)]
  sealKey()
  tail(sorted_degrees, n)
}

get_bot_degrees <- function(n = 5) {
  sealKey()
  scCCnet <- readRDS("./data/scCCnet.rds") # S. cerevisiae Cell Cycle network
  scCCnet_ns <- scCCnet[-3] # Remove scores
  graph <- igraph::graph_from_edgelist(as.matrix(scCCnet_ns), directed = FALSE)
  degrees <- igraph::degree(graph, normalized = FALSE)

  sorted_degrees <- degrees[order(degrees)]
  sealKey()
  head(sorted_degrees, n)
}

