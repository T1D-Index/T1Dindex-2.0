plot_Heatmap <- function(matrix,title,x_name="",y_name="",x_interval=0,y_interval=0)
{
  # matrix <- matrix_raw
  # matrix$y <- rownames(matrix)
  # matrix   <- gather(matrix,key=x,value="z",-y)


  matrix$x <- as.character(matrix$x)
  matrix$z <- as.numeric(matrix$z)
  heat_map <- matrix %>%
    e_charts(x) %>%
    e_heatmap(y, z) %>%
    e_visual_map(z) %>%
    e_title(title) %>%
    e_x_axis(axisLabel = list(interval = x_interval ,  rotate = 90)) %>%
    e_y_axis(axisLabel = list(interval = y_interval))  %>%
    e_axis_labels(
      x = x_name,
      y = y_name
    ) # rotate

  return(heat_map)

}


library(ggraph)
library(igraph)
Tools.tree_func <- function(final_model, tree_num)
{

  # get tree by index
  tree <- randomForest::getTree(final_model,
                                k = tree_num,
                                labelVar = TRUE) %>%
    tibble::rownames_to_column() %>%
    # make leaf split points to NA, so the 0s won't get plotted
    mutate(`split point` = ifelse(!is.na(`split var`), `split point`, NA))
  if(class(tree$prediction)=="character")
  {
  }else
  {
    tree$prediction <- round(tree$prediction,0)
  }
  # prepare data frame for graph
  graph_frame <- data.frame(from = rep(tree$rowname, 2),
                            to = c(tree$`left daughter`, tree$`right daughter`))

  # convert to graph and delete the last node that we don't want to plot
  graph <- graph_from_data_frame(graph_frame) %>%
    delete_vertices("0")

  # set node labels
  V(graph)$node_label <- gsub("_", " ", as.character(tree$`split var`))
  V(graph)$leaf_label <- (tree$prediction)
  if(class(tree$prediction)=="character" )
  {
    V(graph)$leaf_label_color <- tree$prediction
  }else
  {
    V(graph)$leaf_label_color <- as.character(  cut(as.numeric(tree$prediction ), c(quantile(tree$prediction)) ) )
  }
  V(graph)$split <- as.character(round(tree$`split point`, digits = 2))

  # plot
  plot <- ggraph(graph, 'dendrogram',circular = TRUE) +
    theme_bw() +
    geom_edge_link() +
    #geom_edge_loop0() +
    geom_node_point() +
    geom_node_text(aes(label = node_label), na.rm = TRUE, repel = TRUE) +
    geom_node_label(aes(label = split), vjust = 2.5, na.rm = TRUE,fontface = "bold", fill = "white") +
    geom_node_label(aes(label = leaf_label, fill = leaf_label_color), na.rm = TRUE,
                    repel = TRUE, colour = "black", fontface = "bold", show.legend = TRUE) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          panel.background = element_blank(),
          plot.background = element_rect(fill = "white"),
          panel.border = element_blank(),
          axis.line = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          plot.title = element_text(size = 18))
  print(plot)
}



