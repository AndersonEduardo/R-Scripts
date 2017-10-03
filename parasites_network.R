
---
title: "Using \\texttt{helminthR} to query and visualize host-parasite interaction networks"
author: "Tad Dallas"
includes:
  in_header:
    - \usepackage{lmodern}
output:
  pdf_document:
    fig_caption: yes
    fig_height: 6
    fig_width: 6
  html_document:
    fig_caption: yes
    fig_height: 6
    fig_width: 6
    highlight: tango
    theme: journal
---

This is a short guide to obtaining and visualizing host-parasite interaction data  from the London Natural History Museum's Host-Parasite Database using the \texttt{helminthR} R package. The resulting visualization is Figure 1 in the main text of the following article

> T. Dallas 2015. helminthR: An R interface to the London Natural History Museum's Host-Parasite Database - Ecography.


## Obtain unique interactions for a given location

```{r eval=TRUE, echo=TRUE, warning=FALSE, message=FALSE}
# load required packages
if(!require(igraph)){
    install.packages("igraph")
}

if(!require(devtools)){
    install.packages("devtools")
}

devtools::install_github('ropensci/helminthR')

library(igraph)
library(helminthR)
```

Below, I query the database for all host-parasite interactions at the host species level (no common names allowed) for Lake Erie. I then convert the output edgelist (2 column data.frame of host and parasite names) into an interaction matrix, and remove empty rows or columns. This occurs because hosts and parasites are treated as factors, and using table on a factor will display all levels of the factor.

```{r eval=TRUE, echo=TRUE}
 ErieHostPars <- unique(helminthR::findLocation("Lake Erie", speciesOnly=TRUE,
                                            hostState=1)[,1:2])

# Turning edgelist into igraph object.
 ErieAdj <- as.matrix(table(ErieHostPars[,1:2]))
 ErieAdj <- ErieAdj[rowSums(ErieAdj)>0, ]
 ErieAdj <- ErieAdj[ ,colSums(ErieAdj)>0]
```

From this point, we can turn the interaction matrix into an `igraph` object, which allows for analysis and visualization. Below, I create an `igraph` object, color the nodes based on whether they are hosts (blue) or parasites (black), suppress node labels, and alter node sizes of hosts and parasites.

```{r eval=TRUE, echo=TRUE}
ErieIg <- graph.incidence(ErieAdj)

# Some style.
V(ErieIg)$color <-  "#000000"
V(ErieIg)$color[1:nrow(ErieAdj)] <- "#1E90FF"
V(ErieIg)$label <- NA
V(ErieIg)$size <- 5
V(ErieIg)$size[1:nrow(ErieAdj)] <- 10
```

Then I plot the network.

```{r eval=TRUE, echo=TRUE, fig.cap="The host-parasite association network for Lake Erie, one of the Great Lakes located in the Northern United States. Grey lines between boxes represent interactions between hosts (larger blue dots) and helminth parasites (smaller black dots)."}  
plot(ErieIg)
```
