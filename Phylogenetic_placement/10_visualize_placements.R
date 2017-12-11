#!/usr/bin/Rscript
library(tidyverse)
library(stringr)
library(magrittr)
library(RColorBrewer)
library(ggtree)
library(ggpubr)


# functions

jplace2tidytrees <- function(jplace, outgroup_tips) {
  
  jplace@phylo$edge.length <- tibble(node = jplace@phylo$edge[, 2]) %>%
    left_join(attr(jplace@phylo, "edgeNum")) %>%
    pull(edgeNum)
  
  outgroup_node <- MRCA(jplace@phylo, tip = outgroup_tips)
  
  jplace@phylo <- reroot(jplace@phylo, outgroup_node)
  
  nodes <- tibble()
  
  for (taxon_index in 1:length(jplace@placements$n)) {
    
    taxon_name <- jplace@placements$n[[taxon_index]]
    
    nodes_new <- as.data.frame(jplace@phylo, branch.length = "none") %>%
      rename(edge_num = length)
    
    placements <- jplace@placements$p[[taxon_index]] %>%
      as_tibble() %>%
      set_names(jplace@fields)
    
    nodes_new <- left_join(nodes_new, placements) %>%
      mutate(placement = ifelse(is.na(like_weight_ratio), "no", "yes")) %>%
      mutate(asv = !! taxon_name)
    
    nodes <- bind_rows(nodes, nodes_new)
    
  }
  
  list(tidytrees = nodes, tree = jplace@phylo)
  
}

# preparation

if (! dir.exists("results")) dir.create("results")

jplace <- read.jplace(file = "data/tree_placement/RAxML_portableTree.placement.jplace")

outgroup_tips <- c(
  "GCA_001437585.1_1",
  "GCA_001438885.1_1",
  "GCA_001437015.1_1",
  "GCA_001437015.1_2"
)

out <- jplace2tidytrees(jplace, outgroup_tips)
tidytrees <- out$tidytrees
tree <- out$tree

genomes <- read_tsv(file = "data/genomes.tsv") 

# individual ASVs

tidytrees %>%
  separate(asv, into = c("genus", "number"), sep = "_", remove = F) %>%
  mutate(number = as.integer(number)) %>%
  arrange(genus, number) %>%
  mutate(asv = str_replace_all(asv, "_", " ")) %>%
  mutate(asv = factor(asv, levels = unique(asv))) %>%
  ggtree(layout = "circular", aes(col = placement), size = 0.3, branch.length = "none") +
  scale_color_manual(values = c(no = "#D0D0D0", yes = "#d95f02")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  facet_wrap(~ asv)

ggsave(filename = "/mnt/H/Drafts/Sander/Fermentatie/Supplementary/SupFigure4.tiff",
       width = 178,
       height = 210,
       units = "mm",
       dpi = 300)


# large tree without clades

asvs_best_node <- tidytrees %>%
  group_by(asv) %>%
  filter(like_weight_ratio == max(like_weight_ratio, na.rm = T)) %>%
  slice(1) %>%
  ungroup() %>%
  select(asv, best_node = node)

tidytree <- tidytrees %>%
  group_by(node, parent, edge_num, x, y, label, isTip, branch, angle) %>%
  summarize(placement = sum(! is.na(like_weight_ratio)) > 0) %>%
  mutate(placement = ifelse(placement, "yes", "no")) %>%
  mutate(best_placement = ifelse(node %in% asvs_best_node$best_node, "yes", "no")) %>%
  mutate(assembly_accession = str_extract(label, "[A-Z]+_[0-9]+\\.[0-9]")) %>%
  left_join(genomes) %>%
  mutate(species = str_replace_all(species, "_", " ")) %>%
  mutate(species = str_extract(species, "[^ ]+ [^ ]+"))

tidytree %>%
  ggtree(layout = "circular", aes(col = placement), size = 0.5, branch.length = "none") +
  geom_tiplab2(aes(label = species), col = "#000000", size = 1.5) +
  scale_color_manual(values = c("no" = "#D0D0D0", "yes" = "#000000")) +
  theme(plot.title = element_text(hjust = 0.5)) 

ggsave("results/placement_no_clades.png", units = "cm", width = 20, height = 20)

# find the node numbers of clade MRCAs

tidytree  %>%
  ggtree(layout = "circular", aes(col = placement), size = 0.5, branch.length = "none") +
  geom_tiplab2(aes(label = node), col = "#000000", size = 3) +
  scale_color_manual(values = c(no = "#D0D0D0", yes = "#d95f02")) +
  theme(plot.title = element_text(hjust = 0.5)) 

clades <- tibble(
    tip1 = c(85, 61, 39, 31, 185, 170, 165, 120, 114, 90),
    tip2 = c(62, 52, 32, 28, 178, 166, 161, 115, 111, 94),
    clade_name = c(
      "Leuconostoc",
      "L. brevis group", 
      "L. kunkei group",
      "L. coryniformis group",
      "L. salivarius group",
      "L. sakei group",
      "L. casei group",
      "L. plantarum group",
      "L. vaccinostercus group",
      "Weissella"
    )
  ) %>%
  mutate(mrca = map2_int(tip1, tip2, function(x, y) MRCA(tree, c(x, y)))) %>%
  mutate(col = brewer.pal(nrow(.), name = "Paired") %>% as.character())

nodes_clades <- clades %>%
  group_by(mrca, clade_name) %>%
  do({
    tibble(node = c(getDescendants(tree = tree, node = .$mrca), .$mrca))
  }) %>%
  ungroup() %>%
  rename(clade = mrca)

colors <- clades$col %>%
  set_names(clades$clade_name %>% as.character()) %>%
  c(., "other" = "#D0D0D0") %>%
  c(., "unclassified" = "#000000")

# large tree with clades on branches

tree_plot <- tidytree  %>%
  left_join(nodes_clades) %>%
  mutate(clade_name = ifelse(is.na(clade_name), "unclassified", clade_name)) %>%
  mutate(color = ifelse(placement == "yes", clade_name, "other")) %>%
  ggtree(layout = "circular", aes(col = color), size = 1, branch.length = "none") +
  geom_tiplab2(aes(label = species), col = "#303030", size = 1.2, offset = 0.8) +
  scale_color_manual(values = colors,
                     breaks = c("L. brevis group",
                                "L. casei group",
                                "L. coryniformis group",
                                "L. plantarum group",
                                "L. sakei group",
                                "L. salivarius group",
                                "L. vaccinostercus group",
                                "Leuconostoc",
                                "Weissella",
                                "unclassified"),
                     name = "",
                     labels = c(expression(paste(italic("Lactobacillus brevis")," group")),
                                expression(paste(italic("Lactobacillus casei")," group")),
                                expression(paste(italic("Lactobacillus coryniformis")," group")),
                                expression(paste(italic("Lactobacillus plantarum")," group")),
                                expression(paste(italic("Lactobacillus sakei")," group")),
                                expression(paste(italic("Lactobacillus salivarius")," group")),
                                expression(paste(italic("Lactobacillus vaccinostercus")," group")),
                                expression(italic("Leuconostoc")),
                                expression(italic("Weissella")),
                                "Unclassified")) +
  theme(
    legend.position = "none",
    legend.text.align = 0
  )

tree_plot

ggsave("results/placement_clades.png", units = "cm", width = 30, height = 20)

# Figure out to which group the best node per ASV belongs

count <- asvs_best_node %>%
  left_join(nodes_clades, by = c("best_node" =  "node")) %>%
  group_by(clade_name) %>%
  add_count() %>%
  select(clade_name, n) %>%
  distinct() %>%
  ggplot(aes(x = reorder(clade_name, n), y = n, fill = clade_name)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_manual(values = colors,
                     breaks = c("L. brevis group",
                                "L. casei group",
                                "L. coryniformis group",
                                "L. plantarum group",
                                "L. sakei group",
                                "L. salivarius group",
                                "L. vaccinostercus group",
                                "Leuconostoc",
                                "Weissella"),
                     name = "") +
  scale_y_continuous(breaks = c(0,2,4,6,8)) +
  scale_x_discrete(breaks = c("L. brevis group",
                                "L. casei group",
                                "L. coryniformis group",
                                "L. plantarum group",
                                "L. sakei group",
                                "L. salivarius group",
                                "L. vaccinostercus group",
                                "Leuconostoc",
                                "Weissella"),
                     labels = c(expression(paste(italic("Lactobacillus brevis")," group")),
                                expression(paste(italic("Lactobacillus casei")," group")),
                                expression(paste(italic("Lactobacillus coryniformis")," group")),
                                expression(paste(italic("Lactobacillus plantarum")," group")),
                                expression(paste(italic("Lactobacillus sakei")," group")),
                                expression(paste(italic("Lactobacillus salivarius")," group")),
                                expression(paste(italic("Lactobacillus vaccinostercus")," group")),
                                expression(italic("Leuconostoc")),
                                expression(italic("Weissella"))))+
  ylab("Count") +
  xlab("") + 
  theme(
    # Set text size
    axis.title.x = element_text(size = 12, 
                                margin = margin(5,0,0,0)),
    
    axis.text.x = element_text(size = 12, 
                               colour="black",
                               margin = margin(2,0,0,0)),
    axis.text.y = element_text(size = 10,
                               colour="black",
                               margin = margin(0,0,0,0),
                               face = "italic"),

    # Configure lines and axes
    axis.ticks.x = element_line(colour = "black"), 
    axis.ticks.y = element_blank(), 
    
    # Plot background
    panel.background = element_rect(fill = "white"),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    
    # Remove border
    panel.border = element_rect(fill="transparent",color=NA),
    
    # Remove legend
    legend.position =  "none"
    )

ggarrange(tree_plot,
          count,
          nrow = 2,
          heights = c(4, 1),
          labels = c("A", "B"))

ggsave("results/placement_clades_count.png", units = "cm", width = 30, height = 30)

ggsave(filename = "/mnt/H/Drafts/Sander/Fermentatie/Paperpictures/Figure3.tiff",
       width = 178,
       height = 210,
       units = "mm",
       dpi = 300)

# Add livestyles based on Walter

livestyles <- tibble(
  tip1 = c(1, 6, 9, 207, 205, 32, 38, 40, 41, 52, 97, 99, 101, 111, 115, 136, 161, 166, 179, 183, 184, 188, 192, 201, 200, 197),
  tip2 = c(5, 6, 25, 208, 206, 37, 39, 40, 51, 61, 97, 100, 110, 114, 120, 141, 165, 170, 180, 183, 184, 190, 194, 201, 198, 197),
  clade_name = c(
    "Vertebrate-adapted",
    "Vertebrate-adapted",
    "Vertebrate-adapted",
    "Free-living",
    "Free-living",
    "Insect-adapted",
    "Insect-adapted",
    "Vertebrate-adapted",
    "Free-living",
    "Free-living",
    "Vertebrate-adapted",
    "Vertebrate-adapted",
    "Free-living",
    "Free-living",
    "Nomadic",
    "Free-living",
    "Nomadic",
    "Free-living",
    "Free-living",
    "Free-living",
    "Free-living",
    "Vertebrate-adapted",
    "Vertebrate-adapted",
    "Vertebrate-adapted",
    "Vertebrate-adapted",
    "Vertebrate-adapted"
  )
) %>%
  mutate(colour = if_else(clade_name == "Free-living", "#b3e2cd", if_else(
    clade_name == "Nomadic", "#fdcdac", if_else(
      clade_name == "Vertebrate-adapted", "#cbd5e8", "#f4cae4")))) %>%
  mutate(mrca = map2_int(tip1, tip2, function(x, y) {
    if (x == y) return(as.integer(x))
    MRCA(tree, c(x, y))
  }))


tree_plot_clades <- tidytree  %>%
  left_join(nodes_clades) %>%
  mutate(clade_name = ifelse(is.na(clade_name), "unclassified", clade_name)) %>%
  mutate(color = ifelse(placement == "yes", clade_name, "other")) %>%
  ggtree(layout = "circular", aes(col = color), size = 1, branch.length = "none") +
  scale_color_manual(values = colors,
                     breaks = c("L. brevis group",
                                "L. casei group",
                                "L. coryniformis group",
                                "L. plantarum group",
                                "L. sakei group",
                                "L. salivarius group",
                                "L. vaccinostercus group",
                                "Leuconostoc",
                                "Weissella",
                                "unclassified"),
                     name = "") +
  theme(
    legend.position = "none",
    legend.text = element_text(face = "italic")
  ) 

for (row in 1:nrow(livestyles)) {
  tree_plot_clades <- tree_plot_clades +
    geom_cladelabel(node = livestyles$mrca[row], color = livestyles$colour[row], offset = 1, align = T, label = "", barsize = 4)
}

tree_plot_clades 

# WARNING: NOT EVERYTHING IS PLOTTED!!! KAN NOG BETER!!!

ggsave("results/placement_livestyles.png", units = "cm", width = 30, height = 20)

# make legend
legend <- ggplot(livestyles, aes(x = clade_name)) +
  geom_bar(aes(fill = clade_name)) +
  scale_fill_manual(labels = unique(livestyles$clade_name),
                      values = unique(livestyles$colour)) +
  theme(legend.title = element_blank())
legend <- get_legend(legend)

ggarrange(legend)

ggsave("results/legend_livestyles.png", units = "cm", width = 10, height = 10)


