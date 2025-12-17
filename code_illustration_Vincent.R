library(tidyverse)
library(ggimage)

# -----------------------
# 1. Parameters
# -----------------------
set.seed(42)

n_x <- 28
n_y <- 28
n_var <- 5

# Spike PNGs
spike_files <- c("spike_1.png","spike_2.png","spike_3.png","spike_4.png","spike_5.png")

# -----------------------
# 2. Initial lattice
# -----------------------
grid0 <- expand.grid(x=1:n_x, y=1:n_y) %>%
  mutate(variety = sample(1:n_var, n(), replace=TRUE))

# -----------------------
# 3. Interactions
# -----------------------
interaction_matrix <- matrix(runif(n_var^2,-1,1), nrow=n_var)
diag(interaction_matrix) <- 1
interaction <- function(v1,v2) interaction_matrix[v1,v2]
neighbors <- tibble(dx=c(1,-1,0,0), dy=c(0,0,1,-1))

local_score <- function(df,i){
  cell <- df[i,]
  neighbors %>%
    mutate(x=cell$x+dx, y=cell$y+dy) %>%
    inner_join(df, by=c("x","y")) %>%
    summarise(score=sum(interaction(cell$variety,variety))) %>% pull(score)
}

total_yield <- function(df){
  yield <- 0
  for(i in 1:nrow(df)){
    cell <- df[i,]
    neigh <- neighbors %>%
      mutate(x=cell$x+dx, y=cell$y+dy) %>%
      inner_join(df, by=c("x","y"))
    yield <- yield + sum(interaction(cell$variety, neigh$variety))
  }
  ceiling(yield/2)
}

# -----------------------
# 4. Random vs Optimized
# -----------------------
grid_random <- grid0
yield_random <- total_yield(grid_random)

grid_opt <- grid0
for(step in 1:2000){
  i <- sample(nrow(grid_opt),1)
  neigh <- neighbors %>%
    mutate(x=grid_opt$x[i]+dx, y=grid_opt$y[i]+dy) %>%
    inner_join(grid_opt %>% mutate(id=row_number()), by=c("x","y"))
  if(nrow(neigh)==0) next
  j <- sample(neigh$id,1)
  before <- local_score(grid_opt,i)+local_score(grid_opt,j)
  grid_opt$variety[c(i,j)] <- grid_opt$variety[c(j,i)]
  after <- local_score(grid_opt,i)+local_score(grid_opt,j)
  if(after<before) grid_opt$variety[c(i,j)] <- grid_opt$variety[c(j,i)]
}
yield_opt <- total_yield(grid_opt)

# -----------------------
# 5. Plotting data
# -----------------------
plot_data <- bind_rows(
  grid_random %>% mutate(system=factor("Random", levels=c("Random","Optimized"))),
  grid_opt %>% mutate(system=factor("Optimized", levels=c("Random","Optimized")))
) %>%
  mutate(spike = spike_files[variety])   # assign spike PNG by variety

facet_labels <- c(
  Random = paste0("Random\nTotal yield: ", yield_random),
  Optimized = paste0("Optimized\nTotal yield: ", yield_opt)
)

# -----------------------
# 6. Plot
# -----------------------
p <- ggplot(plot_data, aes(x=x, y=y)) +

  # spike PNGs
  geom_image(aes(image=spike), size=0.04) +  # increased size
  facet_wrap(~system, nrow=1, labeller=labeller(system=facet_labels)) +
  coord_fixed(xlim=c(0.5,n_x+0.5), ylim=c(0.5,n_y+0.5), expand=TRUE) +
  labs(title="MOSAIC", subtitle="by ChatGPT prompted by Germain Montazeaud") +
  theme_void() +
  theme(
    panel.background = element_rect(fill="white", color=NA),  # white background
    legend.position="none",
    panel.spacing = unit(3, "cm"),
    strip.text = element_text(size=16, face="bold", lineheight=1.1),
    plot.title = element_text(size=40, face="bold", hjust=0.5),
    plot.subtitle = element_text(size=18, hjust=0.5, color="grey40", margin=margin(b=20)),
    plot.margin = margin(20,20,20,20),
    plot.background  = element_rect(fill="white", color=NA)
  )

p

# -----------------------
# 7. Export
# -----------------------
ggsave(
  "MOSAIC_A5_spikes_clean.png",
  p,
  width=24, height=18, units="cm", dpi=600
)
