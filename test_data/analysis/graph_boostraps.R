library(dplyr)
library(ggplot2)
library(cowplot)

"%||%" <- function(a, b) {
  if (!is.null(a)) a else b
}

geom_flat_violin <- function(mapping = NULL, data = NULL, stat = "ydensity",
                             position = "dodge", trim = TRUE, scale = "area",
                             show.legend = NA, inherit.aes = TRUE, ...) {
  layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = GeomFlatViolin,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      trim = trim,
      scale = scale,
      ...
    )
  )
}

GeomFlatViolin <- ggproto("GeomFlatViolin", Geom,
                          setup_data = function(data, params) {
                            data$width <- data$width %||%
                              params$width %||% (resolution(data$x, FALSE) * 0.9)
                            
                            # ymin, ymax, xmin, and xmax define the bounding rectangle for each group
                            data %>%
                              group_by(group) %>%
                              mutate(ymin = min(y),
                                     ymax = max(y),
                                     xmin = x,
                                     xmax = x + width / 2)
                            
                          },
                          draw_group = function(data, panel_scales, coord) {
                            # Find the points for the line to go all the way around
                            data <- transform(data, xminv = x,
                                              xmaxv = x + violinwidth * (xmax - x))
                            
                            # Make sure it's sorted properly to draw the outline
                            newdata <- rbind(plyr::arrange(transform(data, x = xminv), y),
                                             plyr::arrange(transform(data, x = xmaxv), -y))
                            
                            # Close the polygon: set first and last point the same
                            # Needed for coord_polar and such
                            newdata <- rbind(newdata, newdata[1,])
                            
                            ggplot2:::ggname("geom_flat_violin", GeomPolygon$draw_panel(newdata, panel_scales, coord))
                          },
                          draw_key = draw_key_polygon,
                          default_aes = aes(weight = 1, colour = "grey20", fill = "white", size = 0.5,
                                            alpha = NA, linetype = "solid"),
                          required_aes = c("x", "y")
)



file = '/Users/ulisesrosas/Desktop/qsin/test_data/n15_depth3_v0.9_n0.9_all/compared_nets_n15_depth3_v0.9_n0.9_all.csv'
combs_file = '/Users/ulisesrosas/Desktop/qsin/test_data/n15_depth3_v0.9_n0.9_all/n15_depth3_v0.9_n0.9_overlappedBatches_uniq.txt'
time_file ='/Users/ulisesrosas/Desktop/qsin/test_data/n15_depth3_v0.9_n0.9_all/processed_time_n15_depth3_v0.9_n0.9_all.csv'
pseudo_file ='/Users/ulisesrosas/Desktop/qsin/test_data/n15_depth3_v0.9_n0.9_all/processed_pseudo_n15_depth3_v0.9_n0.9_all.csv'

is_linear = F



pseudo_liks <- read.csv(pseudo_file)
starts <- min(table(pseudo_liks$row))


df_lik = pseudo_liks %>%
  mutate(id = paste(row, boot, sep = ",")) %>%
  dplyr::select(id, likdev)


df_time = read.csv(time_file) %>%
  dplyr::rename( boot = bootstrap)%>%
  mutate(id = paste(row, boot, sep = ",")) %>%
  dplyr::select(id, time)


if (is_linear){
  msubtitle_d = paste('Elastic Net,', starts, 'random starts')
}else{
  msubtitle_d = paste('ISLE,', starts, 'random starts')
}


max_time = 14827.94
CT_combs <- read.csv(combs_file, sep = "\t", header = F) 

lapply(CT_combs$V1, function (x){ length(strsplit(x =  x, split = ",")[[1]]) }) %>%
  unlist() ->  CT_row_count


df <- read.csv(file, header = F) %>%
  dplyr::group_by(V1, V2) %>%
  dplyr::arrange(V2,V1) %>%
  as.data.frame(.)

colnames(df) <- c("row", 'boot', 'dist', 'root')

df <- df %>%
  mutate(id = paste(row, boot, sep = ",")) %>%
  dplyr::left_join(., df_lik, by = 'id') %>%
  dplyr::left_join(., df_time, by = 'id') %>%
  dplyr::group_by(row) %>%
  # dplyr::mutate(rel_likdev = (likdev - mean(likdev))/sd(likdev))
  dplyr::mutate(rel_likdev = (likdev - min(likdev))/(max(likdev) - min(likdev))) %>%
  dplyr::ungroup()


# 
# 
test_rows <- unique(df$row)
# 
lapply(test_rows, function(r){
  tmp_df = df[df$row == r,]
  # r
  tmp_df$CT_row_count <- CT_row_count[r]
  tmp_df
}) %>%
  do.call('rbind', .) -> df


df$CT_row_count <- factor(df$CT_row_count, levels = sort(unique(df$CT_row_count)))
df <- df[!is.na(df$dist),]

df[df$row == 1,]


df %>%
  ggplot(data = ., aes(x = CT_row_count
                       , y = dist)) +
  geom_point(
    aes(color = rel_likdev),
    position = position_jitter(width = 0.08),
    size = 1.5, 
    alpha = 0.6
  ) +
  scale_colour_gradient(low = "blue", high = "red")+
  geom_boxplot(width = 0.2, 
               outlier.shape = NA, 
               alpha = 0,
               position = position_nudge(x = 0, y = 0)) +
  geom_flat_violin(width=0.5, 
                   position = position_nudge(x = 0.14, y = 0),
                   alpha = 2) +
  theme_bw(base_size = 12) +
  labs(title="Distance to the true network",
       subtitle = msubtitle_d,
       color = 'Relative\ndeviance',
       y ="Distance", x = "Rows in CT")+
  geom_hline(yintercept=0, linetype="dashed", color = "red") -> dist_plot
# nboots <- length(unique(df$boot))



df %>%
  ggplot(data = ., 
         aes(x = CT_row_count,
             y = time/60/60)) +
  geom_point(aes(colour = rel_likdev), position = position_jitter(width = 0.10), size = 1, alpha = 0.5, )+
  scale_colour_gradient(low = "blue", high = "red")+
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0, position = position_nudge(x = 0, y = 0)) +
  geom_flat_violin(width=0.5,  position = position_nudge(x = 0.14, y = 0), alpha = 2) +
  theme_bw(base_size = 12) +
  labs(title="Running time",
       subtitle = msubtitle_d,
       y ="Hours", x = "Rows in CT")+
  theme(legend.position="none") +
  geom_hline(yintercept=max_time/60/60, linetype="dashed", color = "red")  -> time_plot


# pdf(
#   file = '../../ABL/software/qsin/test_data/n15_depth3_v0.9_n0.5/dist_time_nonlinear.pdf',
#   width = 6,height = 5
#     )
plot_grid(
  dist_plot,
  time_plot,
  ncol = 1,
  align = "v"
)

lapply(unique(df$row), function(x) wilcox.test(jitter(df[df$row ==x,]$dist)  )$p.value)



