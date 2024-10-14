library(dplyr)
library(ggplot2)
library(cowplot)
source("https://gist.githubusercontent.com/benmarwick/2a1bb0133ff568cbe28d/raw/fb53bd97121f7f9ce947837ef1a4c65a73bffb3f/geom_flat_violin.R")




# file = '/Users/ulises/Desktop/ABL/software/qsin/test_data/n15_linear_tmp/compared_nets_n15_linear_tmp.csv'
# combs_file = '/Users/ulises/Desktop/ABL/software/qsin/test_data/n15_linear_tmp/linear_5e-4_overlappedBatches_uniq.txt'
# time_file ='/Users/ulises/Desktop/ABL/software/qsin/test_data/n15_linear_tmp/processed_time_n15_linear_tmp.csv'
# pseudo_file ='/Users/ulises/Desktop/ABL/software/qsin/test_data/n15_linear_tmp/processed_pseudo_n15_linear_tmp.csv'

# file = '/Users/ulises/Desktop/ABL/software/qsin/test_data/n15_depth3_v0.9_n0.5/compared_nets_n15_depth3_v0.9_n0.5.csv'
# combs_file = '/Users/ulises/Desktop/ABL/software/qsin/test_data/n15_depth3_v0.9_n0.5/nonlinear_depth3_nu0.9_overlappedBatches_uniq.txt'
# time_file ='/Users/ulises/Desktop/ABL/software/qsin/test_data/n15_depth3_v0.9_n0.5/processed_time_n15_depth3_v0.9_n0.5.csv'
# pseudo_file ='/Users/ulises/Desktop/ABL/software/qsin/test_data/n15_depth3_v0.9_n0.5/processed_pseudo_n15_depth3_v0.9_n0.5.csv'




file = '/Users/ulises/Desktop/ABL/software/qsin/test_data/n15_depth3_v0.9_n0.5/compared_nets_n15_depth3_v0.9_n0.5.csv'
combs_file = '/Users/ulises/Desktop/ABL/software/qsin/test_data/n15_depth3_v0.9_n0.5/nonlinear_depth3_nu0.9_overlappedBatches_uniq.txt'
time_file ='/Users/ulises/Desktop/ABL/software/qsin/test_data/n15_depth3_v0.9_n0.5/processed_time_n15_depth3_v0.9_n0.5.csv'
pseudo_file ='/Users/ulises/Desktop/ABL/software/qsin/test_data/n15_depth3_v0.9_n0.5/processed_pseudo_n15_depth3_v0.9_n0.5.csv'



is_linear = F
starts = 15

if (is_linear){
  msubtitle_d = paste('Elastic Net,', starts, 'random starts')
}else{
  msubtitle_d = paste('ISLE,', starts, 'random starts')
}



df_lik = read.csv(pseudo_file)%>%
  mutate(id = paste(row, boot, sep = ",")) %>%
  dplyr::select(id, likdev)


df_time = read.csv(time_file) %>%
  dplyr::rename( boot = bootstrap)%>%
  mutate(id = paste(row, boot, sep = ",")) %>%
  dplyr::select(id, time)


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


# 
df <- df[!is.na(df$dist),]

df %>%
  # filter(row > 1) %>%
  ggplot(data = ., aes(x = CT_row_count, y = dist,
                       colour = rel_likdev)) +
  geom_point(
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

# table(df[df$row == 4, 'dist'])

df %>%
  ggplot(data = ., 
         aes(x = CT_row_count,
             y = time/60/60,
             colour = rel_likdev)) +
  geom_point(position = position_jitter(width = 0.10), size = 1, alpha = 0.5, )+
  scale_colour_gradient(low = "blue", high = "red")+
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0, position = position_nudge(x = 0, y = 0)) +
  geom_flat_violin(width=0.5, lwd=0.5, position = position_nudge(x = 0.14, y = 0), alpha = 2) +
  theme_bw(base_size = 12) +
  labs(title="Running time",
       subtitle = msubtitle_d,
       y ="Hours", x = "Rows in CT")+
  theme(legend.position="none") +
  geom_hline(yintercept=max_time/60/60, linetype="dashed", color = "red") -> time_plot


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
# dev.off()