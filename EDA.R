# load packages
library(ggplot2)
source("PhD/Ogata.R")


# Load data
load("PhD-oversize/cleaned_whale_data.Rdata")
wh <- 53

incr <- cleaned_data[[wh]]$incr
ts <- incr * seq(1, length(cleaned_data[[wh]]$depths))
depths <- cleaned_data[[wh]]$depths
cues <- cleaned_data[[wh]]$times
cue_depths <- numeric(length(cues))

for(i in 1:length(cues)){
  cue_depths[i] = depths[floor(cues[i]/incr)]
}


# FIGURE 2: depths and cue times
depths <- data.frame(ts=ts,depths=depths)
cues <- data.frame(cues=cues,depths=cue_depths)

plot <- ggplot() +
  geom_line(data=depths, aes(x=(1/3600)*ts,y=-100*depths), linewidth=1.5) +
  geom_point(data=cues, aes(x=(1/3600)*cues,y=-100*depths),col="red") +
  ylab("Depth (m)") +
  xlab("Time (hr)") +
  theme_classic() +
  theme(axis.title = element_text(size = 20))

# FIGURE 2A
plotZ <- ggplot() +
  geom_line(data=depths, aes(x=(1/3600)*ts,y=-100*depths), linewidth=1.5) +
  geom_point(data=cues, aes(x=(1/3600)*cues,y=-100*depths),col="red") +
  ylab("Depth (m)") +
  xlab("Time (hr)") +
  xlim(5.4,6.1) +
  theme_classic() +
  theme(axis.title = element_text(size = 20))

# FIGURE 2B
plotN1 <- ggplot() +
  geom_line(data=depths, aes(x=(1/3600)*ts,y=100*depths), linewidth=2) +
  scale_y_reverse() +
  ylab("") +
  xlab("") +
  xlim(5.4,6.2) +
  theme_classic() +
  theme(axis.title = element_text(size = 27),
        axis.text.y=element_text(size=15),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

plotN2 <- ggplot() +
  geom_point(data=cues, aes(x=(1/3600)*cues,y=100*depths),col="red",size=1.33) +
  scale_y_reverse() +
  ylab("Depth (m)") +
  xlab("Time (hr)") +
  xlim(5.4,6.2) +
  theme_classic() +
  theme(axis.title = element_text(size = 27),
        axis.text=element_text(size=15))

gridExtra::grid.arrange(plotN1, plotN2, nrow=2)

# FIGURE 2B
# times as a point process with 1D band (depth) underneath

# times
plotB1 <- ggplot() +
  #geom_point(data=cues, aes(x=(1/3600)*cues,y=1.5),col="red",shape="I") +
  geom_segment(data=cues,aes(x=(1/3600)*cues,y=1.6,yend=1.4,xend=(1/3600)*cues),linewidth=0.01) +
  geom_rect(data=depths, 
            aes(xmin = (1/3600)*ts, xmax = (1/3600)*(ts+1),
                ymin = 0, ymax = 1, fill = depths*100)) +
  scale_fill_viridis_c() +
  labs(fill="Depth (m)") +
  xlab("Time (hr)") +
  #xlim(5.4,6.2) +
  #ylim(0,2) +
  coord_cartesian(xlim=c(5.4,6.2),ylim=c(0,2)) +
  #scale_x_continuous(limits=c(5.4,6.2)) +
  theme_classic() +
  theme(axis.title.x = element_text(size = 27),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size=15),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())


# depth
# use geom_rect
plotB2 <- ggplot() +
  geom_rect(data=depths, 
            aes(xmin = (1/3600)*ts, xmax = (1/3600)*(ts+1),
                ymin = 0, ymax = 1, fill = depths*100)) +
  scale_fill_viridis_c() +
  ylab("") +
  xlab("") +
  xlim(5.4,6.2) +
  theme_classic() +
  theme(axis.title = element_text(size = 27),
        axis.text.y=element_text(size=15),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

gridExtra::grid.arrange(plotB1, plotB2, nrow=2)

# FIGURE 4: IWP and Homogenous Poisson comparison
times = UWOgata(1, 0, 1, 1, 10000) # homogeneous Poisson process
diffs = times[2:length(times)] - times[1:(length(times)-1)]

# scatter plot of extract of Poisson times
data = data.frame(x = times, y = rep(1,length(times)))[100:119,]
library(ggplot2)
plot1 <- ggplot(data, aes(x,y)) +
  geom_point(size=10, col = "blue") +
  xlab("Time") +
  theme_classic() +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

times2 = UWOgata(1, 0, 1, 5, 10000) # homogeneous Poisson process
diffs2 = times2[2:length(times2)] - times2[1:(length(times2)-1)]

# scatter plot of extract of Weibull times
data2 = data.frame(x = times2, y = rep(1,length(times2)))[100:119,]
plot2 <- ggplot(data2, aes(x,y)) +
  geom_point(size=10, col = "orange") +
  xlab("Time") +
  theme_classic() +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

# histogram of Poisson interarrivals
data3 = data.frame(x = diffs)
plot3 <- ggplot(data3, aes(x)) +
  geom_histogram(color="blue", fill = "blue", boundary = 0, binwidth = 0.05) +
  xlim(0,4) +
  theme_classic() +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

# histogram of Weibull interarrivals
data4 = data.frame(x = diffs2)
plot4 <- ggplot(data4, aes(x)) +
  geom_histogram(color="orange", fill = "orange", boundary = 0, binwidth = 0.05) +
  xlim(0,4) +
  theme_classic() +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

library(ggplot2)
require(patchwork)
### npc function
library(grid)
require(gtable)

get_x_y_values <- function(gg_plot)
{
  img_dim      <- grDevices::dev.size("cm") * 10
  gt           <- ggplot2::ggplotGrob(gg_plot)
  to_mm        <- function(x) grid::convertUnit(x, "mm", valueOnly = TRUE)
  n_panel      <- which(gt$layout$name == "panel")
  panel_pos    <- gt$layout[n_panel, ]
  panel_kids   <- gtable::gtable_filter(gt, "panel")$grobs[[1]]$children
  point_grobs  <- panel_kids[[grep("point", names(panel_kids))]]
  from_top     <- sum(to_mm(gt$heights[seq(panel_pos$t - 1)]))
  from_left    <- sum(to_mm(gt$widths[seq(panel_pos$l - 1)]))
  from_right   <- sum(to_mm(gt$widths[-seq(panel_pos$l)]))
  from_bottom  <- sum(to_mm(gt$heights[-seq(panel_pos$t)]))
  panel_height <- img_dim[2] - from_top - from_bottom
  panel_width  <- img_dim[1] - from_left - from_right
  xvals        <- as.numeric(point_grobs$x)
  yvals        <- as.numeric(point_grobs$y)
  yvals        <- yvals * panel_height + from_bottom
  xvals        <- xvals * panel_width + from_left
  data.frame(x = xvals/img_dim[1], y = yvals/img_dim[2])
}
## simulate pp
set.seed(1234)
########################
## Weibull waiting plot
########################
wei <- ggplot(data.frame(x = seq(0,10, length.out = 100),
                         y = dweibull(seq(0,10, length.out = 100),
                                      shape = 5, scale = 1/gamma(1.2))),
              aes(x = x, y = y)) +
  theme_void() + geom_polygon(fill = "grey", alpha = 0.2)
## base plot Weibull
ts <- (times2[100:109] - times2[100]) * 10/(times2[109] - times2[100])
base <- ggplot(data.frame(x = ts, y = rep(0, 10)), aes(x = x, y = y)) +
  geom_point(size = 5) +  theme_classic() +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) + ylab("") + xlab("") +
  xlim(0, 10)
xy <- get_x_y_values(base)
## add Weibull waiting time plots
for(i in 1:nrow(xy)){
  base <- base + inset_element(p = wei, left = xy$x[i],
                               bottom =  xy$y[1], right = xy$x[i] + 0.2,
                               top = 0.75, align_to = "plot")
}
print(base)
########################
## Poisson waiting plot
########################
poi <- ggplot(data.frame(x = seq(0,10, length.out = 100),
                         y = dexp(seq(0,10, length.out = 100),
                                      mean = 1),
              aes(x = x, y = y))) + 
  theme_void() + geom_col(fill = "grey", alpha = 0.2, col = "darkgrey")
## base plot Poisson
base <- ggplot(data.frame(x = times[100:109], y = rep(0, 10)), aes(x = x, y = y)) +
  geom_point(size = 5) +  theme_classic() +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) + ylab("") + xlab("") +
  xlim(0, 10)
xy <- get_x_y_values(base)
## add Poisson waiting time plots
for(i in 1:nrow(xy)){
  base <- base + inset_element(p = poi, left = xy$x[i],
                               bottom =  xy$y[1], right = xy$x[i] + 0.1,
                               top = 0.75, align_to = "plot")
}
print(base)


library(ggpubr)
out <- ggarrange(plot1, 
                 ggarrange(plot3, plot4), 
                 plot2, nrow = 3)


# FIGURE 1: ICI histogram
times = cleaned_data[[wh]]$times; l = length(times)
IAT = times[2:l] - times[1:(l-1)]

data <- data.frame(IAT = IAT)
plot <- ggplot(data = data, aes(x=IAT)) +
  geom_histogram(boundary=0, binwidth=0.01, aes(y = after_stat(density))) +
  stat_function(fun = dexp, args = (mean = 0.774695), color = "red", linewidth = 2) +
  xlim(0,2) +
  xlab("Interclick Interval (s)") +
  theme_minimal() +
  theme(axis.title = element_text(size = 20))
