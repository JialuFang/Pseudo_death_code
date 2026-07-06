library(tidyr)
library(ggplot2)
library(gridExtra)

t_values <- seq(0, 2, by = 0.01)

custom_labels <- c("S_D1" = bquote(italic(S[D1])), 
                   "S_D2" = bquote(italic(S[D0])), 
                   "S_R1" = bquote(italic(S[R1])), 
                   "S_R2" = bquote(italic(S[R0])))
colors <- c("S_D1" = "blue", "S_D2" = "green", "S_R1" = "red", "S_R2" = "purple")

t0 = 1; lambda = 1/4

df_list <- list(
  data.frame(
    t = t_values,
    S_D1 = 1/(1+t_values/2),
    S_D2 = 1/(1+t_values/2)/(1+pmin(t_values, 1+t0/2)-pmin(t_values,1-t0/2)), 
    S_R1 = exp(-3/2*t_values),
    S_R2 = exp(-(1+lambda)*3/2 * t_values) #lambda=1/4
  ),
  data.frame(
    t = t_values,
    S_D1 = 1/(1+t_values/2),
    S_D2 = 1/(1+t_values/2)/(1+pmin(t_values, 1+t0/2)-pmin(t_values,1-t0/2)),
    S_R1 = exp(-3/2*t_values),
    S_R2 = exp(-(3/2+lambda)*t_values + 2*lambda*pmax(t_values-1, 0))
  ),
  data.frame(
    t = t_values,
    S_D1 = 1/(1+t_values/2),
    S_D2 = 1/(1+t_values/2)/(1+pmin(t_values, 1+t0/2)-pmin(t_values,1-t0/2)),
    S_R1 = exp(-3/2*t_values),
    S_R2 = (1/2)*(1/(1+3*t_values))+(1/2)*exp(-3/2*t_values)
  ),
  data.frame(
    t = t_values,
    S_D1 = 1/(1+t_values/2),
    S_D2 = 1/(1+t_values/2)/(1+pmin(t_values, 1+t0/2)-pmin(t_values,1-t0/2)), 
    S_R1 = 1/(1+2*t_values),
    S_R2 = 2/(1+2*t_values)/(2+pmin(t_values, t0))
  ),
  data.frame(
    t = t_values,
    S_D1 = 1/(1+t_values/2),
    S_D2 = 1/(1+t_values/2)/(1+pmin(t_values, 1+t0/2)-pmin(t_values,1-t0/2)), 
    S_R1 = 1/(1+2*t_values),
    S_R2 = 2/(1+2*t_values)/(2+pmin(t_values, 1+t0/2)-pmin(t_values,1-t0/2))
  ),
  data.frame(
    t = t_values,
    S_D1 = 1/(1+t_values/2),
    S_D2 = 1/(1+t_values/2)/(1+pmin(t_values, 1+t0/2)-pmin(t_values,1-t0/2)), 
    S_R1 = 1/(1+2*t_values),
    S_R2 = 2/(1+2*t_values)/(2+pmax(t_values+t0-2, 0))
  ),
  data.frame(
    t = t_values,
    S_D1 = 1/(1+t_values/2),
    S_D2 = 1/(1+t_values/2)/(1+pmin(t_values, 1+t0/2)-pmin(t_values,1-t0/2)), 
    S_R1 = exp(-(1+lambda)*3/2 * t_values), #lambda=1/4
    S_R2 = exp(-3/2*t_values)
  )
)

df_list_three <- list(
  data.frame(
    t = t_values,
    S_1 = 1/(1+t_values/2),
    S_2 = 1/(1+t_values/2)/(1+pmin(t_values, 1+t0/2)-pmin(t_values,1-t0/2)),
    S_3 = exp(-3/2*t_values)
  )
)

custom_labels_three <- list(
  c("S_1" = bquote(italic(S[D1])), 
    "S_2" = bquote(italic(S[D0])), 
    "S_3" = bquote(italic(S[R1]==S[R0]))),
  c("S_1" = bquote(italic(S[D1]==S[D0])), 
    "S_2" = bquote(italic(S[R1])), 
    "S_3" = bquote(italic(S[R0])))  
)

colors_three <- c("S_1" = "blue", "S_2" = "green", "S_3" = "red")
title<- c("Proportional hazards", "Crossing hazards", "Crossing survival functions", 
          "Early difference", "Middle difference", "Late difference",
          "Counteracting effect")

plot_list <- list()

pdf(file = "plots/simu_setting_supp.pdf", width = 8, height = 5.5)

for (i in 1:length(df_list)) {
  df_long <- pivot_longer(df_list[[i]], cols = -t, names_to = "Function", values_to = "Value")
  
  p <- ggplot(df_long, aes(x = t, y = Value, color = Function)) +
    geom_line(linewidth = 0.5) + 
    scale_color_manual(values = colors, labels = custom_labels) +
    labs(x = "t",
         y = "S(t)") + theme_minimal() +
    theme(legend.key.size = unit(0.5, "cm"),
          plot.title = element_text(size = 10,hjust=0.5))+
    guides(color = guide_legend(title = NULL)) + 
    ggtitle(title[i])
  
  plot_list[[i]] <- p
}

for (i in 1:1) {
  df_long <- pivot_longer(df_list_three[[i]], cols = -t, names_to = "Function", values_to = "Value")
  
  p <- ggplot(df_long, aes(x = t, y = Value, color = Function)) +
    geom_line(linewidth = 0.5) + 
    scale_color_manual(values = colors_three, labels = custom_labels_three[[i]]) +
    labs(x = "t",
         y = "S(t)") + theme_minimal() +
    theme(legend.key.size = unit(0.5, "cm"),
          legend.text = element_text(size = 7),
          legend.margin = margin(t = 0, r = 0, b = 0, l = 0),
          plot.title = element_text(size = 10,hjust=0.5))+
    guides(color = guide_legend(title = NULL)) + 
    ggtitle("Absent effect")
  
  plot_list[[length(df_list)+i]] <- p
}

grid.arrange(grobs = plot_list, ncol = 3)

dev.off()
