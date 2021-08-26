library(readxl)
library(readr)
library(reshape2)
library(dplyr)
library(ggplot2)
library(cowplot)
library(gplots)
library(ggpubr)
library(tidyverse)
options(scipen = 999)
options(readr.default_locale=readr::locale(tz="US/Eastern"))
#-----------
# Figure 2
#----------
# modifications
ggdensity(subset(clinicalML.results.long, variable %in% c('AUROC', 'Accuracy')) 
          , "value", fill = "variable", palette = "jco") 
# spider chart
# Color vector
colors_border=c( rgb(0.2,0.5,0.5,0.9), rgb(0.8,0.2,0.5,0.9) , rgb(0.7,0.5,0.1,0.9) )
colors_in=c( rgb(0.2,0.5,0.5,0.4), rgb(0.8,0.2,0.5,0.4) , rgb(0.7,0.5,0.1,0.4) )
# plot with default options:
radarchart( clinicalML.results.wide[c(3:5),c(7,9,12,13)]  , axistype=0 , maxmin = F,
            #custom polygon
            pcol=colors_border , pfcol=colors_in , plwd=4 , plty=1,
            #custom the grid
            cglcol="grey", cglty=1, axislabcol="grey", caxislabels=seq(0,20,5), cglwd=0.8,
            #custom labels
            vlcex=0.8 
)
legend(x=0.7, y=1, legend = c('LR', 'RF', 'XGB'), bty = "n", pch=20 , col=colors_in , text.col = "grey", cex=1.2, pt.cex=3)
# fig.2: bar plots
fig.2 <- subset(clinicalML.results.wide, Folds == 5) %>%
  select(c(Model, AUROC, Accuracy, PPV, NPV, `95_lower`, `95_upper`)) %>%
  mutate(ID = 1:n()) %>%
  gather(key, value, c(2:5)) %>%
  mutate(`95_upper` = ifelse(key != 'AUROC', NA, `95_upper`)) %>%
  mutate(`95_lower` = ifelse(key != 'AUROC', NA, `95_lower`)) %>%
  ggplot(aes(x = Model, y = value, fill = Model)) +
  geom_bar(stat="identity", color = 'black', width = 1, position ='dodge')+
  #geom_errorbar(aes(ymin = `95_lower`, ymax = `95_upper`), color = 'navyblue', width = .5) +
  geom_text(aes(label= format(round(value,digits = 2), nsmall = 2) #round(value,2)
                ), 
            vjust=1.6, color="black", size = 3.5, family = 'serif' ) +
  facet_wrap(key ~.) +
  scale_fill_manual(values = c( 'plum4', 'gray','slategray4')) +
  scale_y_continuous(limits = c(0, 1)) +
  theme_light() +
  theme(
    legend.position = 'none'
    , axis.title = element_text(size = 10, color = 'black', family = 'serif')
    , strip.text.x = element_text(size = 10, color = 'black', family = 'serif')
    , axis.text = element_text(size = 10, color = 'black', family = 'serif')
  ) +
  labs(y = 'Performance Measure', x = 'Model') +
  coord_cartesian(clip = 'off')
fig.3.heatmap <- ggplot(subset(clinicalML.feaimp.long, variable != 'Average'), aes(x = variable, y = Var.abbr))+
  geom_tile(aes(fill = value)) +
  theme_classic()+
  theme(
    legend.position = "top"
    #legend.position = c(.5,.5) , legend.direction = 'vertical'
    , axis.title.y = element_blank()
    , axis.title.x = element_text(size = 10, color = 'black', family = 'serif')
    , axis.text.y = element_text(size = 10, color = 'black', family = 'serif')
    , axis.text.x = element_text(size = 10, color = 'black', family = 'serif')
    , legend.title = element_text(size = 10, color = 'black', family = 'serif')
    , legend.text = element_text(size = 10, color = 'black', family = 'serif')
  )+
  scale_fill_gradient(name = "Feature Importance",
                      low = "lightcoral",high = "steelblue",
                      breaks = c(0, 50, 100),
                      labels = c('Low', 50, 'High'))+
  coord_cartesian(clip = 'off') +
  labs(x = 'Model')
fig.3.barplot <- ggplot(subset(clinicalML.feaimp.long, variable %in% c('Average')),
                        aes(y = Var.abbr, x =value))+
  geom_bar(stat = 'identity', position = 'stack', alpha = .8, fill = 'steelblue') +
  theme(
    axis.title.y = element_blank()
    , axis.title.x = element_text(size = 10, color = 'black', family = 'serif')
    , axis.text.y = element_blank()
    , axis.text.x = element_text(size = 10, color = 'black', family = 'serif')
    , axis.ticks.y = element_blank()
    , panel.background = element_blank()
    , legend.position = 'none'
    #, plot.margin = unit(c(2.15, 0, 0, 0), 'cm') #t,r,b,l
    , plot.margin = unit(c(2.0, 0, 0.20, 0), 'cm') #t,r,b,l
  ) + 
  geom_text(aes(label = format(round(value,digits = 2), nsmall = 2)#round(value,1)
                )
            , size = 3, nudge_y = 0, nudge_x = 3.25, color = "black", family = 'serif') +
  coord_cartesian(clip = 'off', xlim = c(0, 100)) +
  labs(x = 'Average feature importance score')
tiff("figure_2.tiff", units = "in", width = 18, height = 10, res = 300, compression = 'jpeg')
plot_grid(
  fig.2, fig.3.heatmap, fig.3.barplot,
  nrow = 1 , rel_widths = c(.4, .35, .25) ,
  labels = c("(A)", "(B)", "(C)"), label_size = 12, label_fontfamily = "serif",
  label_y = .99
)
dev.off()
#-----------
# Figure 4
#----------
#modified plots
GGally::ggparcoord(NLP.results.wide%>%
                     arrange(Model),
                   columns = 5:8, groupColumn = 4, order = "anyClass",
                   showPoints = TRUE, 
                   title = "Parallel Coordinate Plot for the Iris Data",
                   alphaLines = 0.3
) + 
  theme(
    plot.title = element_text(size=10)
  )
# using ggplot
NLP.results.wide %>%
  mutate(ID = 1:n()) %>%
  #mutate_if(is.numeric, scale) %>%
  gather(key, value, c(5:8)) %>%
  ggplot(aes(key, value, group = ID, colour = Model, fill = Category)) +
  geom_line() +
  geom_point(size=2, shape=21, colour="grey50") +
  scale_color_manual(values = c('grey', 'yellow', 'blue', 'red'))
# radar chart
coord_radar <- function (theta = "x", start = 0, direction = -1) {
  theta <- match.arg(theta, c("x", "y"))
  r <- if (theta == "x") "y" else "x"
  ggproto("CordRadar", CoordPolar, theta = theta, r = r, start = start, 
          direction = sign(direction),
          is_linear = function(coord) TRUE)
}
fig4.triage.precovid <- subset(NLP.results.long, Category == 'triage' & COVID == 'pre' & Splits == '80:20' & `What Controls?` == 'Controls1+Controls3') %>% 
  mutate_if(is.numeric, ~replace(.,is.na(.),0)) %>%
  ggplot(aes(x = variable, y = value, group = `No.of cuis`, fill = `No.of cuis`)) +
  geom_polygon(alpha = .25) +
  facet_wrap(Model ~.) +
  coord_radar() +
  labs(title = 'Triage PRE-COVID') 
model_names <- c(`GLM` = 'GLM', `RF` = 'RF', `SVM` = 'SVM', `XGB` = 'XGB')
notes_names <- c(`provider` = 'Provider Notes', `triage` = 'Triage Notes')
covid_names <- c(`post` = 'Post-COVID', `pre` = 'Pre-COVID')
# ALL models 
 coord_radar_test <- function (theta = "x", start = 0, direction = 1) {
  theta <- match.arg(theta, c("x", "y"))
  r <- if (theta == "x") "y" else "x"
  #dirty
  rename_data <- function(coord, data) {
    if (coord$theta == "y") {
      plyr::rename(data, c("y" = "theta", "x" = "r"), warn_missing = FALSE)
    } else {
      plyr::rename(data, c("y" = "r", "x" = "theta"), warn_missing = FALSE)
    }
  }
  theta_rescale <- function(coord, x, scale_details) {
    rotate <- function(x) (x + coord$start) %% (2 * pi) * coord$direction
    rotate(scales::rescale(x, c(0, 2 * pi), scale_details$theta.range))
  }
  r_rescale <- function(coord, x, scale_details) {
    scales::rescale(x, c(0, 0.4), scale_details$r.range)
  }
  ggproto("CordRadar", CoordPolar, theta = theta, r = r, start = start, 
          direction = sign(direction),
          is_linear = function(coord) TRUE,
          render_bg = function(self, scale_details, theme) {
            scale_details <- rename_data(self, scale_details)

            theta <- if (length(scale_details$theta.major) > 0)
              theta_rescale(self, scale_details$theta.major, scale_details)
            thetamin <- if (length(scale_details$theta.minor) > 0)
              theta_rescale(self, scale_details$theta.minor, scale_details)
            thetafine <- seq(0, 2 * pi, length.out = 100)

            rfine <- c(r_rescale(self, scale_details$r.major, scale_details))

            # This gets the proper theme element for theta and r grid lines:
            #   panel.grid.major.x or .y
            majortheta <- paste("panel.grid.major.", self$theta, sep = "")
            minortheta <- paste("panel.grid.minor.", self$theta, sep = "")
            majorr     <- paste("panel.grid.major.", self$r,     sep = "")

            ggplot2:::ggname("grill", grid::grobTree(
              ggplot2:::element_render(theme, "panel.background"),
              if (length(theta) > 0) ggplot2:::element_render(
                theme, majortheta, name = "angle",
                x = c(rbind(0, 0.45 * sin(theta))) + 0.5,
                y = c(rbind(0, 0.45 * cos(theta))) + 0.5,
                id.lengths = rep(2, length(theta)),
                default.units = "native"
              ),
              if (length(thetamin) > 0) ggplot2:::element_render(
                theme, minortheta, name = "angle",
                x = c(rbind(0, 0.45 * sin(thetamin))) + 0.5,
                y = c(rbind(0, 0.45 * cos(thetamin))) + 0.5,
                id.lengths = rep(2, length(thetamin)),
                default.units = "native"
              ),

              ggplot2:::element_render(
                theme, majorr, name = "radius",
                x = rep(rfine, each = length(thetafine)) * sin(thetafine) + 0.5,
                y = rep(rfine, each = length(thetafine)) * cos(thetafine) + 0.5,
                id.lengths = rep(length(thetafine), length(rfine)),
                default.units = "native"
              )
            ))
          })
}
tiff("figure_4_supp.tiff", units = "in", width = 8, height = 10, res = 300, compression = 'jpeg')
NLP.results.long$`Percentage of All CUIs used:` <- if_else(
  NLP.results.long$`No.of cuis` == 31, '5%', if_else(
    NLP.results.long$`No.of cuis` == 62, '10%', if_else(
      NLP.results.long$`No.of cuis` == 124, '20%', if_else(
        NLP.results.long$`No.of cuis` %in% c(291, 332), '50%', '100%'
      )
    )
  )
)
NLP.results.long$`Percentage of All CUIs used:` <- factor(NLP.results.long$`Percentage of All CUIs used:`, 
                                                          levels = c('5%', '10%', '20%', '50%', '100%'))
NLP.results.long %>%
  ggplot(aes(x = variable, y = value)) +
  geom_polygon(aes(group = `No.of cuis`, 
                   color = `Percentage of All CUIs used:`),
               size = .5, alpha = .15) +
  facet_wrap(Model ~ Category + COVID, 
             labeller = label_wrap_gen(multi_line=T)) +
  coord_radar_test(start = -pi) +
  theme_light() +
  theme(
    axis.title.x = element_blank() 
    , panel.spacing = unit(c(-0.5,0-0.5,0), "lines")
    , strip.text = element_text(size = 10, color = 'black', family = 'serif')
    , legend.text = element_text(size = 10, color = 'black', family = 'serif')
    , legend.title = element_text(size = 10, color = 'black', family = 'serif')
    , axis.title = element_text(size = 10, color = 'black', family = 'serif')
    , legend.position = 'top'
  ) +
  labs(y = 'Performance Measure' , fill = 'Percentage of All CUIs used:') +
  scale_y_continuous(limits = c(0,1), breaks = seq(0,1,.25), expand = c(0,0)) +
  scale_fill_gradient2(limits = c(0,1), midpoint = .5)
dev.off()
# Best performers
# Best performers
melt(subset(NLP.results.wide, AUROC >=.8 & Accuracy >= .8 & PPV >= .8 & NPV >= .8)) %>%
  ggplot(aes(x = variable, y = value, group = `No.of cuis`, fill = `No.of cuis`)) +
  geom_polygon(alpha = .25, size = 1) +
  coord_polar(start = -pi) +
  facet_wrap(Category+ Model ~ COVID,
              labeller = labeller(
                Model = as_labeller(c(`rf` = 'RF', `xgbDART` = 'XGB')) ,
                Category = as_labeller(c(`provider` = 'Provider Notes')) ,
                COVID = as_labeller(c(`pre` = 'Pre-COVID', `post` = 'Post-COVID')),
                groupwrap = label_wrap_gen(width = 5))
  ) +
  theme_light() +
  theme(
    axis.title.x = element_blank() 
    , strip.text = element_text(face = 'bold')
    , legend.position = 'top'
  ) +
  labs(y = 'Performance Measure' , fill = 'No.of CUIS') +
  guides(fill = guide_legend(nrow = 1))

tmp <- read_excel("~/Genentech/tmp/figuresData.xlsx"
           , sheet = 'NLP_v3Results_04202021'
           , na=c("NA")
           , col_names = TRUE) %>%
  #filter(`What Controls?` %in% 'Controls1+Controls3' & Splits == '80:20' & AUROC >=.9 & Accuracy >=.9 & Sensitivity >= .9 & Specificity >= .9) %>%
  filter(`What Controls?` %in% 'Controls1+Controls3' & Splits == '80:20' & AUROC >=.95) %>%
  select(c(Category, COVID, desc, `No.of cuis`, Model, AUROC, Accuracy, PPV, NPV)) %>%
  mutate_at(vars(Category, COVID, desc, `No.of cuis`, Model), list(factor))
tmp$COVID <- factor(tmp$COVID, levels = c('pre', 'post'))
melt(tmp) %>%
  ggplot(aes(x = variable, y = value, group = `No.of cuis`, fill = `No.of cuis`)) +
  geom_polygon(alpha = .25, size = 1) +
  coord_polar(start = -pi) +
  facet_wrap(Category ~ COVID + Model, nrow = 2,
             labeller = labeller(
               Model = as_labeller(c(`rf` = 'RF', `svmRadial` = 'SVM',`xgbDART` = 'XGB')) ,
               Category = as_labeller(c(`provider` = 'Provider Notes', `triage` = 'Triage Notes')) ,
               COVID = as_labeller(c(`pre` = 'Pre-COVID', `post` = 'Post-COVID')),
               groupwrap = label_wrap_gen(width = 5))
  ) +
  theme_light() +
  theme(
    axis.title.x = element_blank() 
    , strip.text = element_text(face = 'bold')
    , legend.position = 'top'
  ) +
  labs(y = 'Performance Measure' , fill = 'No.of CUIS') +
  guides(fill = guide_legend(nrow = 1))
tiff("figure_4_main.tiff", units = "in", width = 8, height = 8, res = 300, compression = 'jpeg')
#fig.4.main <- 
  tmp %>%
  mutate(ID = 1:n()) %>%
  gather(key, value, c(6:9)) %>%
  ggplot(aes(key, value, group = ID, colour = Model#, fill = `No.of cuis`
             )) +
  geom_line(aes(linetype = Model),  size = .5) +
  geom_point(#aes(shape = `No.of cuis`),
             size = 1.5, colour = "black", alpha = .35) +
  facet_wrap(Category ~ COVID,
             labeller = labeller(
               Category = as_labeller(c(`provider` = 'Provider Notes', `triage` = 'Triage Notes')) ,
               COVID = as_labeller(c(`pre` = 'Pre-COVID', `post` = 'Post-COVID')),
               groupwrap = label_wrap_gen(width = 5)
             )
  ) +
  scale_linetype_manual(values = c('dotted' ,'solid', 'twodash')) +
  scale_color_manual(values = c('black', 'blue', 'green')) +
  theme_light() +
  theme(
    legend.position = 'top'
    , axis.title.x = element_blank()
    , axis.title.y =  element_text(size = 10, color = 'black', family = 'serif')
    , strip.text.x = element_text(size = 10, color = 'black', family = 'serif')
    , axis.text = element_text(size = 10, color = 'black', family = 'serif')
  ) +
  labs(y = 'Performance Measure', shape = 'No.of CUIs', fill = 'No.of CUIs') +
  coord_cartesian(clip = 'off')
dev.off()
#----------
# Figure 5a
#----------
tmp <- TRIAGE.binary %>%
  group_by(label) %>%
  sample_n(25)
tmp1 <- melt(tmp) %>%
  #melt(TRIAGE.binary) %>%
  mutate(
    labelfactor = cut(value, breaks = c(-2, -1, -.5, 0, .5, 1)
                      ,labels = c('-1', '-.5', '0', '.5', '1'))
  ) %>%
  mutate(labelfactor = factor(as.character(labelfactor), 
                              levels = rev(levels(labelfactor))))
tmp1$label <- if_else(
  tmp1$label %in% c('CONTROLS1'), 'CONTROLS-1', if_else(
    tmp1$label %in% c('CONTROLS3'), 'CONTROLS-3', if_else(
      tmp1$label %in% c('v.cases'), 'TEST-CASES', if_else(
        tmp1$label %in% c('v.controls'), 'TEST-CONTROLS', 'CASES'
      )
    )
  )
)
tmp2 <- unique(tmp1)
cui.to.text <- read_excel('fig.5a.xlsx') 
tmp3 <- left_join(tmp2, cui.to.text, by = c('variable','label', 'value'))
tmp3$displayText <- left_join(tmp2, subset(cui.to.text, value==1), by = c('variable','label', 'value'))$cuiText
tmp3$variable <- factor(tmp3$variable, levels(tmp2$variable))
tmp3$label <- if_else(
  tmp3$label %in% c('CONTROLS-1'), 'CONTROLS-1', if_else(
    tmp3$label %in% c('CONTROLS-3'), 'CONTROLS-3', if_else(
      tmp3$label %in% c('TEST-CASES'), 'TEST-CASES', if_else(
        tmp3$label %in% c('TEST-CONTROLS'), 'TEST-CONTROLS', 'CASES'
      )
    )
  )
)
x.axis.breaks <- unique(tmp3$positiveText)[-1]
tiff("fig.5a.test.tiff", units = "in", width = 20, height = 20, res = 500, compression = 'jpeg')
fig.5a <- ggplot(tmp3, aes(y = variable, x = label, fill = value))+
  geom_tile(colour = "white", size = 0.5)+
  guides(fill = guide_legend(title = "CUI strength:"))+
  labs(x = "", y = "Concept Unique Identifiers (CUIs)")+
  theme_grey(base_size = 10)+
  theme(
    legend.position = "top", legend.direction = "horizontal"
    , legend.title = element_text(size = 20, color = 'black', family = 'serif')
    , legend.margin = margin(grid::unit(0,"cm"))
    , legend.text = element_text(size = 22, color = 'black', family = 'serif')
    , axis.text.x = element_text(size = 18,colour = 'black', vjust = .5, angle = 0, family = 'serif')
    , axis.text.y = element_blank()
    , axis.ticks = element_blank()
    , axis.title = element_text(size = 24, color = 'black', family = 'serif')
    , plot.background = element_blank()
    , panel.border = element_blank()
    , plot.margin = margin(0.7, 0.4, 0.1, 0.2, "cm")
  )  +
  scale_x_discrete(labels = c(
    'CASES' = 'Development Cohort\nCASES',
    'CONTROLS-1' = 'Development Cohort\nCONTROLS-1',
    'CONTROLS-3' = 'Development Cohort\nCONTROLS-2',
    'TEST-CASES' = 'Validation Cohort\nCASES',
    'TEST-CONTROLS' = 'Validation Cohort\nCONTROLS'
  )) +
  scale_fill_gradient(name = "CUI strength",
                      low = "lightcoral",high = "navyblue",
                      breaks = c(-1, 0, 1),
                      labels = c('-1', '0', '1')) +
  coord_cartesian(expand = FALSE)
dev.off()
#-----------
# Figure 5b
#----------
tmp <- NLP.feaimp.long %>%
  select(-c(displayonly)) %>%
  group_by(model, short_category) %>%
  arrange(model, short_category, desc(score)) %>%
  filter(score > quantile(score, .85))
fig.5b.2 <- subset(tmp, short_category =='PROVIDER') %>%
  ggplot(aes(y = text, x = model, fill = score)) +
  geom_tile() +
  #coord_fixed(ratio = .1) +
  theme(
    axis.text = element_text(hjust = .5, size = 10, color = 'black', family = 'serif')
    , axis.title = element_text(size = 10, color = 'black', family = 'serif')
    , axis.ticks = element_blank()
    , legend.position = 'top'
    , legend.title.align = .5
    , legend.title = element_text(size = 10, color = 'black', family = 'serif')
    , legend.text = element_text(size = 10, color = 'black', family = 'serif')
    , strip.text = element_text(size = 10, color = 'black', family = 'serif')
    , strip.background = element_rect(fill = 'gray')
    #, panel.spacing.x = unit(1, 'mm') 
    #, plot.margin = unit(c(0, 0, 0, 0), "cm")
    , plot.margin = margin(0,0,0,0) #t,r,b,l
  )  +
  scale_fill_gradient(name = "Feature Importance",
                      low = "lightcoral",high = "navyblue",
                      breaks = c(0, 50, 100),
                      labels = c('Low', 50, 'High'))+
  scale_y_discrete(position = 'right') +
  labs(x = 'Model', y = 'Text')
fig.5b.1 <- subset(tmp, short_category =='TRIAGE') %>%
  ggplot(aes(y = text, x = model, fill = score)) +
  geom_tile() +
  #coord_fixed(ratio = .1) +
  theme(
    axis.text = element_text(hjust = .5, size = 10, color = 'black', family = 'serif')
    , axis.title = element_text(size = 10, color = 'black', family = 'serif')
    , axis.ticks = element_blank()
    , legend.position = 'top'
    , legend.title.align = .5
    , legend.title = element_text(size = 10, color = 'black', family = 'serif')
    , legend.text = element_text(size = 10, color = 'black', family = 'serif')
    , strip.text = element_text(size = 10, color = 'black', family = 'serif')
    , strip.background = element_rect(fill = 'gray')
    #, panel.spacing.x = unit(1, 'mm') 
    #, plot.margin = unit(c(0, 0, 0, 0), "cm")
    , plot.margin = margin(0,0,0,0)
  )  +
  scale_fill_gradient(name = "Feature Importance",
                      low = "lightcoral",high = "navyblue",
                      breaks = c(0, 50, 100),
                      labels = c('Low', 50, 'High'))+
  scale_y_discrete(position = 'left') +
  labs(x = 'Model', y = 'Text')
plot_grid(fig.5b.1, fig.5b.2)
# averaging the score
tmp <- read_excel('figuresData.xlsx', sheet = 'displayonly') %>%
  select(c('model', 'short_category', 'formatted_text', 'score')) %>%
  group_by(short_category, formatted_text) %>%
  summarise(score = mean(score)) %>%
  group_by(short_category) %>%
  arrange(short_category, desc(score)) %>%
  filter(score > quantile(score, .93))

fig.5b.1 <- subset(tmp, short_category == 'TRIAGE') %>%
  ggplot(aes(y = reorder(formatted_text, -score), 
             score)) +
  geom_col(fill = 'lightskyblue1') +
  geom_text(aes(label = format(round(score,digits = 2), nsmall = 2)),
            vjust = .5, size = 6, hjust = 1.0) +
  geom_vline(xintercept = c(10, 20, 30, 40, 50), linetype = 3, color = 'dimgray') +
  labs(x ='', y = 'Concept Unique Identifiers (Triage Notes)') +
  theme_classic() %+replace%
  theme(
      axis.text = element_text(hjust = .5, size = 20, color = 'black', family = 'serif')
      #, axis.text.x = element_text(size = 20,colour = 'black', vjust = .5, angle = 0, family = 'serif')
    , axis.title = element_text(size = 24, color = 'black', family = 'serif')
    , axis.ticks = element_blank()
    , plot.margin = margin(50,0,0,0)
  ) +
  scale_x_continuous(breaks = seq(0, 55, by = 5),
                     labels = seq(0, 55, by = 5),
                     expand = expansion(mult = c(0,0.075)), limits = c(0,NA)) 

fig.5b.2 <-  subset(tmp, short_category == 'PROVIDER') %>%
  ggplot(aes(y = reorder(formatted_text, -score), 
             -score)) +
  geom_col(fill = 'lightskyblue3') +
  geom_text(aes(label = format(round(score,digits = 2), nsmall = 2)), 
            vjust = .5, size = 6, hjust = -.15) + #.01
  geom_vline(xintercept = c(-10, -20, -30, -40, -50), linetype = 3, color = 'dimgray') +
  labs(x='', y = 'Concept Unique Identifiers (Provider Notes)') +
  theme_classic() %+replace%
  theme(
      axis.text = element_text(hjust = .5, size = 20, color = 'black', family = 'serif')
    , axis.title = element_text(size = 24, color = 'black', family = 'serif')
    , axis.ticks = element_blank()
    , plot.margin = margin(50,0,0,0)
  ) +
  scale_y_discrete(position = 'right') +
  scale_x_continuous(breaks = seq(0, -55, by = -5),
                     labels = seq(0, 55, by = 5)
                     , expand = expansion(mult = c(0.1, 0)), limits = c(NA,0)
  ) 

tiff("figure_5b.tiff", units = "in", width = 10, height = 8, res = 300, compression = 'jpeg')
fig.5b <- ggdraw(add_sub(
  plot_grid(fig.5b.1, fig.5b.2),
  "Average Feature Importance Score",                       # to make single axis title
  vpadding=grid::unit(0,"lines"),y=6, x=0.5, vjust= 5.5, hjust = .5,
  fontfamily = 'serif', fontface = 'plain', color = 'black', size =10))
dev.off()
#----------
# Figure 5: combined 5a and 5b
#----------
tiff("fig.5.tiff", units = "in", width = 35, height = 20, res = 500, compression = 'jpeg')
plot_grid(fig.5a , NULL,
          ggdraw(add_sub(
            plot_grid(fig.5b.1, fig.5b.2),
            "Average Feature Importance Score",                       # to make single axis title
            vpadding=grid::unit(0,"lines"),y=6, x=0.5, vjust= 5.5, hjust = .5,
            fontfamily = 'serif', fontface = 'plain', color = 'black', size =24)
            )
          , rel_widths = c(.75, .05, 1.20)
          , nrow = 1
          , labels = c("A", "", "B")
          , label_size = 26
          )
dev.off()
