library(flexsurv)
library(ggplot2)
library(survminer)
library(magrittr)
library(dplyr)
library(gifski)
library(gganimate)

save_path = "outputs/"

# will use the ovarian dataset from flexsurv
data("ovarian")
ova_surv = survfit(Surv(futime, fustat) ~ 1, data = ovarian)
plot(ova_surv)

# cut the follow-up time so that data are less mature
ova_short =
  ovarian %>%
  mutate(
    fustat = case_when(
      futime > 400 ~ 0,
      TRUE ~ fustat
    ),
    futime = case_when(
      futime > 400 ~ 400,
      TRUE ~ futime
    )
  )

# set maximum survival time for the extrapolations
tmax = 2000
ova_surv = survfit(Surv(futime, fustat) ~ 1, data = ova_short)

# weibull and loglogistic parametric model fits
ova_wei = flexsurvreg(Surv(futime, fustat) ~ 1, data = ova_short, dist = "weibull")
ova_llo = flexsurvreg(Surv(futime, fustat) ~ 1, data = ova_short, dist = "llogis")
plot(ova_wei, t = seq(0, tmax), xlim = c(0, tmax), col = "red")
plot(ova_llo, t = seq(0, tmax), add = T, col = "blue")

### First plot: predictive uncertainty

# build ancillary dataset for the extrapolated curves
df = data.frame(
  pred_survival = c(
    1 - pweibull(seq(0, tmax), ova_wei$coefficients[1] %>% exp, ova_wei$coefficients[2] %>% exp),
    1 - pllogis(seq(0, tmax), ova_llo$coefficients[1] %>% exp, ova_llo$coefficients[2] %>% exp)
  ),
  pred_time = rep(
    seq(0, tmax),
    2
  ),
  labs = factor(c(
    rep("Weibull", tmax + 1),
    rep("Log-logistic", tmax + 1)
  ))
)

# test the plot
ggsurvplot(ova_surv, xlim = c(0, tmax), legend.title = "")$plot +
  geom_line(aes(x = pred_time, y = pred_survival, colour = labs), data = df, show.legend = TRUE) +
  scale_colour_discrete(breaks = unique(df$labs)) +
  scale_fill_discrete(guide = "none") +
  theme(axis.text.x = element_blank(), axis.ticks.length.x = unit(0, "mm"))

# do and save the gif
animate(
  ggsurvplot(ova_surv, xlim = c(0, tmax), legend.title = "")$plot +
    geom_line(aes(x = pred_time, y = pred_survival, colour = labs), data = df, show.legend = TRUE) +
    scale_colour_discrete(breaks = unique(df$labs)) +
    scale_fill_discrete(guide = "none") +
    theme(axis.text.x = element_blank(), axis.ticks.length.x = unit(0, "mm")) +
    transition_reveal(pred_time),
  renderer = gifski_renderer(),
  res = 300, width = 1600, height = 1000
)
anim_save(paste0(save_path,"curves.gif"))

### Second plot: parametric uncertainty

# bootstrap weibull parameters and calculate predictions
nBoot = 10
df2 = data.frame(
  pred_survival = apply(
    normboot.flexsurvreg(ova_wei, nBoot), 1, function(x) {1 - pweibull(seq(0, tmax), x[1], x[2])}
  ) %>% c(),
  pred_time = rep(seq(0, tmax), nBoot),
  label = as.factor(sort(rep(1:nBoot, tmax + 1)))
)

# test the plot
ggsurvplot(ova_surv, xlim = c(0, tmax), legend.title = "")$plot +
  geom_line(aes(x = pred_time, y = pred_survival, colour = label), data = df2) +
  scale_fill_discrete(guide = "none") +
  scale_colour_discrete(guide = "none") +
  theme(axis.text.x = element_blank(), axis.ticks.length.x = unit(0, "mm"))

# plot and save the gif
animate(
  ggsurvplot(ova_surv, xlim = c(0, tmax), legend.title = "")$plot +
    geom_line(aes(x = pred_time, y = pred_survival, colour = label), data = df2) +
    scale_fill_discrete(guide = "none") +
    scale_colour_discrete(guide = "none") +
    theme(axis.text.x = element_blank(), axis.ticks.length.x = unit(0, "mm")) +
    transition_reveal(pred_time),
  renderer = gifski_renderer(),
  res = 300, width = 1600, height = 1000
)
anim_save(paste0(save_path,"uncertainty.gif"))

### Third plot: compare against somehow less effective dummy real world data

# ancillary dataset - reduce efficacy and add noise to times 
df3 =
  ppp$data %>%
  mutate(
    surv = surv ^ 2.5,
    strata = "RWE",
    time2 = c(0, sort(abs(time[-1] + rnorm(9, 0, 50))))
  )

# test the plot
ggsurvplot(ova_surv, legend.title = "")$plot +
  geom_step(aes(x = time2, y = surv, colour = strata), data = df3, lwd = 1) +
  scale_colour_discrete(labels = c("RCT", "RWD")) +
  scale_fill_discrete(guide = "none")

# make and save the gif
animate(
  ggsurvplot(ova_surv, legend.title = "")$plot +
    geom_step(aes(x = time2, y = surv, colour = strata), data = df3, lwd = 1) +
    scale_colour_discrete(labels = c("RCT", "RWD")) +
    scale_fill_discrete(guide = "none") +
    transition_reveal(time2),
  renderer = gifski_renderer(),
  res = 300, width = 1600, height = 1000
)
anim_save(paste0(save_path,"RWD.gif"))
