
# ---- packages ----
library(readr)
library(dplyr)
library(rstan)
library(ggplot2)
library(cowplot)

supdata<-read_csv('supplemental_data.csv') #Source Data is in the manuscript.
model3 <- stan_model("stan/M3.stan") #M1.stan and M2.stan are interchangeable

d3_stan<- supdata %>%  filter(infected=="uninfected")
data_vax <-list(N=nrow(d3_stan),Y=d3_stan$anti_s,X=d3_stan$intervals,K=4,VAX=d3_stan$vax)

fit_vax <- sampling(model3,
                        data=data_vax,
                        iter = 9000,
                        warmup = 1000,
                        thin = 5,
                        chains = 4, core=4,
                        seed=1235,
                        control = list(adapt_delta = 0.99, max_treedepth=15))

#visualize
# ---- posterior extract ----
post <- rstan::extract(fit_vax)

x_days <- 0:730
q_probs <- c(0.025, 0.25, 0.50, 0.75, 0.975)

ci_df <- function(mat, x, probs = q_probs) {
  mat <- as.matrix(mat)
  qs  <- apply(mat, 2, quantile, probs = probs)
  data.frame(X = x, t(qs), check.names = FALSE)
}

draw_exp_curve <- function(c, a, b, s_y, x, seed = 1234) {
  c   <- drop(c)
  a   <- drop(a)
  b   <- drop(b)
  s_y <- drop(s_y)
  
  n   <- length(c)
  nx  <- length(x)
  
  stopifnot(length(a) == n, length(b) == n)
  
  eta  <- exp(-outer(b, x, `*`))
  base <- sweep(eta, 1, a, `*`)
  base <- sweep(base, 1, c, `+`)
  
  set.seed(seed)
  
  sd_vec <- rep(s_y, times = nx)
  sd_vec <- rep(sd_vec, length.out = n * nx)
  
  eps <- rnorm(n * nx, mean = 0, sd = sd_vec)
  y   <- base + matrix(eps, nrow = n, ncol = nx)
  
  list(
    base_draws = as.data.frame(base),
    pred_draws = as.data.frame(y),
    base_ci    = ci_df(base, x),
    pred_ci    = ci_df(y, x)
  )
}


# ---- generate per-vax and total ----
vax_ids <- 1:4

vax_res <- setNames(
  lapply(vax_ids, function(j) {
    draw_exp_curve(
      c   = post$c[, j],
      a   = post$a[, j],
      b   = post$b[, j],
      s_y = post$s_y,
      x   = x_days,
      seed = 1234
    )
  }),
  paste0("vax", vax_ids)
)

total_res <- draw_exp_curve(
  c   = post$c0,
  a   = post$a0,
  b   = post$b0,
  s_y = post$s_y,
  x   = x_days,
  seed = 1234
)

ci_base    <- lapply(vax_res, `[[`, "base_ci")     
ci_pred    <- lapply(vax_res, `[[`, "pred_ci")     
pred_draws <- lapply(vax_res, `[[`, "pred_draws")  

d_estTotal <- total_res$base_ci   
d_est      <- total_res$pred_ci   
d_est.1    <- total_res$pred_draws

# ---- plotting helpers ----
colp <- c(
  "#cccccc50", "#9966b450", "#39568CFF", "#53B647",
  "#d3c848",   "#440154FF", "#39568CFF", "#55C667FF", "#C9B848"
)

plot_ci <- function(ci, vax_id = NULL, title = "", fill = "#7744a4",
                    points_df = d3_stan, overall_ci = d_est) {
  pts <- if (is.null(vax_id)) points_df else filter(points_df, vax == vax_id)
  
  ggplot() +
    theme_bw(base_size = 12) +
    geom_point(
      data = pts,
      aes(x = intervals, y = Ylin, colour = as.factor(vax)),
      shape = 1, alpha = .6, size = .5
    ) +
    geom_ribbon(data = ci, aes(x = X, ymin = `2.5%`, ymax = `97.5%`), fill = fill, alpha = 1/6) +
    geom_ribbon(data = ci, aes(x = X, ymin = `25%`,  ymax = `75%`),  fill = fill, alpha = 2/6) +
    geom_line(  data = ci, aes(x = X, y = `50%`), colour = fill, size = 1, alpha = .8) +
    geom_line(  data = overall_ci, aes(x = X, y = `50%`), linetype = 3, colour = "gray20", size = .8) +
    coord_cartesian(ylim = c(1, 6)) +
    scale_color_manual(values = colp) +
    labs(
      y = expression(paste("Anti S titer (log"["10"], " BAU/ml)")),
      x = "Days post vaccination",
      color = "Number of vaccination",
      title = title
    ) +
    theme(
      axis.ticks = element_line(colour = "black"),
      axis.ticks.length = grid::unit(-3, "mm"),
      axis.title.x = element_text(colour = "black"),
      axis.title.y = element_text(colour = "black"),
      axis.text = element_text(colour = "black", size = 12),
      plot.title = element_text(colour = "black", size = 12),
      legend.position = "none",
      plot.tag = element_text(face = "bold", size = 12)
    )
}

# ---- colors/titles ----
fill_cols  <- c("#7744a4", "#39568CFF", "#53B647", "#b3a828")
titles     <- c("Twice vaccinated", "3 times vaccinated", "4 times vaccinated", "5 times vaccinated")
line_sizes <- c(1, 0.9, 1, 0.8)

# ---- CI plots (mean curve) ----
Vc_plots <- Map(
  function(ci, id, col, ttl) plot_ci(ci, vax_id = id, title = ttl, fill = col),
  ci_base, vax_ids, fill_cols, titles
)
Vc1 <- Vc_plots[[1]]; Vc2 <- Vc_plots[[2]]; Vc3 <- Vc_plots[[3]]; Vc4 <- Vc_plots[[4]]
Vc  <- plot_grid(plotlist = Vc_plots, ncol = 2)
Vc

# ---- predictive plots ----
Vp_plots <- Map(
  function(ci, id, col, ttl) plot_ci(ci, vax_id = id, title = ttl, fill = col),
  ci_pred, vax_ids, fill_cols, titles
)
Vp1 <- Vp_plots[[1]]; Vp2 <- Vp_plots[[2]]; Vp3 <- Vp_plots[[3]]; Vp4 <- Vp_plots[[4]]

# ---- overlay mean plot  ----
vpt <- ggplot() +
  theme_bw(base_size = 12) +
  coord_cartesian(ylim = c(1, 6)) +
  labs(
    y = expression(paste("Anti S titer (log"["10"], " BAU/ml)")),
    x = "Days post vaccination",
    title = "Mean",
    tag = ""
  ) +
  theme(
    axis.ticks = element_line(colour = "black"),
    axis.ticks.length = grid::unit(-3, "mm"),
    axis.text = element_text(colour = "black", size = 12),
    axis.title.y = element_text(colour = "black", face = "bold"),
    plot.title = element_text(colour = "black", size = 12),
    legend.position = "none",
    plot.tag = element_text(face = "bold", size = 12)
  )

# vax1 ribbon + line first
vpt <- vpt +
  geom_ribbon(data = ci_base[[1]], aes(x = X, ymin = `2.5%`, ymax = `97.5%`),
              fill = fill_cols[1], alpha = 2/6) +
  geom_line(data = ci_base[[1]], aes(x = X, y = `50%`),
            colour = fill_cols[1], size = line_sizes[1], alpha = 0.6)

# ribbons for vax2-4
for (j in 2:4) {
  vpt <- vpt + geom_ribbon(
    data = ci_base[[j]],
    aes(x = X, ymin = `2.5%`, ymax = `97.5%`),
    fill = fill_cols[j], alpha = 2/6
  )
}

for (j in c(3, 2, 4)) {
  vpt <- vpt + geom_line(
    data = ci_base[[j]],
    aes(x = X, y = `50%`),
    colour = fill_cols[j], size = line_sizes[j], alpha = 0.6
  )
}

Vpp <- plot_grid(plotlist = c(Vp_plots, list(vpt)), ncol = 3, align = "hv")
Vpp
