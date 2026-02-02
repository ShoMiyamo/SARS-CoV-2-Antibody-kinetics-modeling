library(readr)
library(dplyr)
library(rstan)
library(ggplot2)
library(cowplot)

supdata<-read_csv('supplemental_data.csv')
Hmodel3 <- stan_model('hierarchical M3.stan')

d3_stan<- supdata

#data curation
d3_stan$infection_status0=case_when(d3_stan$infectionstatus3=="uninfected"~1,d3_stan$infectionstatus3=="Omicron"~3,d3_stan$infectionstatus3=="Pre-Omicron"~2)
d3_stan<-d3_stan %>% mutate(
  vax_infection=as.numeric(case_when(
    vax==2&infection_status0==1~1,vax==3&infection_status0==1~2,vax==4&infection_status0==1~3,vax==5&infection_status0==1~3,
    vax==2&infection_status0==2~4,vax==3&infection_status0==2~5,vax==4&infection_status0==2~6,vax==5&infection_status0==2~6,
    vax==2&infection_status0==3~7,vax==3&infection_status0==3~8,vax==4&infection_status0==3~9,vax==5&infection_status0==3~9,T~10
  )),
  vax_mixed=factor(               
    case_when(
      vax == 2     ~ "Twice",
      vax == 3     ~ "3 times",
      vax >= 4     ~ "4-5 times"
    ),
    levels = c("Twice", "3 times", "4-5 times")
  ))

infection_status<-unique(d3_stan[ , c('vax_infection','infection_status0')])$infection_status0

data_infectionst <-list(N=nrow(d3_stan),Y=d3_stan$anti_s,X=d3_stan$general_intervals,K=9,S=3,VAX=as.numeric(d3_stan$vax_infection),AGE=infection_status,AGE0=d3_stan$infection_status0)

#fitting
fit_inf <- sampling(Hmodel3,
                    data=data_infectionst,
                    iter = 9000,
                    warmup = 1000,
                    thin = 5,
                    chains = 4, core=4,
                    seed=1235,
                    control = list(adapt_delta = 0.99, max_treedepth=15)
                    )

#visualise
library(dplyr)
library(ggplot2)
library(cowplot)
library(readr)   # write_excel_csv

ms <- rstan::extract(fit_infectionstexp2)

x_days  <- 0:730
q_probs <- c(0.025, 0.25, 0.50, 0.75, 0.975)

d3_stan <- d3_stan %>% mutate(Ylin = anti_s)

# ---- quantile summary ----
ci_df <- function(mat, x, probs = q_probs) {
  mat <- as.matrix(mat)
  qs  <- apply(mat, 2, quantile, probs = probs)
  
  out <- data.frame(X = x, t(qs), check.names = FALSE)
  names(out) <- c("X", "q025", "q25", "q50", "q75", "q975")
  out
}

# ---- exp-decay curve draws: base = c + a*exp(-b*x) ----
draw_exp_curve <- function(c, a, b, s_y, x, seed = 1234, probs = q_probs) {
  c   <- drop(c)
  a   <- drop(a)
  b   <- drop(b)
  s_y <- drop(s_y)
  
  n  <- length(a)
  nx <- length(x)
  
  if (length(c) == 1) c <- rep(c, n)
  stopifnot(length(c) == n, length(b) == n)
  
  sd_row <- if (length(s_y) == 1) rep(s_y, n) else rep(s_y, length.out = n)
  
  eta  <- exp(-outer(b, x, `*`))      # n x nx
  base <- sweep(eta, 1, a, `*`)       # a*exp(-b*x)
  base <- sweep(base, 1, c, `+`)      # + c
  
  set.seed(seed)
  eps <- rnorm(n * nx, mean = 0, sd = rep(sd_row, times = nx))
  y   <- base + matrix(eps, nrow = n, ncol = nx)
  
  list(
    base_draws = as.data.frame(base),
    pred_draws = as.data.frame(y),
    base_ci    = ci_df(base, x, probs),
    pred_ci    = ci_df(y, x, probs)
  )
}

# ---- 9 groups ----
vax_levels <- c("Twice", "3 times", "4–5 times")
inf_levels <- c("Uninfected", "Pre-Omicron", "Omicron")

group_map <- expand.grid(
  infection_ix = 1:3,
  vax_ix       = 1:3
) %>%
  as_tibble() %>%
  mutate(
    group      = (infection_ix - 1) * 3 + vax_ix,  
    sd_col     = infection_ix,                     
    vaccination = vax_levels[vax_ix],
    infection   = inf_levels[infection_ix]
  ) %>%
  arrange(group)

grp_res <- setNames(
  lapply(seq_len(nrow(group_map)), function(i) {
    g <- group_map$group[i]
    w <- group_map$sd_col[i]
    draw_exp_curve(
      c   = ms$c[, g],
      a   = ms$a[, g],
      b   = ms$b[, g],
      s_y = ms$s_y[, w],
      x   = x_days,
      seed = 1234
    )
  }),
  paste0("g", group_map$group)
)

base_ci <- lapply(grp_res, `[[`, "base_ci")   
pred_ci <- lapply(grp_res, `[[`, "pred_ci")   

# ---- base = a1*exp(-b1*x) ----
age_res <- setNames(
  lapply(1:3, function(z) {
    draw_exp_curve(
      c   = 0,
      a   = ms$a1[, z],
      b   = ms$b1[, z],
      s_y = ms$s_y[, z],
      x   = x_days,
      seed = 1234
    )
  }),
  paste0("age", 1:3)
)

d_infection1.2 <- age_res$age1$pred_ci
d_infection2.2 <- age_res$age2$pred_ci
d_infection3.2 <- age_res$age3$pred_ci

# ---- total: base = c0 + a0*exp(-b0*x) ----
total_res <- draw_exp_curve(
  c   = ms$c0,
  a   = ms$a0,
  b   = ms$b0,
  s_y = ms$s_y,
  x   = x_days,
  seed = 1234
)

d_estTotal <- total_res$base_ci
d_est      <- total_res$pred_ci

colp<-c("#cccccc50","#53B647","#7744a4")

plot_infection_curves <- function(ci_uninf, ci_pre, ci_omi, g_ids,
                                  cols = c("#99999950", "#53B647", "#7744a4"),
                                  title = "",
                                  points_df = d3_stan,
                                  ylim = c(1, 6)) {
  
  pts <- points_df %>% filter(vax_infection %in% g_ids)
  
  ggplot() +
    theme_bw(base_size = 12) +
    geom_point(
      data = pts,
      aes(x = general_intervals, y = Ylin, colour = factor(infection_status0), fill = factor(infection_status0)),
      shape = 1, alpha = .6, size = .5
    ) +
    geom_ribbon(data = ci_uninf, aes(x = X, ymin = q025, ymax = q975), fill = cols[1], alpha = 3/6) +
    geom_line(  data = ci_uninf, aes(x = X, y = q50),  colour = cols[1], size = 1, alpha = .8) +
    geom_ribbon(data = ci_pre,   aes(x = X, ymin = q025, ymax = q975), fill = cols[2], alpha = 3/6) +
    geom_line(  data = ci_pre,   aes(x = X, y = q50),  colour = cols[2], size = 1, alpha = .8) +
    geom_ribbon(data = ci_omi,   aes(x = X, ymin = q025, ymax = q975), fill = cols[3], alpha = 3/6) +
    geom_line(  data = ci_omi,   aes(x = X, y = q50),  colour = cols[3], size = 1, alpha = .8) +
    coord_cartesian(ylim = ylim) +
    scale_colour_manual(values = colp) +
    
    labs(
      y = expression(paste("Anti S titer (log"["10"], " BAU/ml)")),
      x = "Days post exposure",
      colour = "Infection Status",
      title = title
    ) +
    theme(
      axis.ticks = element_line(colour = "black"),
      axis.ticks.length = grid::unit(-3, "mm"),
      axis.text = element_text(colour = "black", size = 12),
      axis.title.x = element_text(colour = "black"),
      axis.title.y = element_text(colour = "black"),
      plot.title = element_text(colour = "black", size = 12),
      legend.position = "none",
      plot.tag = element_text(face = "bold")
    )
}


Vp1 <- plot_infection_curves(base_ci$g1, base_ci$g4, base_ci$g7, g_ids = c(1, 4, 7), title = "Twice vaccinated")
Vp2 <- plot_infection_curves(base_ci$g2, base_ci$g5, base_ci$g8, g_ids = c(2, 5, 8), title = "3 times vaccinated")
Vp3 <- plot_infection_curves(base_ci$g3, base_ci$g6, base_ci$g9, g_ids = c(3, 6, 9), title = "4 - 5 times")

VpiO <- plot_grid(Vp1, Vp2, Vp3, ncol = 2)
VpiO


#make csv for infection risk estimation
est_all_v2 <- bind_rows(
  lapply(seq_len(nrow(group_map)), function(i) {
    key <- paste0("g", group_map$group[i])
    base_ci[[key]] %>%
      transmute(
        vaccination = group_map$vaccination[i],
        infection   = group_map$infection[i],
        days        = X,
        q025, q25, q50, q75, q975
      )
  })
) %>%
  mutate(
    vaccination = factor(vaccination, levels = vax_levels),
    infection   = factor(infection, levels = inf_levels)
  ) %>%
  relocate(vaccination, infection, days)

write_excel_csv(est_all_v2, "anti_s_est_all_by_group.csv")

