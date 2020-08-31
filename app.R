


# Packages ----------------------------------------------------------------


library(shiny)
library(shinythemes)
library(tidyverse)
library(stats4)
library(scales)
library(plotly)
library(shinycssloaders)
library(googledrive)
library(googlesheets4)
library(DT)
library(tidybayes)
library(ggtext)
library(cowplot)


# Load googlesheet wihtout authorisation --------------------------------------------------------


gs4_deauth()
ss <- "1VcOl0hN_DRiq9BZyn3rZQ1oa9Bg6en1HSvqGeGPx7HU"

sheet <-
  read_sheet(ss) %>% mutate(Link = paste0("<a href='", Link, "' target='_blank'>", Link, "</a>"))


#N
n_rep <- 1e5

####Functions -----

# https://stats.stackexchange.com/questions/12232/calculating-the-parameters-of-a-beta-distribution-using-the-mean-and-variance/12239#12239
normalbeta <- function(mu, var) {
  alpha0 <- ((1 - mu) / var ^ 2 - 1 / mu) * mu ^ 2
  beta0 <- alpha0 * (1 / mu - 1)
  list(alpha0 = alpha0, beta0 = beta0)
}


# exact solution function, based on "Introduction to empirical Bayes", Robinson
h <- function(alpha_a, beta_a,
              alpha_b, beta_b) {
  j <- seq.int(0, round(alpha_b) - 1)
  log_vals <-
    (
      lbeta(alpha_a + j, beta_a + beta_b) - log(beta_b + j) -
        lbeta(1 + j, beta_b) - lbeta(alpha_a, beta_a)
    )
  1 - sum(exp(log_vals))
}
#Absolute difference closed solution
ebd_fun <- function(alpha0, beta0, alphaA, betaA, alphaB, betaB) {
  tibble(
    absolute = rbeta(n_rep, shape1 = alphaA + alpha0, shape2 = betaA + beta0) - rbeta(n_rep, shape1 = alphaB + alpha0, shape2 = betaB + beta0)
  )
}

#Relative risk  closed solution
ebrr_fun <- function(alpha0, beta0, alphaA, betaA, alphaB, betaB) {
  tibble(
    rr = rbeta(n_rep, shape1 = alphaA + alpha0, shape2 = betaA + beta0) / rbeta(n_rep, shape1 = alphaB + alpha0, shape2 = betaB + beta0)
  )
}
#Odds ratio   closed solution

ebor_fun <- function(alpha0, beta0, alphaA, betaA, alphaB, betaB) {
  a <- rbeta(n_rep, shape1 = alphaA + alpha0, shape2 = betaA + beta0)
  b <-
    rbeta(n_rep, shape1 = alphaB + alpha0, shape2 = betaB + beta0)
  
  tibble(or = (a / (1 - a)) / (b / (1 - b)))
}
###### UI ########

ui <-
  navbarPage(
    title = "Empirical Bayes (beta version)",
    theme = shinytheme("lumen"),
    tabPanel('Introduction',
             includeMarkdown('Introduction.Rmd')),
    tabPanel("Single study",
             fluidPage(
               sidebarLayout(sidebarPanel(
                 tabsetPanel(
                   id = "sidebar",
                   tabPanel(
                     "Main controls",
                     br(),
                     selectInput(
                       "dist",
                       "Prior",
                       choices = c(
                         "Mean and SD" = "normal",
                         "Data derived beta" = "empirical.beta",
                         "Custom beta" = "beta",
                         "Uniform" = "neutral"
                       ),
                       selected = "empirical.beta"
                     ),
                     conditionalPanel(
                       condition = "input.dist == 'normal'",
                       numericInput(
                         "mu",
                         "Mean",
                         min = 0,
                         max = 1,
                         value = 0.25
                       ),
                       numericInput(
                         "sigma",
                         "Standard deviation",
                         min = 0,
                         value = 0.1
                       )
                     ),
                     conditionalPanel(
                       condition = "input.dist== 'beta'",
                       h4("Beta parameters"),
                       numericInput("shape1", "Alpha (events)", value = 1),
                       numericInput("shape2", "Beta (non-events)", value = 1)
                     ),
                     br(),
                     h4("Events and group sizes"),
                     numericInput(
                       "eventA",
                       "Number of events in intevention/exposed group",
                       step = 1,
                       value = 30
                     ),
                     numericInput(
                       "sizeA",
                       "Size of intervention/exposed group",
                       step = 1,
                       value = 100
                     ),
                     numericInput(
                       "eventB",
                       "Number of events in control group",
                       step = 1,
                       value = 40
                     ),
                     numericInput(
                       "sizeB",
                       "Size of control group",
                       step = 1,
                       value = 100
                     ),
                     radioButtons(
                       inputId = "direction",
                       label = "Difference of interest",
                       choices = c(
                         "Lower eventrate in intervention" = "lower",
                         "Higher eventrate in intervention" = "higher"
                       ),
                       selected = "lower"
                     ),
                     radioButtons(
                       inputId = "difftype",
                       label = "Absolute, RR or OR?",
                       choices = c(
                         "Absolute" = "absolute",
                         "RR" = "rr",
                         "OR" = "or"
                       ),
                       selected = "absolute",
                       inline = TRUE
                     ),
                     conditionalPanel(
                       "input.difftype=='absolute'",
                       numericInput(
                         inputId = "mcid_a",
                         label = "Minimum clinically important difference, absolute",
                         step = .01,
                         value = 0
                       )
                     ),
                     conditionalPanel(
                       "input.difftype=='rr'",
                       numericInput(
                         inputId = "mcid_rr",
                         label = "Minimum clinically important difference, RR",
                         step = .01,
                         value = 1
                       )
                     ),
                     conditionalPanel(
                       "input.difftype=='or'",
                       numericInput(
                         inputId = "mcid_or",
                         label = "Minimum clinically important difference, OR",
                         step = .01,
                         value = 1
                       )
                     ),
                     numericInput(
                       "cred",
                       "Credibility interval",
                       step = 0.01,
                       value = 0.95,
                       min = 0.01,
                       max = 0.99
                     ),
                     actionButton("goB", "Go!")
                   ),
                   tabPanel(
                     "Options",
                     radioButtons(
                       inputId = "pointest",
                       label = "Point estimate",
                       choices = c(
                         "Mean" = "mean",
                         "Median" = "median",
                         "Mode" = "mode"
                       ),
                       selected = "mean"
                     ),
                     radioButtons(
                       inputId = "interval",
                       label = "Credible interval definition?",
                       choices = c(
                         "Highest posterior density interval" = "hpdi",
                         "Equal-tailed interval" = "pi"
                       ),
                       selected = "hpdi"
                     )
                   )
                 )
               ),
               mainPanel(
                 tabsetPanel(
                   id = "viewpane",
                   tabPanel(
                     "Probability density plots",
                     value = "probplots",
                     plotOutput("posteriorplot") %>% withSpinner(),
                     plotOutput("diffplot") %>% withSpinner(),
                     textOutput("HPDI"),
                     h6(htmlOutput("exact"))
                   ),
                   tabPanel("Heatmap",
                            plotlyOutput("heat") %>% withSpinner()),
                   tabPanel("Select data from table", value = "datatable", DTOutput("dt"))
                 )
               ))
             )),
    tabPanel(
      "Comparison of multiple studies",
      fluidPage(
        title = "Empirical Bayes",
        titlePanel("Empirical Bayes"),
        sidebarLayout(sidebarPanel(
          tabsetPanel(
            id = "sidebar",
            tabPanel(
              "Main controls",
              br(),
              radioButtons(
                inputId = "difftypem",
                label = "Absolute, RR or OR?",
                choices = c(
                  "Absolute" = "absolute",
                  "RR" = "rr",
                  "OR" = "or"
                ),
                selected = "absolute",
                inline = TRUE
              ),
              actionButton("goM", "Go!")
            ),
            tabPanel("Options",
                     actionButton("goR", "Refresh"))
          )
        ), mainPanel(
          tabsetPanel(
            id = "viewpanem",
            tabPanel("Data table", value = "datatablem", DTOutput("dtm") %>% withSpinner()),
            tabPanel("Probability density plots",
                     value = "multiplot", withSpinner(
                       plotOutput("multiplot", click = "point_click", brush = "brush")
                     )),
            tabPanel(
              "Empirical prior and likelihoods",
              plotOutput("multipriorplot") %>% withSpinner()
            )
          )
        ))
      )
    ),
    tabPanel('About',
             fluidPage(title = 'About',includeMarkdown ('Background.Rmd'))
             )
  )





###### SERVER #####


server <- function(input, output, session) {
  # sheet <- eventReactive(input$goB, ignoreNULL = FALSE, {
  #   read_sheet(ss)
  # })
  observeEvent(input$dt_rows_selected, {
    updateNumericInput(session, inputId = "eventA", value = as.numeric(sheet[input$dt_rows_selected, "Intervention events"]))
    updateNumericInput(session, inputId = "sizeA", value = as.numeric(sheet[input$dt_rows_selected, "Intervention n"]))
    updateNumericInput(session, inputId = "eventB", value = as.numeric(sheet[input$dt_rows_selected, "Control events"]))
    updateNumericInput(session, inputId = "sizeB", value = as.numeric(sheet[input$dt_rows_selected, "Control n"]))
    updateActionButton(session, inputId = "goB")
  })
  
  observeEvent(input$goB, {
    if (input$viewpane == "datatable") {
      updateTabsetPanel(session = session,
                        inputId = "viewpane",
                        selected = "probplots")
    }
  })
  
  
  # Calculate prior parameters------
  prior_param <- eventReactive(input$goB, ignoreNULL = FALSE, {
    if (input$dist == "normal") {
      return(normalbeta(input$mu, input$sigma))
    }
    else {
      if (input$dist == "empirical.beta") {
        return(tibble(
          alpha0 = (input$eventA + input$eventB),
          beta0 = (input$sizeA + input$sizeB) - (input$eventA + input$eventB)
        ))
      }
      else {
        if (input$dist == "beta") {
          return(tibble(alpha0 = input$shape1,
                        beta0 = input$shape2))
        } else {
          if (input$dist == "neutral") {
            return(tibble(alpha0 = 1, beta0 = 1))
          }
        }
      }
    }
  })
  
  
  
  # Sampling of prior, posterior and difference distributions----
  
  
  df <- eventReactive(input$goB, ignoreNULL = FALSE, {
    prior <- rbeta(1e5,
                   shape1 = prior_param()$alpha0,
                   shape2 = prior_param()$beta0)
    alphaA <- input$eventA + prior_param()$alpha0
    betaA <- prior_param()$beta0 + input$sizeA - input$eventA
    alphaB <- input$eventB + prior_param()$alpha0
    betaB <- prior_param()$beta0 + input$sizeB - input$eventB
    likelihood_A <- rbeta(1e5,
                          shape1 = input$eventA,
                          shape2 = input$sizeA - input$eventA)
    likelihood_B <- rbeta(1e5,
                          shape1 = input$eventB,
                          shape2 = input$sizeB - input$eventB)
    prior <-
      rbeta(1e5,
            shape1 = prior_param()$alpha0,
            shape2 = prior_param()$beta0)
    post_A <- rbeta(1e5, shape1 = alphaA, shape2 = betaA)
    post_B <- rbeta(1e5, alphaB, betaB)
    
    tibble(
      prior,
      post_A,
      likelihood_A,
      likelihood_B,
      post_B = post_B,
      d = post_A - post_B,
      or = (post_A / (1 - post_A)) / (post_B / (1 - post_B)),
      rr = post_A / post_B
    )
  })
  
  
  observeEvent(input$goB, ignoreNULL = FALSE, {
    
  })
  
  
  
  
  diff <- reactive({
    if (input$difftype == "absolute") {
      return(df()$d)
    } else {
      if (input$difftype == "rr") {
        return(df()$rr)
      } else {
        return(df()$or)
      }
    }
  })
  
  
  pest <- reactive({
    if (input$pointest == "mean") {
      mean(diff())
    } else {
      if (input$pointest == "median") {
        median(diff())
      } else {
        Mode(diff())
      }
    }
  })
  
  # Plot of prior and posterior distributions ------
  post_plot <-
    reactive({
      df() %>%
        dplyr::select(prior, likelihood_A, likelihood_B) %>%
        pivot_longer(cols = c(1:3),
                     names_to = "distribution",
                     values_to = "density") %>%
        ggplot(aes(x = density, colour = distribution)) +
        geom_density(adjust =
                       2, size = 2) +
        labs(title = "Prior and likelihood distributions") +
        theme_classic() +
        ylab("") +
        xlab("") +
        theme(
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.y = element_blank(),
          legend.title = element_blank(),
          legend.position = c(0.8, 0.8),
          legend.text = element_text(size = 15)
        ) +
        scale_x_continuous() +
        scale_colour_discrete(name =
                                "Probability distributions",
                              labels = c("Intervention", "Control", "Prior"))
    })
  
  
  # exact solution for probability of difference between groups, ----
  exact_solution <- reactive({
    alphaA <- input$eventA + prior_param()$alpha0
    betaA <- prior_param()$beta0 + input$sizeA - input$eventA
    alphaB <- input$eventB + prior_param()$alpha0
    betaB <- prior_param()$beta0 + input$sizeB - input$eventB
    if (input$direction == "lower") {
      (1 - h(alphaA, betaA, alphaB, betaB))
    } else if (input$direction == "higher") {
      h(alphaA, betaA, alphaB, betaB)
    }
  })
  
  # calculate probability of difference exceeding mcid ----
  over_mcid <- reactive({
    df <- df()
    if (input$direction == "lower") {
      if (input$difftype == "absolute") {
        (sum(df$d < input$mcid_a) / 1e5)
      }
      else {
        if (input$difftype == "rr") {
          (sum(df$rr < input$mcid_rr) / 1e5)
        }
        else {
          if (input$difftype == "or") {
            (sum(df$or < input$mcid_or) / 1e5)
          }
        }
      }
    } else {
      if (input$direction == "higher") {
        if (input$difftype == "absolute") {
          (sum(df$d > input$mcid_a) / 1e5)
        }
        else {
          if (input$difftype == "rr") {
            (sum(df$rr > input$mcid_rr) / 1e5)
          }
          else {
            if (input$difftype == "or") {
              (sum(df$or > input$mcid_or) / 1e5)
            }
          }
        }
      }
    }
  })
  
  
  # Probability density plot for difference in eventrates between groups
  plot_diff <-
    reactive({
      req(df())
      df <- df()
      isolate({
        alphaA <- input$eventA + prior_param()$alpha0
        betaA <- prior_param()$beta0 + input$sizeA - input$eventA
        alphaB <- input$eventB + prior_param()$alpha0
        betaB <- prior_param()$beta0 + input$sizeB - input$eventB
      })
      
      
      if (input$difftype == "absolute") {
        ci <- if (input$interval == "hpdi") {
          hdi(df()$d, input$cred)
        } else {
          qi(df()$d, .width = input$cred)
        }
        pest <- pest()
        
        dens <- density(df$d, adjust = 2)
        df <- tibble(x = dens$x, y = dens$y)
        df <-
          df %>% mutate(shade = if (input$direction == "lower") {
            case_when((x <
                         input$mcid_a) ~ "on",
                      TRUE ~ NA_character_)
          } else {
            case_when((x >
                         input$mcid_a) ~ "on",
                      TRUE ~ NA_character_)
          })
        return(
          ggplot(aes(x, y), data = df) +
            geom_line(
              size = 2.5,
              colour = "black",
              alpha = 0.8
            ) +
            geom_area(
              data = filter(df, shade == "on"),
              fill = "brown2",
              alpha = 0.2
            ) +
            geom_vline(
              xintercept = isolate({
                input$mcid_a
              }),
              color = "red",
              linetype = "dashed",
              size = 2
            ) +
            geom_point(aes(x = pest, y = 0), size = 6) +
            geom_segment(aes(
              x = ci[1],
              y = 0,
              xend = ci[2],
              yend = 0
            ), size = 3) +
            theme_classic() +
            xlab(NULL) +
            ylab(NULL) +
            scale_x_continuous(labels = scales::percent) +
            theme(
              axis.text.y = element_blank(),
              axis.ticks.y = element_blank(),
              axis.line.y = element_blank()
            )
        )
      }
      else if (input$difftype == "rr") {
        ci <- if (input$interval == "hpdi") {
          hdi(df()$rr, input$cred)
        } else {
          qi(df()$rr, input$cred)
        }
        pest <- pest()
        
        dens <- density(df$rr, adjust = 2)
        df <- tibble(x = dens$x, y = dens$y)
        df <-
          df %>% mutate(shade = if (input$direction == "lower") {
            case_when((x <
                         input$mcid_rr) ~ "on",
                      TRUE ~ NA_character_)
          } else {
            case_when((x >
                         input$mcid_rr) ~ "on",
                      TRUE ~ NA_character_)
          })
        
        return(
          ggplot(aes(x, y), data = df) +
            geom_line(
              size = 2.5,
              colour = "black",
              alpha = 0.8
            ) +
            geom_area(
              data = filter(df, shade == "on"),
              fill = "brown2",
              alpha = 0.2
            ) +
            geom_vline(
              xintercept = isolate({
                input$mcid_rr
              }),
              color = "red",
              linetype = "dashed",
              size = 2
            ) +
            geom_point(aes(x = pest, y = 0), size = 6) +
            geom_segment(aes(
              x = ci[1],
              y = 0,
              xend = ci[2],
              yend = 0
            ), size = 3) +
            theme_classic() +
            xlab(NULL) +
            ylab(NULL) +
            scale_x_continuous() +
            theme(
              axis.text.y = element_blank(),
              axis.ticks.y = element_blank(),
              axis.line.y = element_blank()
            )
        )
      }
      else if (input$difftype == "or") {
        ci <- if (input$interval == "hpdi") {
          hdi(df()$or, input$cred)
        } else {
          qi(df()$or, .width = input$cred)
        }
        pest <- pest()
        
        dens <- density(df$or, adjust = 2)
        df <- tibble(x = dens$x, y = dens$y)
        df <-
          df %>% mutate(shade = if (input$direction == "lower") {
            case_when((x <
                         input$mcid_or) ~ "on",
                      TRUE ~ NA_character_)
          } else {
            case_when((x >
                         input$mcid_or) ~ "on",
                      TRUE ~ NA_character_)
          })
        
        return(
          ggplot(aes(x, y), data = df) +
            geom_line(
              size = 2.5,
              colour = "black",
              alpha = 0.8
            ) +
            geom_area(
              data = filter(df, shade == "on"),
              fill = "brown2",
              alpha = 0.2
            ) +
            geom_vline(
              xintercept = isolate({
                input$mcid_or
              }),
              color = "red",
              linetype = "dashed",
              size = 2
            ) +
            geom_point(aes(x = pest, y = 0), size = 6) +
            geom_segment(aes(
              x = ci[1],
              y = 0,
              xend = ci[2],
              yend = 0
            ), size = 3) +
            theme_classic() +
            xlab(NULL) +
            ylab(NULL) +
            scale_x_continuous() +
            theme(
              axis.text.y = element_blank(),
              axis.ticks.y = element_blank(),
              axis.line.y = element_blank()
            )
        )
      }
    })
  
  # heatmap of posterior probability of difference in event rate depending on mu and sigma of prior
  heatmap <- eventReactive(input$goB, ignoreNULL = FALSE, {
    a1 <- input$eventA
    b1 <- input$sizeA - input$eventA
    a2 <- input$eventB
    b2 <- input$sizeB - input$eventB
    
    heat <- function(mu, sigma) {
      nb <- normalbeta(mu, sigma)
      alpha0 <- nb[[1]]
      beta0 <- nb[[2]]
      
      alpha_a <- alpha0 + a1
      beta_a <- beta0 + b1
      alpha_b <- alpha0 + a2
      beta_b <- beta0 + b2
      
      if (input$direction == "lower") {
        h(alpha_a, beta_a, alpha_b, beta_b)
      } else if (input$direction == "higher") {
        1 - h(alpha_a, beta_a, alpha_b, beta_b)
      }
    }
    df <-
      tibble(
        Mu = seq(0.1, 0.9, length.out = 100),
        Sigma = seq(0.01, 0.05, length.out = 100)
      )
    
    df <- expand.grid(df)
    df %>% mutate(Probability = map2(Mu, Sigma, heat))
  })
  
  
  
  # Multianalysis -------------------------------------------------------------------
  
  diff_m <- eventReactive(input$goM, ignoreNULL = FALSE, {
    req(input$dtm_rows_selected)
    
    sheet %>%
      as_tibble() %>%
      slice(input$dtm_rows_selected) %>%
      janitor::clean_names() %>%
      mutate(
        alpha0 = intervention_events + control_events,
        beta0 = intervention_n + control_n - alpha0,
        alphaA = intervention_events,
        betaA = intervention_n - alphaA,
        alphaB = control_events,
        betaB = control_n - alphaB,
        comb = paste0(
          intervention,
          " vs.",
          control,
          ",\n ",
          population,
          ",\n ",
          endpoint
        )
      ) %>%
      select(-c(field, link)) %>%
      mutate(
        absolute = pmap(list(
          alpha0, beta0, alphaA, betaA, alphaB, betaB
        ), ebd_fun),
        rr = pmap(list(
          alpha0, beta0, alphaA, betaA, alphaB, betaB
        ), ebrr_fun),
        or = pmap(list(
          alpha0, beta0, alphaA, betaA, alphaB, betaB
        ), ebor_fun),
        higher = pmap_dbl(list(
          alpha_a = (alphaA + alpha0),
          beta_a = (betaA + beta0),
          alpha_b = (alphaB + alpha0),
          beta_b = (betaB + beta0)
        ), h),
        lower = pmap_dbl(list(
          alpha_a = (alphaB + alpha0),
          beta_a = (betaB + beta0),
          alpha_b = (alphaA + alpha0),
          beta_b = (betaA + beta0)
        ), h)
      ) %>%
      unnest(input$difftypem) %>%
      select(
        endpoint,
        population,
        comb,
        diff = input$difftypem,
        higher,
        lower,
        alpha0,
        beta0,
        alphaA,
        betaA,
        alphaB,
        betaB
      )
  })
  
  # Multiplot ---------------------------------------------------------------
  
  
  mplot <- eventReactive(input$goM,
                         {
                           req(input$dtm_rows_selected)
                           ebd <- diff_m()
                           
                           clim <-
                             tibble(
                               xmn = switch(
                                 input$difftypem,
                                 absolute = -max(abs(ebd$diff)),
                                 rr = (1 - max(abs((
                                   ebd$diff - 1
                                 )))),
                                 or = (1 - max(abs((
                                   ebd$diff - 1
                                 ))))
                               ),
                               xmx = switch(
                                 input$difftypem,
                                 absolute = max(abs(ebd$diff)),
                                 rr = 1 + max(abs((ebd$diff - 1))),
                                 or = 1 + max(abs((ebd$diff - 1)))
                               ),
                               ymn = 1,
                               ymx = (length(unique(ebd$comb)) + 0.5),
                               xcmn = xmn + (xmx - xmn) * 0.2,
                               xcmx = xmx - (xmx - xmn) * 0.2
                             )
                           cords <- ebd %>%
                             group_by(comb) %>%
                             nest() %>%
                             mutate(hpdi = map(.x = data, .f = ~ hdi(.x$diff, .width = 0.95)),
                                    point = map(data, ~ mean(.x$diff))) %>%
                             unnest() %>%
                             bind_cols(clim) %>%
                             select(-diff) %>%
                             distinct() %>%
                             arrange(higher) %>%
                             mutate(
                               higher = scales::percent(higher),
                               lower = scales::percent(lower),
                               hpdimin = ifelse(
                                 isolate({
                                   input$difftypem
                                 }) == "absolute",
                                 scales::percent(hpdi[, 1], 0.1),
                                 round(hpdi[, 1], 2)
                               ),
                               hpdimax = ifelse(
                                 isolate({
                                   input$difftypem
                                 }) == "absolute",
                                 scales::percent(hpdi[, 2], 0.1),
                                 round(hpdi[, 2], 2)
                               ),
                               point = ifelse(
                                 isolate({
                                   input$difftypem
                                 }) == "absolute",
                                 scales::percent(point, 0.1),
                                 round(point, 2)
                               )
                             ) %>%
                             select(-hpdi)
                           
                           
                           ebd %>%
                             ggplot(aes(
                               x = diff,
                               y = fct_reorder(comb, higher),
                               fill = stat(x > (ifelse(
                                 isolate({
                                   input$difftypem
                                 })
                                 == "absolute", 0, 1
                               )))
                             )) +
                             stat_halfeye() +
                             geom_vline(xintercept = (switch(
                               isolate({
                                 input$difftypem
                               }),
                               absolute = 0,
                               rr = 1,
                               or = 1
                             )), linetype = "dotted") +
                             geom_text(aes(x = xmx,
                                           y = comb,
                                           label = higher),
                                       data = cords,
                                       inherit.aes = FALSE) +
                             geom_text(aes(x = xmn,
                                           y = comb,
                                           label = lower),
                                       data = cords,
                                       inherit.aes = FALSE) +
                             geom_richtext(
                               aes(
                                 x = xcmx,
                                 y = ymx,
                                 label = "Posterior probability of *any* negative effect",
                                 fill = NA,
                                 label.color = NA
                               ),
                               inherit.aes = FALSE,
                               data = cords,
                               nudge_y = 0.2
                             ) +
                             geom_richtext(
                               aes(
                                 x = xcmn,
                                 y = ymx,
                                 label = ("Posterior probability of *any* positive effect"),
                                 fill = NA,
                                 label.color = NA
                               ),
                               inherit.aes = FALSE,
                               data = cords,
                               nudge_y = 0.2
                             ) +
                             labs(
                               x = NULL,
                               y = NULL,
                               fill = NULL,
                               title = "Empirical Bayes analysis of results",
                               subtitle = "Empirical bayesian analysis of relative risk between treatment groups.\nData derived beta prior based on the sum of events in the two groups",
                               caption = "@load_dependent"
                             ) +
                             scale_fill_manual(values = c("skyblue", "pink")) +
                             # coord_cartesian(xlim = c(clim$xmn, clim$xmx), ylim = c(clim$ymn, clim$ymx)) +
                             scale_y_discrete(
                               position = "left",
                               labels = paste0(
                                 cords$comb,
                                 "\n",
                                 
                                 cords$point,
                                 " (",
                                 
                                 cords$hpdimin,
                                 " to ",
                                 
                                 cords$hpdimax,
                                 ")"
                               )
                             ) +
                             theme_tidybayes() +
                             theme(legend.position = "none") +
                             panel_border() +
                             scale_x_continuous(
                               n.breaks = 9,
                               limits = c(clim$xmn, clim$xmx),
                               labels = switch(
                                 input$difftypem,
                                 absolute = scales::percent,
                                 rr = waiver(),
                                 or = waiver()
                               )
                             )
                         })
  
  
  multipriorplot <- eventReactive(input$goM, {
    req(input$dtm_rows_selected)
    sheet %>%
      as_tibble() %>%
      slice(input$dtm_rows_selected) %>%
      janitor::clean_names() %>%
      mutate(
        alpha0 = intervention_events + control_events,
        beta0 = intervention_n + control_n - alpha0,
        alphaA = intervention_events,
        betaA = intervention_n - alphaA,
        alphaB = control_events,
        betaB = control_n - alphaB,
        comb = as.factor(
          paste0(intervention, ",\n ", population, ",\n Endpoint:", endpoint)
        )
      ) %>%
      mutate(
        Prior = pmap(list(n = n_rep, alpha0, beta0), rbeta),
        Intervention = pmap(list(n = n_rep, alphaA, betaA), rbeta),
        Control = pmap(list(n = n_rep, alphaB, betaB), rbeta)
      ) %>%
      unnest(Prior, Intervention, Control) %>%
      select(comb, Prior, Intervention, Control) %>%
      pivot_longer(-comb) %>%
      ggplot(aes(x = value, colour = name)) +
      geom_density(size = 1.5, adjust = 2) +
      facet_wrap( ~ comb, ncol = 1, scales = "free") +
      theme_minimal() +
      xlab(NULL) +
      ylab(NULL) +
      scale_x_continuous() +
      theme(
        legend.title = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank()
      )
  })
  observeEvent(input$goM, {
    if (input$viewpanem == "datatablem") {
      req(input$dtm_rows_selected)
      updateTabsetPanel(session = session,
                        inputId = "viewpanem",
                        selected = "multiplot")
    }
  })
  
  observeEvent(input$goR, {
    session$reload()
  })
  
  # Output --------
  output$posteriorplot <- renderPlot(post_plot())
  output$diffplot <- renderPlot(plot_diff())
  
  output$exact <- renderUI(HTML(paste(
    paste0
    (
      "Exact solution for probability of ",
      input$direction,
      " eventrate in intervention group is ",
      percent(exact_solution(), accuracy = 0.1)
    ),
    paste0(" Mean ",
           if (input$difftype == "absolute") {
             paste0(
               "absolute difference is ",
               percent(pest(), .01),
               " (",
               percent(input$cred),
               " credible interval ",
               percent(ifelse(
                 input$interval == "hpdi",
                 hdi(df()$d, input$cred)[1],
                 qi(df()$d, input$cred)[1]
               ), .01),
               " to ",
               percent(ifelse(
                 input$interval == "hpdi",
                 hdi(df()$d, input$cred)[2],
                 qi(df()$d, input$cred)[2]
               ), .01),
               ")."
             )
           } else {
             if (input$difftype == "rr") {
               paste0(
                 "relative risk is ",
                 round(pest(), 2),
                 " (",
                 percent(input$cred),
                 " credible interval ",
                 round(ifelse(
                   input$interval == "hpdi",
                   hdi(df()$rr, input$cred)[1],
                   qi(df()$rr, input$cred)[1]
                 ), 2),
                 " to ",
                 round(ifelse(
                   input$interval == "hpdi",
                   hdi(df()$rr, input$cred)[2],
                   qi(df()$rr, input$cred)[2]
                 ), 2),
                 ")."
               )
             } else {
               if (input$difftype == "or") {
                 paste0(
                   "odds ratio is ",
                   round(pest(), 2),
                   " (",
                   percent(input$cred),
                   ", credible interval ",
                   round(ifelse(
                     input$interval == "hpdi",
                     hdi(df()$or, input$cred)[1],
                     qi(df()$or, input$cred)[1]
                   ), 2),
                   " to ",
                   round(ifelse(
                     input$interval == "hpdi",
                     hdi(df()$or, input$cred)[2],
                     qi(df()$or, input$cred)[2]
                   ), 2),
                   ")."
                 )
               }
             }
           }),
    paste
    (
      "Probability of",
      if (input$difftype == "absolute") {
        "absolute event rate difference"
      } else {
        if (input$difftype == "rr") {
          "relative risk"
        } else {
          "odds ratio"
        }
      },
      "between the groups surpassing MCID:",
      percent(over_mcid(), accuracy = 0.1)
    ),
    sep = "<br/>"
  )))
  
  output$heat <- # heatmap output
    renderPlotly(
      plot_ly(
        heatmap(),
        x =  ~ Mu,
        y =  ~ Sigma,
        z =  ~ Probability,
        type = "heatmap",
        colorscale='Reds'
        )
      )
    
  output$dt <-
    renderDataTable(datatable(sheet, selection = "single", escape = FALSE)) # datatable output and selection single
  output$dtm <-
    renderDataTable(datatable(
      sheet,
      selection = list(mode = "multiple", target = "row"),
      escape = FALSE
    )) # datatable output and selection multiple
  output$multiplot <- renderPlot(
    mplot(),
    height = function() {
      length(input$dtm_rows_selected) * 100
    }
  )
  output$multipriorplot <- renderPlot(multipriorplot())
  outputOptions(output, "dt", suspendWhenHidden = TRUE) # datatable output options
  # outputOptions(output, "posteriorplot", suspendWhenHidden = FALSE)
  # outputOptions(output, "diffplot", suspendWhenHidden = FALSE)
  # outputOptions(output, "exact", suspendWhenHidden = FALSE)
}



# Run the application
shinyApp(ui = ui, server = server)
