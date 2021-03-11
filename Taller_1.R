rm(list = ls()) ## Borra los objetos anteriores

library(ggplot2)
library(cowplot)
library(plotly)
library(dplyr)

#### physical constants ####
K_B <- 1.3806488e-23 # Boltzmann constant
h <- 6.62607015e-34
#### Microstates function ####
Microstates <- function(x, N, Thermodynamic_Limit = F){
  if(Thermodynamic_Limit == F){
    Omega <- factorial(x+N-1)/(factorial(x)*factorial(N-1))
  } else {
    Omega <- factorial(x+N)/(factorial(x)*factorial(N))
  }
  return(Omega)
}
#### Entropy function ####
Entropy <- function(N, E, Nu, Thermodynamic_Limit = F){
  x <- E/(N*h*Nu)
  if(Thermodynamic_Limit == F){
    S <- K_B*N*((x+1)*log(x+1)-x*log(x))
  } else {
    S <- K_B*N*(1+log(1+x))
  }
  return(S)
}
#### Temperature ####
Temperature <- function(N, E, Nu, Thermodynamic_Limit = F){
  x <- E/(N*h*Nu)
  if(Thermodynamic_Limit == F){
    Temp <- (h*Nu/K_B)/log((1+x)/x)
  } else {
    Temp <- E/(N*K_B)
  }
  return(Temp)
}
#### Plot Microstate's number ####
Omega <- function(x, N){
  Data <- data.frame()
  Matrix1 <- data.frame()
  Matrix2 <- data.frame()
  for(i in c(1:length(x))){
    Aux <- data.frame(
      "x" = x[i], "N" = N,
      "M_1" = Microstates(x = x[i], N = N, Thermodynamic_Limit = F),
      "M_2" = Microstates(x = x[i], N = N, Thermodynamic_Limit = T))
    Data <- rbind(Data, Aux)
    if(i == 1){
      Matrix1 <- Microstates(x = x[i], N = N, Thermodynamic_Limit = F)
      Matrix2 <- Microstates(x = x[i], N = N, Thermodynamic_Limit = T)
    } else{
      Matrix1 <- 
        cbind(Matrix1, Microstates(x = x[i], N = N, Thermodynamic_Limit = F))
      Matrix2 <- 
        cbind(Matrix2, Microstates(x = x[i], N = N, Thermodynamic_Limit = T))
    }
  }
  
  fig1 <- plot_ly(x = unique(Data$x), y = Data$N, z = ~Matrix1, 
                  color = 'Caso no asintótico',
                  legendgroup = 'Caso no asintótico') %>% 
    add_surface(
      contours = list(
        z = list(
          show=TRUE,
          usecolormap=TRUE,
          highlightcolor="#ff0000",
          project=list(z=TRUE)
        )
      ),
      colorscale = list(c(0, 1), c("green", "Firebrick2")),
      name = 'Caso no asintótico', showscale = F
    ) %>%
    layout(
      scene = list(
        xaxis = list(nticks = 20, title = "Número de paquetes de energía x"),
        yaxis = list(nticks = 20, title = "Número de partículas N"),
        zaxis = list(nticks = 10, title = "Número de microestados", type = 'log',
                     showexponent = "all", exponentformat = "e")
      ),
      legend = list(
        title = 'Caso no asintótico'
      ),
      showlegend = TRUE
    )
  
  fig2 <- plot_ly(x = unique(Data$x), y = Data$N, z = ~Matrix2, 
                  color = 'Caso asintótico',
                  legendgroup = 'Caso asintótico') %>% add_surface(
    contours = list(
      z = list(
        show=TRUE,
        usecolormap=TRUE,
        highlightcolor="#ff0000",
        project=list(z=TRUE)
      )
    ),
    colorscale = list(c(0, 1), c("midnightblue", "darkorange")),
    name = 'Caso asintótico'
  ) %>%
    layout(
      scene = list(
        xaxis = list(nticks = 20, title = "Número de paquetes de energía x"),
        yaxis = list(nticks = 20, title = "Número de partículas N"),
        zaxis = list(nticks = 10, title = "Número de microestados", type = 'log',
                     showexponent = "all", exponentformat = "e")
      ),
      legend = list(
        title = 'Caso asintótico'
      ),
      showlegend = TRUE
    )
  
  ## Update menu component
  updatemenus <- list(
    list(
      active = -1,
      type= 'buttons',
      buttons = list(
        list(
          label = "Caso no asintótico",
          method = "update",
          args = list(list(visible = c(TRUE, FALSE)),
                      list(title = TeX("\\text{Número de microestados no asintótico }\\Omega_{0}(N,E)\\text{ en escala semilogarítmica}")))),
        list(
          label = "Caso asintótico",
          method = "update",
          args = list(list(visible = c(FALSE, TRUE)),
                      list(title = TeX("\\text{Número de microestados asintótico }\\Omega_{asymp}(N,E)\\text{ en escala semilogarítmica}")))),
        list(
          label = "Todos los casos",
          method = "update",
          args = list(list(visible = c(TRUE, TRUE)),
                      list(title = TeX("\\text{Número de microestados }\\Omega(N,E)\\text{ en escala semilogarítmica}")))),
        list(
          label = "Reset",
          method = "update",
          args = list(list(visible = c(TRUE, TRUE)),
                      list(title = TeX("\\text{Número de microestados }\\Omega(N,E)\\text{ en escala semilogarítmica}"),
                           annotations = list(c(), c())))))
    )
  )
  
  fig <- subplot(fig1, fig2)
  fig <- fig %>% layout(
    title = TeX("\\text{Número de microestados }\\Omega(N,E)\\text{ en escala semilogarítmica}"),
    updatemenus=updatemenus)
  fig <- fig %>% config(mathjax = 'cdn')
  
  MyList <- list("Data" = Data, "Omega_Normal" = Matrix1, 
                 "Omega_asymp" = Matrix2, "Figure" = fig)
  
  return(MyList)
}

N <- seq(0.64 , 64, by = 0.64)
x <- seq(0.4 , 40, by = 0.4)
OmegaW <- Omega(x = x, N = N)
OmegaW$Figure
htmlwidgets::saveWidget(as_widget(OmegaW$Figure), "Taller_1_Microestados.html")

#### Plot Entropy and temperature ####

S_T <- function(N, E, Nu_vector){
  Data <- data.frame()
  MatrixList1 <- list()
  MatrixList2 <- list()
  MatrixList3 <- list()
  MatrixList4 <- list()
  for(j in c(1:length(Nu_vector))){
    Matrix1 <- data.frame()
    Matrix2 <- data.frame()
    Matrix3 <- data.frame()
    Matrix4 <- data.frame()
    for(i in c(1:length(N))){
      Aux1 <- data.frame(
        "N" = N[i], "E" = E,
        "S" = Entropy(N = N[i], E = E, Nu = Nu_vector[j], 
                        Thermodynamic_Limit = F),
        "T" = Temperature(N = N[i], E = E, Nu = Nu_vector[j], 
                            Thermodynamic_Limit = F),
        "Frequency" = paste0("Caso no asintótico ", Nu_vector[j]))
      Aux1 <- data.frame(
        "N" = N[i], "E" = E,
        "S" = Entropy(N = N[i], E = E, Nu = Nu_vector[j], 
                      Thermodynamic_Limit = F),
        "T" = Temperature(N = N[i], E = E, Nu = Nu_vector[j], 
                          Thermodynamic_Limit = F),
        "Frequency" = Nu_vector[j],
        "Type" = "Caso no asintótico")
      Aux2 <- data.frame(
        "N" = N[i], "E" = E,
        "S" = Entropy(N = N[i], E = E, Nu = Nu_vector[j], 
                      Thermodynamic_Limit = T),
        "T" = Temperature(N = N[i], E = E, Nu = Nu_vector[j], 
                          Thermodynamic_Limit = T),
        "Frequency" = Nu_vector[j],
        "Type" = "Caso asintótico")
      Data <- rbind(Data, Aux1)
      Data <- rbind(Data, Aux2)
      if(i == 1){
        Matrix1 <- Entropy(N = N[i], E = E, Nu = Nu_vector[1], 
                           Thermodynamic_Limit = F)
        Matrix2 <- Entropy(N = N[i], E = E, Nu = Nu_vector[1], 
                           Thermodynamic_Limit = T)
        Matrix3 <- Temperature(N = N[i], E = E, Nu = Nu_vector[1], 
                               Thermodynamic_Limit = F)
        Matrix4 <- Temperature(N = N[i], E = E, Nu = Nu_vector[1], 
                               Thermodynamic_Limit = T)
      } else{
        Matrix1 <- 
          cbind(Matrix1, 
                Entropy(N = N[i], E = E, Nu = Nu_vector[1], 
                        Thermodynamic_Limit = F))
        Matrix2 <- 
          cbind(Matrix2, 
                Entropy(N = N[i], E = E, Nu = Nu_vector[1], 
                        Thermodynamic_Limit = T))
        Matrix3 <- 
          cbind(Matrix3, 
                Temperature(N = N[i], E = E, Nu = Nu_vector[1], 
                            Thermodynamic_Limit = F))
        Matrix4 <- 
          cbind(Matrix4, 
                Temperature(N = N[i], E = E, Nu = Nu_vector[1], 
                            Thermodynamic_Limit = T))
      }
    }
    MatrixList1 <- list.append(MatrixList1, Matrix1)
    MatrixList2 <- list.append(MatrixList2, Matrix2)
    MatrixList3 <- list.append(MatrixList3, Matrix3)
    MatrixList4 <- list.append(MatrixList4, Matrix4)
  }
  
  fig1 <- plot_ly(x = unique(Data$N), y = Data$E, 
                  color = 'Caso no asintótico',
                  z = ~MatrixList1[[1]],
                  legendgroup = 'Caso no asintótico') %>% 
    add_surface(
      contours = list(
        z = list(
          show=TRUE,
          usecolormap=TRUE,
          highlightcolor=list(c(0, 1, 1), c("midnightblue", "Steelblue", "white")),
          project=list(z=TRUE)
        )
      ),
      colorscale = list(c(0, 1, 1), c("midnightblue", "Steelblue", "white")),
      name = 'Caso no asintótico', 
      showscale = F
    ) %>%
    layout(
      scene = list(
        xaxis = list(nticks = 20, title = "Número de partículas N"),
        yaxis = list(nticks = 20, title = "Energía E (J)"),
        zaxis = list(nticks = 10, title = "Entropía S (J/K)", type = 'linear',
                     showexponent = "all", exponentformat = "e")
      ),
      legend = list(title = 'Caso no asintótico'),
      showlegend = TRUE
    )
  
  fig2 <- plot_ly(x = unique(Data$N), y = Data$E, 
                  color = 'Caso asintótico',
                  z = ~MatrixList2[[1]],
                  legendgroup = 'Caso asintótico') %>% 
    add_surface(
      contours = list(
        z = list(
          show=TRUE,
          usecolormap=TRUE,
          highlightcolor = list(c(1, 1, 1), c("darkorange", "yellow", "red")),
          project=list(z=TRUE)
        )
      ),
      colorscale = list(c(1, 1 , 1), c("darkorange", "yellow", "red")),
      name = 'Caso asintótico', 
      showscale = F
    ) %>%
    layout(
      scene = list(
        xaxis = list(nticks = 20, title = "Número de partículas N"),
        yaxis = list(nticks = 20, title = "Energía E (J)"),
        zaxis = list(nticks = 10, title = "Entropía S (J/K)", type = 'linear',
                     showexponent = "all", exponentformat = "e")
      ),
      legend = list(title = 'Caso asintótico'),
      showlegend = TRUE
    )
  
  Names <- paste0(Data$Type, " - Freq: ", as.character(Data$Frequency))
  
  fig3 <- plot_ly(x = ~Data$E, y = ~Data$S, frame = ~Data$N, 
                  size = ~Data$Frequency,
                  color = ~Data$Type,
                  name = ~Names,
                  type = 'scatter', mode = 'markers')
  
  ## Update menu component
  updatemenus <- list(
    list(
      active = -1,
      type= 'buttons',
      buttons = list(
        list(
          label = "Caso no asintótico",
          method = "update",
          args = list(list(visible = rep(c(TRUE, FALSE),length(Nu_vector))),
                      list(title = TeX("\\text{Entropía no asintótica }S_{0}(N,E), \\nu=400 Hz")))),
        list(
          label = "Caso asintótico",
          method = "update",
          args = list(list(visible = rep(c(FALSE, TRUE),length(Nu_vector))),
                      list(title = TeX("\\text{Entropía asintótica }S_{asymp}(N,E), \\nu=400 Hz")))),
        list(
          label = "Todos los casos",
          method = "update",
          args = list(list(visible = rep(c(TRUE, TRUE),length(Nu_vector))),
                      list(title = TeX("\\text{Entropía }S(N,E), \\nu=400 Hz")))),
        list(
          label = "Reset",
          method = "update",
          args = list(list(visible = rep(c(TRUE, TRUE),length(Nu_vector))),
                      list(title = TeX("\\text{Entropía }S(N,E), \\nu=400 Hz"),
                           annotations = list(c(), c())))))
    )
  )
  
  fig_T1 <- subplot(fig1,fig2)
  fig_T1 <- fig_T1 %>% layout(
    title = TeX("\\text{Entropía }S(N,E)"),
    updatemenus=updatemenus)
  fig_T1 <- fig_T1 %>% config(mathjax = 'cdn')
  
  fig_T2 <- fig3
  fig_T2 <- fig_T2 %>% layout(
    title = TeX("\\text{Entropía }S(N,E)\\text{ variando la frecuencia }\\nu"),
    xaxis = list(title = "Energía E (J)", nticks = 20),
    yaxis = list(title = "Entropía S (J/K)", nticks = 10, type = 'linear',
                 showexponent = "all", exponentformat = "e")) %>%
    animation_opts(transition = 0, easing = "elastic", redraw = FALSE) %>%
    animation_slider(
      currentvalue = list(prefix = "Número de partículas N: ")
    )
  fig_T2 <- fig_T2 %>% config(mathjax = 'cdn')
  
  MyList <- list("Data" = Data, 
                 "S_Normal" = MatrixList1, "S_asymp" = MatrixList2, 
                 "T_Normal" = MatrixList3, "T_asymp" = MatrixList4, 
                 "Figure3D" = fig_T1, "Figure2D" = fig_T2)
  
  return(MyList)
}

N <- seq(0.5e30 , 1e31, by = 0.5e30)
E <- seq(1 , 100, by = 1)
Nu <- c(400, 1600, 8000)
EntropyW <- S_T(N = N, E = E, Nu_vector = Nu)
EntropyW$Figure3D
EntropyW$Figure2D
htmlwidgets::saveWidget(as_widget(EntropyW$Figure3D), "Taller_1_Entropía3D.html")
htmlwidgets::saveWidget(as_widget(EntropyW$Figure2D), "Taller_1_Entropía2D.html")

