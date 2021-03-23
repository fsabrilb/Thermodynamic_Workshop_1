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
#### Helmholtz Energy function ####
Helmholtz <- function(N, E, Nu, Thermodynamic_Limit = F){
  x <- E/(N*h*Nu)
  if(Thermodynamic_Limit == F){
    A <- -(E/x)*log(1+x)/(log(1+x)-log(x))
  } else {
    A <- -E*log(1+x)
  }
  return(A)
}
#### Gibbs Energy function ####
Gibbs <- function(N, E, Nu, Thermodynamic_Limit = F){
  G <- Helmholtz(N = N, E = E, Nu = Nu, Thermodynamic_Limit = Thermodynamic_Limit)
  return(G)
}
#### Enthalpy function ####
Enthalpy <- function(N, E, Nu, Thermodynamic_Limit = F){
  H <- Helmholtz(N = N, E = E, Nu = Nu, Thermodynamic_Limit = Thermodynamic_Limit)
  return(H)
}
#### Specific heat capacity function ####
Specific_Heat <- function(N, E, Nu, Thermodynamic_Limit = F){
  x <- E/(N*h*Nu)
  if(Thermodynamic_Limit == F){
    Cv <- N*K_B*x*(1+x)*(log(1+1/x))^2
  } else {
    Cv <- N*K_B
  }
  return(Cv)
}
#### Temperature function ####
Temperature <- function(N, E, Nu, Thermodynamic_Limit = F){
  x <- E/(N*h*Nu)
  if(Thermodynamic_Limit == F){
    Temp <- (h*Nu/K_B)/log((1+x)/x)
  } else {
    Temp <- E/(N*K_B)
  }
  return(Temp)
}
#### Plot microstate's number ####
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

#### Plot state's functions ####
States_Function <- function(N, E, Nu_vector){
  #### Data ####
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
        "A" = Helmholtz(N = N[i], E = E, Nu = Nu_vector[j], 
                        Thermodynamic_Limit = F),
        "Frequency" = Nu_vector[j],
        "Type" = "Caso no asintótico")
      Aux2 <- data.frame(
        "N" = N[i], "E" = E,
        "S" = Entropy(N = N[i], E = E, Nu = Nu_vector[j], 
                      Thermodynamic_Limit = T),
        "A" = Helmholtz(N = N[i], E = E, Nu = Nu_vector[j], 
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
        Matrix3 <- Helmholtz(N = N[i], E = E, Nu = Nu_vector[1], 
                             Thermodynamic_Limit = F)
        Matrix4 <- Helmholtz(N = N[i], E = E, Nu = Nu_vector[1], 
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
                Helmholtz(N = N[i], E = E, Nu = Nu_vector[1], 
                          Thermodynamic_Limit = F))
        Matrix4 <- 
          cbind(Matrix4, 
                Helmholtz(N = N[i], E = E, Nu = Nu_vector[1], 
                          Thermodynamic_Limit = T))
      }
    }
    MatrixList1 <- list.append(MatrixList1, Matrix1)
    MatrixList2 <- list.append(MatrixList2, Matrix2)
    MatrixList3 <- list.append(MatrixList3, Matrix3)
    MatrixList4 <- list.append(MatrixList4, Matrix4)
  }
  #### Plot 1 Entropy ####
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
  #### Plot 2 Entropy ####
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
  #### Plot 1 Helmholtz energy ####
  fig3 <- plot_ly(x = unique(Data$N), y = Data$E, 
                  color = 'Caso no asintótico',
                  z = ~MatrixList3[[1]], scene = 'scene1',
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
      scene1 = list(
        xaxis = list(nticks = 20, title = "Número de partículas N"),
        yaxis = list(nticks = 20, title = "Energía E (J)"),
        zaxis = list(nticks = 10, title = "Energía libre de Helmholtz A (J)", 
                     type = 'linear',
                     showexponent = "all", exponentformat = "e")
      ),
      legend = list(title = 'Caso no asintótico'),
      showlegend = TRUE
    )
  #### Plot 2 Helmholtz energy ####
  fig4 <- plot_ly(x = unique(Data$N), y = Data$E, 
                  color = 'Caso asintótico', scene = 'scene1',
                  z = ~MatrixList4[[1]],
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
      scene1 = list(
        xaxis = list(nticks = 20, title = "Número de partículas N"),
        yaxis = list(nticks = 20, title = "Energía E (J)"),
        zaxis = list(nticks = 10, title = "Energía libre de Helmholtz A (J)", 
                     type = 'linear', 
                     showexponent = "all", exponentformat = "e")
      ),
      legend = list(title = 'Caso asintótico'),
      showlegend = TRUE
    )
  #### Plot 3 Entropy (Contour level) ####
  Names <- paste0(Data$Type, " - Freq: ", as.character(Data$Frequency))
  
  fig5 <- plot_ly(x = ~Data$E, y = ~Data$S, frame = ~Data$N, 
                  size = ~Data$Frequency,
                  color = ~Data$Type,
                  name = ~Names,
                  type = 'scatter', mode = 'markers')
  #### Plot 3 Helmholtz energy (Contour level) ####
  fig6 <- plot_ly(x = ~Data$E, y = ~Data$A, frame = ~Data$N, 
                  size = ~Data$Frequency,
                  color = ~Data$Type,
                  name = ~Names,
                  type = 'scatter', mode = 'markers')
  
  #### Update menu component Entropy ####
  updatemenus <- list(
    list(
      active = -1,
      type= 'buttons',
      buttons = list(
        list(
          label = "Caso no asintótico",
          method = "update",
          args = list(list(visible = c(TRUE, FALSE)),
                      list(title = TeX("\\text{Entropía no asintótica }S_{0}(N,E), \\nu=400 Hz")))),
        list(
          label = "Caso asintótico",
          method = "update",
          args = list(list(visible = c(FALSE, TRUE)),
                      list(title = TeX("\\text{Entropía asintótica }S_{asymp}(N,E), \\nu=400 Hz")))),
        list(
          label = "Todos los casos",
          method = "update",
          args = list(list(visible = c(TRUE, TRUE)),
                      list(title = TeX("\\text{Entropía }S(N,E), \\nu=400 Hz")))),
        list(
          label = "Reset",
          method = "update",
          args = list(list(visible = c(TRUE, TRUE)),
                      list(title = TeX("\\text{Entropía }S(N,E), \\nu=400 Hz"),
                           annotations = list(c(), c())))))
    )
  )
  #### Update menu component Helmholtz energy ####
  updatemenus1 <- list(
    list(
      active = -1,
      type= 'buttons',
      buttons = list(
        list(
          label = "Caso no asintótico",
          method = "update",
          args = list(list(visible = c(TRUE, FALSE)),
                      list(title = TeX("\\text{Energía libre de Helmholtz no asintótica }A_{0}(N,E), \\nu=400 Hz")))),
        list(
          label = "Caso asintótico",
          method = "update",
          args = list(list(visible = c(FALSE, TRUE)),
                      list(title = TeX("\\text{Energía libre de Helmholtz asintótica }A_{asymp}(N,E), \\nu=400 Hz")))),
        list(
          label = "Todos los casos",
          method = "update",
          args = list(list(visible = c(TRUE, TRUE)),
                      list(title = TeX("\\text{Energía libre de Helmholtz }A(N,E), \\nu=400 Hz")))),
        list(
          label = "Reset",
          method = "update",
          args = list(list(visible = c(TRUE, TRUE)),
                      list(title = TeX("\\text{Energía libre de Helmholtz }A(N,E), \\nu=400 Hz"),
                           annotations = list(c(), c())))))
    )
  )
  #### Final plots ####
  ## Entropy
  fig_T1 <- subplot(fig1,fig2)
  fig_T1 <- fig_T1 %>% layout(
    title = TeX("\\text{Entropía }S(N,E)"),
    updatemenus=updatemenus)
  fig_T1 <- fig_T1 %>% config(mathjax = 'cdn')
  ## Helmholtz energy
  fig_T2 <- subplot(fig3,fig4)
  fig_T2 <- fig_T2 %>% layout(
    title = TeX("\\text{Energía libre de Helmholtz }A(N,E)"),
    updatemenus=updatemenus1)
  fig_T2 <- fig_T2 %>% config(mathjax = 'cdn')
  ## Entropy contour level
  fig_T3 <- fig5
  fig_T3 <- fig_T3 %>% layout(
    title = TeX("\\text{Entropía }S(N,E)\\text{ variando la frecuencia }\\nu"),
    xaxis = list(title = "Energía E (J)", nticks = 20),
    yaxis = list(title = "Entropía S (J/K)", nticks = 10, type = 'linear',
                 showexponent = "all", exponentformat = "e")) %>%
    animation_opts(transition = 0, easing = "elastic", redraw = FALSE) %>%
    animation_slider(
      currentvalue = list(prefix = "Número de partículas N: ")
    )
  fig_T3 <- fig_T3 %>% config(mathjax = 'cdn')
  ## Helmholtz energy contour level
  fig_T4 <- fig6
  fig_T4 <- fig_T4 %>% layout(
    title = TeX("\\text{Energía libre de Helmholtz }A(N,E)\\text{ variando la frecuencia }\\nu"),
    xaxis = list(title = "Energía E (J)", nticks = 20),
    yaxis = list(title = "Energía libre de Helmholtz A(J)", nticks = 10, type = 'linear',
                 showexponent = "all", exponentformat = "e")) %>%
    animation_opts(transition = 0, easing = "elastic", redraw = FALSE) %>%
    animation_slider(
      currentvalue = list(prefix = "Número de partículas N: ")
    )
  fig_T4 <- fig_T4 %>% config(mathjax = 'cdn')
  #### Function return ####
  MyList <- list("Data" = Data, 
                 "S_Normal" = MatrixList1, "S_asymp" = MatrixList2, 
                 "A_Normal" = MatrixList3, "A_asymp" = MatrixList4, 
                 "Figure3D_S" = fig_T1, "Figure2D_S" = fig_T3,
                 "Figure3D_A" = fig_T2, "Figure2D_A" = fig_T4)
  
  return(MyList)
}

N <- seq(0.5e30 , 1e31, by = 0.5e30)
E <- seq(1 , 100, by = 1)
Nu <- c(400, 1600, 8000)
StatesW <- States_Function(N = N, E = E, Nu_vector = Nu)
StatesW$Figure3D_S
StatesW$Figure2D_S
StatesW$Figure3D_A
StatesW$Figure2D_A
htmlwidgets::saveWidget(as_widget(StatesW$Figure3D_S), "Taller_1_Entropía3D.html")
htmlwidgets::saveWidget(as_widget(StatesW$Figure2D_S), "Taller_1_Entropía2D.html")
htmlwidgets::saveWidget(as_widget(StatesW$Figure3D_A), "Taller_1_Helmholtz3D.html")
htmlwidgets::saveWidget(as_widget(StatesW$Figure2D_A), "Taller_1_Helmholtz2D.html")

#### Plot temperature and specific heat ####
Thermodynamic_Measuring <- function(N, E, Nu_vector){
  #### Data ####
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
        "Te" = Temperature(N = N[i], E = E, Nu = Nu_vector[j], 
                          Thermodynamic_Limit = F),
        "C" = Specific_Heat(N = N[i], E = E, Nu = Nu_vector[j], 
                            Thermodynamic_Limit = F),
        "Frequency" = Nu_vector[j],
        "Type" = "Caso no asintótico")
      Aux2 <- data.frame(
        "N" = N[i], "E" = E,
        "Te" = Temperature(N = N[i], E = E, Nu = Nu_vector[j], 
                          Thermodynamic_Limit = T),
        "C" = Specific_Heat(N = N[i], E = E, Nu = Nu_vector[j], 
                            Thermodynamic_Limit = T),
        "Frequency" = Nu_vector[j],
        "Type" = "Caso asintótico")
      Data <- rbind(Data, Aux1)
      Data <- rbind(Data, Aux2)
      if(i == 1){
        Matrix1 <- Temperature(N = N[i], E = E, Nu = Nu_vector[1], 
                               Thermodynamic_Limit = F)
        Matrix2 <- Temperature(N = N[i], E = E, Nu = Nu_vector[1], 
                               Thermodynamic_Limit = T)
        Matrix3 <- Specific_Heat(N = N[i], E = E, Nu = Nu_vector[1], 
                                 Thermodynamic_Limit = F)
        Matrix4 <- Specific_Heat(N = N[i], E = E, Nu = Nu_vector[1], 
                                 Thermodynamic_Limit = T)
      } else{
        Matrix1 <- 
          cbind(Matrix1, 
                Temperature(N = N[i], E = E, Nu = Nu_vector[1], 
                            Thermodynamic_Limit = F))
        Matrix2 <- 
          cbind(Matrix2, 
                Temperature(N = N[i], E = E, Nu = Nu_vector[1], 
                            Thermodynamic_Limit = T))
        Matrix3 <- 
          cbind(Matrix3, 
                Specific_Heat(N = N[i], E = E, Nu = Nu_vector[1], 
                              Thermodynamic_Limit = F))
        Matrix4 <- 
          cbind(Matrix4, 
                Specific_Heat(N = N[i], E = E, Nu = Nu_vector[1], 
                              Thermodynamic_Limit = T))
      }
    }
    MatrixList1 <- list.append(MatrixList1, Matrix1)
    MatrixList2 <- list.append(MatrixList2, Matrix2)
    MatrixList3 <- list.append(MatrixList3, Matrix3)
    MatrixList4 <- list.append(MatrixList4, Matrix4)
  }
  #### Plot 1 Temperature ####
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
        zaxis = list(nticks = 10, title = "Temperatura T(K)", type = 'linear',
                     showexponent = "all", exponentformat = "e")
      ),
      legend = list(title = 'Caso no asintótico'),
      showlegend = TRUE
    )
  #### Plot 2 Temperature ####
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
        zaxis = list(nticks = 10, title = "Temperatura T(K)", type = 'linear',
                     showexponent = "all", exponentformat = "e")
      ),
      legend = list(title = 'Caso asintótico'),
      showlegend = TRUE
    )
  #### Plot 1 Specific heat ####
  fig3 <- plot_ly(x = unique(Data$N), y = Data$E, 
                  color = 'Caso no asintótico',
                  z = ~MatrixList3[[1]], scene = 'scene1',
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
      scene1 = list(
        xaxis = list(nticks = 20, title = "Número de partículas N"),
        yaxis = list(nticks = 20, title = "Energía E (J)"),
        zaxis = list(nticks = 10, title = "Capacidad calorífica específica C(J/K)", 
                     type = 'linear',
                     showexponent = "all", exponentformat = "e")
      ),
      legend = list(title = 'Caso no asintótico'),
      showlegend = TRUE
    )
  #### Plot 2 Specific heat ####
  fig4 <- plot_ly(x = unique(Data$N), y = Data$E, 
                  color = 'Caso asintótico', scene = 'scene1',
                  z = ~MatrixList4[[1]],
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
      scene1 = list(
        xaxis = list(nticks = 20, title = "Número de partículas N"),
        yaxis = list(nticks = 20, title = "Energía E (J)"),
        zaxis = list(nticks = 10, title = "Capacidad calorífica específica C(J/K)", 
                     type = 'linear', 
                     showexponent = "all", exponentformat = "e")
      ),
      legend = list(title = 'Caso asintótico'),
      showlegend = TRUE
    )
  #### Plot 3 Temperature (Contour level) ####
  Names <- paste0(Data$Type, " - Freq: ", as.character(Data$Frequency))
  
  fig5 <- plot_ly(x = ~Data$E, y = ~Data$Te, frame = ~Data$N, 
                  size = ~Data$Frequency,
                  color = ~Data$Type,
                  name = ~Names,
                  type = 'scatter', mode = 'markers')
  #### Plot 3 Specific heat (Contour level) ####
  fig6 <- plot_ly(x = ~Data$E, y = ~Data$C, frame = ~Data$N, 
                  size = ~Data$Frequency,
                  color = ~Data$Type,
                  name = ~Names,
                  type = 'scatter', mode = 'markers')
  #### Update menu component Temperature ####
  updatemenus <- list(
    list(
      active = -1,
      type= 'buttons',
      buttons = list(
        list(
          label = "Caso no asintótico",
          method = "update",
          args = list(list(visible = c(TRUE, FALSE)),
                      list(title = TeX("\\text{Temperatura no asintótica }T_{0}(N,E), \\nu=400 Hz")))),
        list(
          label = "Caso asintótico",
          method = "update",
          args = list(list(visible = c(FALSE, TRUE)),
                      list(title = TeX("\\text{Temperatura asintótica }T_{asymp}(N,E), \\nu=400 Hz")))),
        list(
          label = "Todos los casos",
          method = "update",
          args = list(list(visible = c(TRUE, TRUE)),
                      list(title = TeX("\\text{Temperatura }T(N,E), \\nu=400 Hz")))),
        list(
          label = "Reset",
          method = "update",
          args = list(list(visible = c(TRUE, TRUE)),
                      list(title = TeX("\\text{Temperatura }T(N,E), \\nu=400 Hz"),
                           annotations = list(c(), c())))))
    )
  )
  #### Update menu component Specific heat ####
  updatemenus1 <- list(
    list(
      active = -1,
      type= 'buttons',
      buttons = list(
        list(
          label = "Caso no asintótico",
          method = "update",
          args = list(list(visible = c(TRUE, FALSE)),
                      list(title = TeX("\\text{Capacidad calorífica específica no asintótica }C_{V,0}(N,E), \\nu=400 Hz")))),
        list(
          label = "Caso asintótico",
          method = "update",
          args = list(list(visible = c(FALSE, TRUE)),
                      list(title = TeX("\\text{Capacidad calorífica específica asintótica }C_{V,asymp}(N,E), \\nu=400 Hz")))),
        list(
          label = "Todos los casos",
          method = "update",
          args = list(list(visible = c(TRUE, TRUE)),
                      list(title = TeX("\\text{Capacidad calorífica específica }C_{V}(N,E), \\nu=400 Hz")))),
        list(
          label = "Reset",
          method = "update",
          args = list(list(visible = c(TRUE, TRUE)),
                      list(title = TeX("\\text{Capacidad calorífica específica }C_{V}(N,E), \\nu=400 Hz"),
                           annotations = list(c(), c())))))
    )
  )
  #### Final plots ####
  ## Temperature
  fig_T1 <- subplot(fig1,fig2)
  fig_T1 <- fig_T1 %>% layout(
    title = TeX("\\text{Temperatura }T(N,E)"),
    updatemenus=updatemenus)
  fig_T1 <- fig_T1 %>% config(mathjax = 'cdn')
  ## Specific heat energy
  fig_T2 <- subplot(fig3,fig4)
  fig_T2 <- fig_T2 %>% layout(
    title = TeX("\\text{Capacidad calorífica específica }C_{V}(N,E)"),
    updatemenus=updatemenus1)
  fig_T2 <- fig_T2 %>% config(mathjax = 'cdn')
  ## Temperature contour level
  fig_T3 <- fig5
  fig_T3 <- fig_T3 %>% layout(
    title = TeX("\\text{Temperatura }T(N,E)\\text{ variando la frecuencia }\\nu"),
    xaxis = list(title = "Energía E (J)", nticks = 20),
    yaxis = list(title = "Temperatura T(K)", nticks = 10, type = 'linear',
                 showexponent = "all", exponentformat = "e")) %>%
    animation_opts(transition = 0, easing = "elastic", redraw = FALSE) %>%
    animation_slider(
      currentvalue = list(prefix = "Número de partículas N: ")
    )
  fig_T3 <- fig_T3 %>% config(mathjax = 'cdn')
  ## Specific heat contour level
  fig_T4 <- fig6
  fig_T4 <- fig_T4 %>% layout(
    title = TeX("\\text{Capacidad calorífica específica }C_{V}(N,E)\\text{ variando la frecuencia }\\nu"),
    xaxis = list(title = "Energía E (J)", nticks = 20),
    yaxis = list(title = "Capacidad calorífica específica C(J/K)", nticks = 10, type = 'linear',
                 showexponent = "all", exponentformat = "e")) %>%
    animation_opts(transition = 0, easing = "elastic", redraw = FALSE) %>%
    animation_slider(
      currentvalue = list(prefix = "Número de partículas N: ")
    )
  fig_T4 <- fig_T4 %>% config(mathjax = 'cdn')
  #### Function return ####
  MyList <- list("Data" = Data, 
                 "T_Normal" = MatrixList1, "T_asymp" = MatrixList2, 
                 "C_Normal" = MatrixList3, "C_asymp" = MatrixList4, 
                 "Figure3D_T" = fig_T1, "Figure2D_T" = fig_T3,
                 "Figure3D_C" = fig_T2, "Figure2D_C" = fig_T4)
  
  return(MyList)
}

N <- seq(0.5e30 , 1e31, by = 0.5e30)
E <- seq(1 , 100, by = 1)
Nu <- c(400, 1600, 8000)
ThermoW <- Thermodynamic_Measuring(N = N, E = E, Nu_vector = Nu)
ThermoW$Figure3D_T
ThermoW$Figure2D_T
ThermoW$Figure3D_C
ThermoW$Figure2D_C
htmlwidgets::saveWidget(as_widget(ThermoW$Figure3D_T), "Taller_1_Temperatura3D.html")
htmlwidgets::saveWidget(as_widget(ThermoW$Figure2D_T), "Taller_1_Temperatura2D.html")
htmlwidgets::saveWidget(as_widget(ThermoW$Figure3D_C), "Taller_1_CalorEspecifico3D.html")
htmlwidgets::saveWidget(as_widget(ThermoW$Figure2D_C), "Taller_1_CalorEspecifico2D.html")
