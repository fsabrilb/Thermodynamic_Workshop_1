rm(list = ls()) ## Borra los objetos anteriores

library(ggplot2)
library(cowplot)
library(plotly)
library(dplyr)

#### physical constants ####
K_B <- 1.3806488e-23 # Boltzmann constant
h <- 6.62607015e-34
c <- 299792458
#### Microstates function ####
Microstates <- function(N, L, E){
  x <- 2*L*E/(h*c)
  Omega <- x^N/(factorial(N)*factorial(N))
  return(Omega)
}
N <- seq(1e0, 1e1, by = 1e0)
E <- seq(1e-12, 1e-11, by = 1e-12)
L <- seq(1e-12, 1e-11, by = 1e-12)

A <- Microstates(N = 3*N, L = L, E = E)
#### Entropy function ####
Entropy <- function(N, L, E){
  x <- 2*L*E/(h*c)
  S <- K_B*N*(2+log(x/(N*N)))
  return(S)
}
N <- seq(1e25, 1e26, by = 1e25)
E <- seq(1e13, 1e14, by = 1e13)
L <- seq(1e13, 1e14, by = 1e13)

A <- Entropy(N = 3*N, L = L, E = E)
#### Helmholtz Energy function ####
Helmholtz <- function(N, L, Temp){
  Lambda <- h*c/(K_B*Temp)
  A <- -K_B*N*Temp*(1+log(2*L/(N*Lambda)))
  return(A)
}
N <- seq(1e25, 1e26, by = 1e25)
L <- seq(1e13, 1e14, by = 1e13)
Te <- seq(1e7, 1e8, by = 1e7)

A <- Helmholtz(N = 3*N, L = L, Temp = Te)
#### Gibbs Energy function ####
Gibbs <- function(N, P, Temp){
  Lambda <- h*c/(K_B*Temp)
  G <- -K_B*N*Temp*log(2*K_B*Temp/(P*Lambda)) 
  return(G)
}
N <- seq(1e25, 1e26, by = 1e25)
P <- seq(1e13, 1e14, by = 1e13)
Te <- seq(1e7, 1e8, by = 1e7)

A <- Gibbs(N = 3*N, P = P, Temp = Te)
#### Enthalpy function ####
Enthalpy <- function(N, P, S){
  H <- N*sqrt(2*h*c*P)*exp(S/(6*N*K_B)-1) 
  return(H)
}
N <- seq(1e22, 1e23, by = 1e22)
P <- seq(1e3, 1e4, by = 1e3)
S <- seq(1e1, 1e2, by = 1e1)

A <- Enthalpy(N = 3*N, P = P, S = S)
#### Specific heat capacity at constant volume (pressure) function ####
Specific_Heat_V <- function(N){
  CV <- N*K_B
  return(CV)
}
Specific_Heat_P <- function(N){
  CP <- 2*N*K_B
  return(CP)
}
N <- seq(1e22, 1e23, by = 1e22)
Gamma <- Specific_Heat_P(N = 3*N)/Specific_Heat_V(N = 3*N)
#### Temperature function ####
Temperature <- function(N, L, E){
  Temp <- E/(N*K_B)
  return(Temp)
}
N <- seq(1e25, 1e26, by = 1e25)
L <- seq(1e13, 1e14, by = 1e13)
E <- seq(1e13, 1e14, by = 1e13)

A <- Temperature(N = 3*N, L = L, E = E)
#### Pressure function ####
Pressure <- function(N, L, E){
  P <- E/L
  return(P)
}
N <- seq(1e25, 1e26, by = 1e25)
L <- seq(1e13, 1e14, by = 1e13)
E <- seq(1e13, 1e14, by = 1e13)

A <- Pressure(N = 3*N, L = L, E = E)
#### Chemical Potential function ####
ChemicalPotential <- function(N, L, E){
  x <- 2*L*E/(h*c)
  Mu <- -(E/N)*log(x/(N*N))
  return(Mu)
}
N <- seq(1e25, 1e26, by = 1e25)
L <- seq(1e13, 1e14, by = 1e13)
E <- seq(1e18, 1e19, by = 1e18)

A <- ChemicalPotential(N = 3*N, L = L, E = E)
#### Plot intensive properties T, P, Mu ####
IntensiveProperties <- function(N, L, E){
  #### Data ####
  Data <- data.frame()
  Matrix1 <- data.frame()
  Matrix2 <- data.frame()
  Matrix3 <- data.frame()
  for(i in c(1:length(E))){
    for(j in c(1:length(L))){
      Aux <- data.frame(
        "E" = E[i], "L" = L[j], "N" = N,
        "Temp" = Temperature(N = N, L = L[j], E = E[i]),
        "P" = Pressure(N = N, L = L[j], E = E[i]),
        "Mu" = ChemicalPotential(N = N, L = L[j], E = E[i]))
      Data <- rbind(Data, Aux)
    }
    if(i == 1){
      Matrix1 <- Temperature(N = N, L = L[1], E = E[i])
      Matrix2 <- Pressure(N = N, L = L[1], E = E[i])
      Matrix3 <- ChemicalPotential(N = N, L = L[1], E = E[i])
    } else{
      Matrix1 <- cbind(Matrix1, Temperature(N = N, L = L[1], E = E[i]))
      Matrix2 <- cbind(Matrix2, Pressure(N = N, L = L[1], E = E[i]))
      Matrix3 <- cbind(Matrix3, ChemicalPotential(N = N, L = L[1], E = E[i]))
    }
  }
  #### Plot Temperature ####
  fig1 <- plot_ly(x = unique(Data$N), y = unique(Data$E), z = ~Matrix1, 
                  color = 'Temperatura', scene = 'scene1',
                  legendgroup = 'Temperatura') %>% 
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
      name = 'Temperatura', showscale = F
    ) %>%
    layout(
      scene1 = list(
        xaxis = list(nticks = 20, title = "Número de partículas N",
                     showexponent = "all", exponentformat = "e"),
        yaxis = list(nticks = 20, title = "Energía E (J)",
                     showexponent = "all", exponentformat = "e"),
        zaxis = list(nticks = 10, title = "Temperatura T (K)", type = 'linear',
                     showexponent = "all", exponentformat = "e")
      ),
      legend = list(
        title = 'Temperatura'
      ),
      showlegend = TRUE
    )
  #### Plot Pressure ####
  fig2 <- plot_ly(x = unique(Data$N), y = unique(Data$E), z = ~Matrix2, 
                  color = 'Presión', scene = 'scene2',
                  legendgroup = 'Presión') %>% 
    add_surface(
      contours = list(
        z = list(
          show=TRUE,
          usecolormap=TRUE,
          highlightcolor="#ff0000",
          project=list(z=TRUE)
        )
      ),
      colorscale = list(c(0, 1), c("midnightblue", "darkorange")),
      name = 'Presión', showscale = F
    ) %>%
    layout(
      scene2 = list(
        xaxis = list(nticks = 20, title = "Número de partículas N",
                     showexponent = "all", exponentformat = "e"),
        yaxis = list(nticks = 20, title = "Energía E (J)",
                     showexponent = "all", exponentformat = "e"),
        zaxis = list(nticks = 10, title = "Presión longitudinal P (J/m)", type = 'linear',
                     showexponent = "all", exponentformat = "e")
      ),
      legend = list(
        title = 'Presión'
      ),
      showlegend = TRUE
    )
  #### Plot Chemical Potential ####
  fig3 <- plot_ly(x = unique(Data$N), y = unique(Data$E), z = ~Matrix3, 
                  color = 'Potencial químico', scene = 'scene3',
                  legendgroup = 'Potencial químico') %>% 
    add_surface(
      contours = list(
        z = list(
          show=TRUE,
          usecolormap=TRUE,
          highlightcolor="#ff0000",
          project=list(z=TRUE)
        )
      ),
      colorscale = list(c(0, 1), c("darkgoldenrod1", "darkmagenta")),
      name = 'Potencial químico', showscale = F
    ) %>%
    layout(
      scene3 = list(
        xaxis = list(nticks = 20, title = "Número de partículas N",
                     showexponent = "all", exponentformat = "e"),
        yaxis = list(nticks = 20, title = "Energía E (J)",
                     showexponent = "all", exponentformat = "e"),
        zaxis = list(nticks = 10, title = "Potencial químico (J)", type = 'linear',
                     showexponent = "all", exponentformat = "e")
      ),
      legend = list(
        title = 'Potencial químico'
      ),
      showlegend = TRUE
    )
  #### Plot Temperature (Contour level) ####
  Names <- paste0("L: ", as.character(Data$L), " m")
  
  fig4 <- plot_ly(x = ~Data$E, y = ~Data$Temp, frame = ~Data$N, 
                  color = ~Data$L,
                  name = ~Names,
                  type = 'scatter', mode = 'markers')
  #### Plot Pressure (Contour level) ####
  fig5 <- plot_ly(x = ~Data$E, y = ~Data$P, frame = ~Data$N, 
                  color = ~Data$L,
                  name = ~Names,
                  type = 'scatter', mode = 'markers')
  #### Plot Chemical Potential (Contour level) ####
  fig6 <- plot_ly(x = ~Data$E, y = ~Data$Mu, frame = ~Data$N, 
                  color = ~Data$L,
                  name = ~Names,
                  type = 'scatter', mode = 'markers')
  #### Update menu component ####
  updatemenus <- list(
    list(
      active = -1,
      type= 'buttons',
      buttons = list(
        list(
          label = "Temperatura",
          method = "update",
          args = list(list(visible = c(TRUE, FALSE)),
                      list(title = TeX("\\text{Temperatura gas ultrarelativista unidimensional }T(N,L,E)\\text{ a longitud constante}")))),
        list(
          label = "Potencial químico",
          method = "update",
          args = list(list(visible = c(FALSE, TRUE)),
                      list(title = TeX("\\text{Potencial químico gas ultrarelativista unidimensional }\\mu(N,L,E)\\text{ a longitud constante}")))),
        list(
          label = "Reset",
          method = "update",
          args = list(list(visible = c(TRUE, TRUE)),
                      list(title = paste0('Propiedades intensivas gas ultrarelativista unidimensional a longitud constante L = ', L[1], " m"),
                           annotations = list(c(), c())))))
    )
  )
  #### Final plots ####
  fig <- subplot(fig1, fig3)
  fig <- fig %>% layout(
    title = paste0('Propiedades intensivas gas ultrarelativista unidimensional a longitud constante L = ', L[1], " m"),
    scene = list(
      xaxis = list(nticks = 20, title = "Número de partículas N",
                   showexponent = "all", exponentformat = "e"),
      yaxis = list(nticks = 20, title = "Energía E (J)",
                   showexponent = "all", exponentformat = "e"),
      zaxis = list(nticks = 10, title = "Temperatura T (K)", type = 'linear',
                   showexponent = "all", exponentformat = "e")
    ),
    legend = list(
      title = 'Variables intensivas'
    ),
    showlegend = TRUE,
    updatemenus=updatemenus)
  fig <- fig %>% config(mathjax = 'cdn')
  ## Temperature contour level
  fig_T4 <- fig4
  fig_T4 <- fig_T4 %>% layout(
    title = paste0("Temperatura gas ultrarelativista unidimensional T(N,L,E) variando la longitud"),
    xaxis = list(title = "Energía E (J)", nticks = 20,
                 showexponent = "all", exponentformat = "e"),
    yaxis = list(title = "Temperatura T (K)", nticks = 10, type = 'linear',
                 showexponent = "all", exponentformat = "e"),
    legend = list(
      title = 'Longitud'
    ), showlegend = TRUE) %>%
    animation_opts(transition = 0, easing = "elastic", redraw = FALSE) %>%
    animation_slider(
      currentvalue = list(prefix = "Número de partículas N: ")
    )
  fig_T4 <- fig_T4 %>% config(mathjax = 'cdn')
  ## Pressure contour level
  fig_T5 <- fig5
  fig_T5 <- fig_T5 %>% layout(
    title = paste0("Presión gas ultrarelativista unidimensional P(N,L,E) variando la longitud"),
    xaxis = list(title = "Energía E (J)", nticks = 20,
                 showexponent = "all", exponentformat = "e"),
    yaxis = list(title = "Presión longitudinal P (J/m)", nticks = 10, type = 'linear',
                 showexponent = "all", exponentformat = "e")) %>%
    animation_opts(transition = 0, easing = "elastic", redraw = FALSE) %>%
    animation_slider(
      currentvalue = list(prefix = "Número de partículas N: ")
    )
  fig_T5 <- fig_T5 %>% config(mathjax = 'cdn')
  ##  Chemical Potential contour level
  fig_T6 <- fig6
  fig_T6 <- fig_T6 %>% layout(
    title = paste0("Potencial químico gas ultrarelativista unidimensional variando la longitud"),
    xaxis = list(title = "Energía E (J)", nticks = 20,
                 showexponent = "all", exponentformat = "e"),
    yaxis = list(title = "Potencial químico (J)", nticks = 10, type = 'linear',
                 showexponent = "all", exponentformat = "e")) %>%
    animation_opts(transition = 0, easing = "elastic", redraw = FALSE) %>%
    animation_slider(
      currentvalue = list(prefix = "Número de partículas N: ")
    )
  fig_T6 <- fig_T6 %>% config(mathjax = 'cdn')
  #### Function return ####
  MyList <- list("Data" = Data, "Temperature" = Matrix1, "Pressure" = Matrix2,
                 "ChemicalPotential" = Matrix3, "Figure" = fig,
                 "Figure3D_T" = fig1, "Figure3D_P" = fig2, "Figure3D_Mu" = fig3,
                 "FigureT" = fig_T4, "FigureP" = fig_T5, "FigureMu" = fig_T6)
  
  return(MyList)
}

N <- seq(0.5e25, 1e26, by = 0.5e25)
L <- seq(2e13, 1e14, by = 2e13)
E <- seq(0.1e15, 1e16, by = 0.1e15)

IP <- IntensiveProperties(N = 3*N, L = L, E = E)
IP$Figure
IP$Figure3D_T
IP$Figure3D_Mu
IP$FigureT
IP$FigureP
IP$FigureMu
htmlwidgets::saveWidget(as_widget(IP$Figure), "Taller 2/Taller_2_PropiedadesIntensivas.html")
htmlwidgets::saveWidget(as_widget(IP$Figure3D_T), "Taller 2/Taller_2_Temperatura_3D.html")
htmlwidgets::saveWidget(as_widget(IP$FigureT), "Taller 2/Taller_2_Temperatura_2D.html")
htmlwidgets::saveWidget(as_widget(IP$FigureP), "Taller 2/Taller_2_Presion_2D.html")
htmlwidgets::saveWidget(as_widget(IP$Figure3D_Mu), "Taller 2/Taller_2_PotencialQuimico_3D.html")
htmlwidgets::saveWidget(as_widget(IP$FigureMu), "Taller 2/Taller_2_PotencialQuimico_2D.html")

#### Plot state's functions ####
States_Function <- function(N, L, E, S, P, Temp){
  #### Data ####
  Data1 <- data.frame()
  Data2 <- data.frame()
  Data3 <- data.frame()
  Data4 <- data.frame()
  Matrix1 <- data.frame()
  Matrix2 <- data.frame()
  Matrix3 <- data.frame()
  Matrix4 <- data.frame()
  for(i in c(1:length(N))){
    ## Entropy and Helmholtz energy
    for(j in c(1:length(L))){
      Aux1 <- data.frame(
        "E_Temp" = E, "L" = L[j], "N" = N[i],
        "S_A" = Entropy(N = N[i], L = L[j], E = E))
      Aux2 <- data.frame(
        "E_Temp" = Temp, "L" = L[j], "N" = N[i],
        "S_A" = Helmholtz(N = N[i], L = L[j], Temp = Temp))
      Data1 <- rbind(Data1, Aux1)
      Data2 <- rbind(Data2, Aux2)
    }
    if(i == 1){
      Matrix1 <- Entropy(N = N[i], L = L[1], E = E)
      Matrix2 <- Helmholtz(N = N[i], L = L[1], Temp = Temp)
    } else{
      Matrix1 <- cbind(Matrix1, Entropy(N = N[i], L = L[1], E = E))
      Matrix2 <- cbind(Matrix2, Helmholtz(N = N[i], L = L[1], Temp = Temp))
    }
    ## Enthalpy and Gibbs energy
    for(j in c(1:length(P))){
      Aux3 <- data.frame(
        "S_Temp" = S, "P" = P[j], "N" = N[i],
        "H_G" = Enthalpy(N = N[i], P = P[j], S = S))
      Aux4 <- data.frame(
        "S_Temp" = Temp, "P" = P[j], "N" = N[i],
        "H_G" = Gibbs(N = N[i], P = P[j], Temp = Temp))
      Data3 <- rbind(Data3, Aux3)
      Data4 <- rbind(Data4, Aux4)
    }
    if(i == 1){
      Matrix3 <- Enthalpy(N = N[i], P = P[j], S = S)
      Matrix4 <- Gibbs(N = N[i], P = P[j], Temp = Temp)
    } else{
      Matrix3 <- cbind(Matrix3, Enthalpy(N = N[i], P = P[j], S = S))
      Matrix4 <- cbind(Matrix4, Gibbs(N = N[i], P = P[j], Temp = Temp))
    }
  }
  #### Plot Entropy ####
  fig1 <- plot_ly(x = unique(Data1$N), y = unique(Data1$E_Temp), z = ~Matrix1, 
                  color = 'Entropía', scene = 'scene1',
                  legendgroup = 'Entropía') %>% 
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
      name = 'Entropía', showscale = F
    ) %>%
    layout(
      scene1 = list(
        xaxis = list(nticks = 20, title = "Número de partículas N",
                     showexponent = "all", exponentformat = "e"),
        yaxis = list(nticks = 20, title = "Energía E (J)",
                     showexponent = "all", exponentformat = "e"),
        zaxis = list(nticks = 10, title = "Entropía S (J/K)", type = 'linear',
                     showexponent = "all", exponentformat = "e")
      ),
      legend = list(
        title = 'Entropía'
      ),
      showlegend = TRUE
    )
  #### Plot Helmholtz Energy ####
  fig2 <- plot_ly(x = unique(Data2$N), y = unique(Data2$E_Temp), z = ~Matrix2, 
                  color = 'Energía de Helmholtz', scene = 'scene2',
                  legendgroup = 'Energía de Helmholtz') %>% 
    add_surface(
      contours = list(
        z = list(
          show=TRUE,
          usecolormap=TRUE,
          highlightcolor="#ff0000",
          project=list(z=TRUE)
        )
      ),
      colorscale = list(c(0, 1), c("midnightblue", "darkorange")),
      name = 'Energía de Helmholtz', showscale = F
    ) %>%
    layout(
      scene2 = list(
        xaxis = list(nticks = 20, title = "Número de partículas N",
                     showexponent = "all", exponentformat = "e"),
        yaxis = list(nticks = 20, title = "Temperatura (K)",
                     showexponent = "all", exponentformat = "e"),
        zaxis = list(nticks = 10, title = "Energía de Helmholtz A (J)", type = 'linear',
                     showexponent = "all", exponentformat = "e")
      ),
      legend = list(
        title = 'Energía de Helmholtz'
      ),
      showlegend = TRUE
    )
  #### Plot Enthalpy ####
  fig3 <- plot_ly(x = unique(Data3$N), y = unique(Data3$S_Temp), z = ~Matrix3, 
                  color = 'Entalpía', scene = 'scene3',
                  legendgroup = 'Entalpía') %>% 
    add_surface(
      contours = list(
        z = list(
          show=TRUE,
          usecolormap=TRUE,
          highlightcolor="#ff0000",
          project=list(z=TRUE)
        )
      ),
      colorscale = list(c(0, 1), c("darkgoldenrod1", "darkmagenta")),
      name = 'Entalpía', showscale = F
    ) %>%
    layout(
      scene3 = list(
        xaxis = list(nticks = 20, title = "Número de partículas N",
                     showexponent = "all", exponentformat = "e"),
        yaxis = list(nticks = 20, title = "Entropía S (J/K)",
                     showexponent = "all", exponentformat = "e"),
        zaxis = list(nticks = 10, title = "Entalpía H (J)", type = 'linear',
                     showexponent = "all", exponentformat = "e")
      ),
      legend = list(
        title = 'Entalpía'
      ),
      showlegend = TRUE
    )
  #### Plot Gibbs Energy ####
  fig4 <- plot_ly(x = unique(Data4$N), y = unique(Data4$S_Temp), z = ~Matrix4, 
                  color = 'Energía de Gibbs', scene = 'scene4',
                  legendgroup = 'Energía de Gibbs') %>% 
    add_surface(
      contours = list(
        z = list(
          show=TRUE,
          usecolormap=TRUE,
          highlightcolor="#ff0000",
          project=list(z=TRUE)
        )
      ),
      colorscale = list(c(0, 1), c("darkgoldenrod1", "darkmagenta")),
      name = 'Energía de Gibbs', showscale = F
    ) %>%
    layout(
      scene4 = list(
        xaxis = list(nticks = 20, title = "Número de partículas N",
                     showexponent = "all", exponentformat = "e"),
        yaxis = list(nticks = 20, title = "Temperatura T (K)",
                     showexponent = "all", exponentformat = "e"),
        zaxis = list(nticks = 10, title = "Energía de Gibbs G (J)", type = 'linear',
                     showexponent = "all", exponentformat = "e")
      ),
      legend = list(
        title = 'Energía de Gibbs'
      ),
      showlegend = TRUE
    )
  #### Plot Entropy (Contour level) ####
  Names1 <- paste0("L: ", as.character(Data1$L), " m")
  
  fig5 <- plot_ly(x = ~Data1$E_Temp, y = ~Data1$S_A, frame = ~Data1$N, 
                  color = ~Data1$L,
                  name = ~Names1,
                  type = 'scatter', mode = 'markers')
  #### Plot Helmholtz Energy (Contour level) ####
  Names2 <- paste0("L: ", as.character(Data2$L), " m")
  
  fig6 <- plot_ly(x = ~Data2$E_Temp, y = ~Data2$S_A, frame = ~Data2$N, 
                  color = ~Data2$L,
                  name = ~Names2,
                  type = 'scatter', mode = 'markers')
  #### Plot Enthalpy (Contour level) ####
  Names3 <- paste0("Presión longitudinal: ", as.character(Data3$P), " J/m")
  
  fig7 <- plot_ly(x = ~Data3$S_Temp, y = ~Data3$H_G, frame = ~Data3$N, 
                  color = ~Data3$P,
                  name = ~Names3,
                  type = 'scatter', mode = 'markers')
  #### Plot Gibbs Energy (Contour level) ####
  Names4 <- paste0("Presión longitudinal: ", as.character(Data4$P), " J/m")
  
  fig8 <- plot_ly(x = ~Data4$S_Temp, y = ~Data4$H_G, frame = ~Data4$N, 
                  color = ~Data4$P,
                  name = ~Names4,
                  type = 'scatter', mode = 'markers')
  #### Update menu component surfaces ####
  updatemenus <- list(
    list(
      active = -1,
      type= 'buttons',
      buttons = list(
        list(
          label = "Entropía",
          method = "update",
          args = list(list(visible = c(TRUE, FALSE, FALSE, FALSE)),
                      list(title = TeX("\\text{Entropía gas ultrarelativista unidimensional }S(N,L,E)\\text{ a longitud constante}")))),
        list(
          label = "Energía libre de Helmholtz",
          method = "update",
          args = list(list(visible = c(FALSE, TRUE, FALSE, FALSE)),
                      list(title = TeX("\\text{Energía libre de Helmholtz gas ultrarelativista unidimensional }A(N,L,T)\\text{ a longitud constante}")))),
        list(
          label = "Entalpía",
          method = "update",
          args = list(list(visible = c(FALSE, FALSE, TRUE, FALSE)),
                      list(title = TeX("\\text{Entalpía gas ultrarelativista unidimensional }H(N,P,S)\\text{ a presión constante}")))),
        list(
          label = "Energía libre de Gibbs",
          method = "update",
          args = list(list(visible = c(FALSE, FALSE, FALSE, TRUE)),
                      list(title = TeX("\\text{Energía libre de Gibbs gas ultrarelativista unidimensional }G(N,P,T)\\text{ a presión constante}")))),
        list(
          label = "Reset",
          method = "update",
          args = list(list(visible = c(TRUE, TRUE, TRUE, TRUE)),
                      list(title = paste0('Funciones de estado gas ultrarelativista unidimensional'),
                           annotations = list(c(), c())))))
    )
  )
  #### Final plots ####
  fig <- subplot(fig1, fig2, fig3, fig4, nrows = 2#, 
                 #titleX = TRUE, titleY = TRUE,
                 #shareX = TRUE, shareY = TRUE,
                 #heights = c(0.5, 0.5), widths = c(0.5, 0.5)
                 )
  fig <- fig %>% layout(
    title = paste0('Funciones de estado gas ultrarelativista unidimensional'),
    scene = list(
      domain=list(x=c(0,0.5),y=c(0.5,1)),
      xaxis = list(nticks = 20, title = "Número de partículas N",
                   showexponent = "all", exponentformat = "e"),
      yaxis = list(nticks = 20, title = "Energía E (J)",
                   showexponent = "all", exponentformat = "e"),
      zaxis = list(nticks = 10, title = "Entropía S (J/K)", type = 'linear',
                   showexponent = "all", exponentformat = "e")
    ),
    scene2 = list(
      domain=list(x=c(0.5,1),y=c(0.5,1)),
      xaxis = list(nticks = 20, title = "Número de partículas N",
                   showexponent = "all", exponentformat = "e"),
      yaxis = list(nticks = 20, title = "Temperatura (K)",
                   showexponent = "all", exponentformat = "e"),
      zaxis = list(nticks = 10, title = "Energía de Helmholtz A (J)", type = 'linear',
                   showexponent = "all", exponentformat = "e")
    ),
    scene3 = list(
      domain=list(x=c(0,0.5),y=c(0,0.5)),
      xaxis = list(nticks = 20, title = "Número de partículas N",
                   showexponent = "all", exponentformat = "e"),
      yaxis = list(nticks = 20, title = "Entropía S (J/K)",
                   showexponent = "all", exponentformat = "e"),
      zaxis = list(nticks = 10, title = "Entalpía H (J)", type = 'linear',
                   showexponent = "all", exponentformat = "e")
    ),
    scene4 = list(
      domain=list(x=c(0.5,1),y=c(0,0.5)),
      xaxis = list(nticks = 20, title = "Número de partículas N",
                   showexponent = "all", exponentformat = "e"),
      yaxis = list(nticks = 20, title = "Temperatura T (K)",
                   showexponent = "all", exponentformat = "e"),
      zaxis = list(nticks = 10, title = "Energía de Gibbs G (J)", type = 'linear',
                   showexponent = "all", exponentformat = "e")
    ),
    updatemenus=updatemenus)
  fig <- fig %>% config(mathjax = 'cdn')
  ## Entropy contour level
  fig_T5 <- fig5
  fig_T5 <- fig_T5 %>% layout(
    title = paste0("Entropía gas ultrarelativista unidimensional S(N,L,E) variando la longitud"),
    xaxis = list(title = "Energía E (J)", nticks = 20,
                 showexponent = "all", exponentformat = "e"),
    yaxis = list(title = "Entropía S (J/K)", nticks = 10, type = 'linear',
                 showexponent = "all", exponentformat = "e")) %>%
    animation_opts(transition = 0, easing = "elastic", redraw = FALSE) %>%
    animation_slider(
      currentvalue = list(prefix = "Número de partículas N: ")
    )
  fig_T5 <- fig_T5 %>% config(mathjax = 'cdn')
  ## Helmholtz energy contour level
  fig_T6 <- fig6
  fig_T6 <- fig_T6 %>% layout(
    title = paste0("Energía de Helmholtz gas ultrarelativista unidimensional A(N,L,T) variando la longitud"),
    xaxis = list(title = "Temperatura T (K)", nticks = 20,
                 showexponent = "all", exponentformat = "e"),
    yaxis = list(title = "Energía libre de Helmholtz A (J)", nticks = 10, type = 'linear',
                 showexponent = "all", exponentformat = "e")) %>%
    animation_opts(transition = 0, easing = "elastic", redraw = FALSE) %>%
    animation_slider(
      currentvalue = list(prefix = "Número de partículas N: ")
    )
  fig_T6 <- fig_T6 %>% config(mathjax = 'cdn')
  ##  Enthalpy contour level
  fig_T7 <- fig7
  fig_T7 <- fig_T7 %>% layout(
    title = paste0("Entalpía gas ultrarelativista unidimensional H(N,P,S) variando la presión"),
    xaxis = list(title = "Entropía S (J/K)", nticks = 20,
                 showexponent = "all", exponentformat = "e"),
    yaxis = list(title = "Entalpía H (J)", nticks = 10, type = 'linear',
                 showexponent = "all", exponentformat = "e")) %>%
    animation_opts(transition = 0, easing = "elastic", redraw = FALSE) %>%
    animation_slider(
      currentvalue = list(prefix = "Número de partículas N: ")
    )
  fig_T7 <- fig_T7 %>% config(mathjax = 'cdn')
  ##  Gibbs energy contour level
  fig_T8 <- fig8
  fig_T8 <- fig_T8 %>% layout(
    title = paste0("Energía libre de Gibbs gas ultrarelativista unidimensional G(N,P,T) variando la presión"),
    xaxis = list(title = "Temperatura T (K)", nticks = 20,
                 showexponent = "all", exponentformat = "e"),
    yaxis = list(title = "Energía libre de Gibbs G (J)", nticks = 10, type = 'linear',
                 showexponent = "all", exponentformat = "e")) %>%
    animation_opts(transition = 0, easing = "elastic", redraw = FALSE) %>%
    animation_slider(
      currentvalue = list(prefix = "Número de partículas N: ")
    )
  fig_T8 <- fig_T8 %>% config(mathjax = 'cdn')
  #### Function return ####
  MyList <- list("Entropy_Data" = Data1, "Helmholtz_Data" = Data2, 
                 "Enthalpy_Data" = Data3, "Gibbs_Data" = Data4,
                 "Entropy" = Matrix1, "Helmholtz" = Matrix2,
                 "Enthalpy" = Matrix3, "Gibbs" = Matrix4,
                 "Figure" = fig,
                 "Figure3D_S" = fig1, "Figure3D_A" = fig2, 
                 "Figure3D_H" = fig3, "Figure3D_G" = fig4,
                 "Figure2D_S" = fig_T5, "Figure2D_A" = fig_T6, 
                 "Figure2D_H" = fig_T7, "Figure2D_G" = fig_T8)
  
  return(MyList)
}

N <- seq(0.5e25, 1e26, by = 0.5e25)
L <- seq(1e13, 1e14, by = 1e13)
E <- seq(0.1e15, 1e16, by = 0.1e15)
S <- seq(0.1e1, 1e2, by = 0.1e1)
P <- seq(1e3, 1e4, by = 1e3)
Te <- seq(0.5e7, 1e8, by = 0.5e7)

StatesW <- States_Function(N = 3*N, L = L, E = E, S = S, P = P, Temp = Te)
StatesW$Figure3D_S
StatesW$Figure3D_A
StatesW$Figure3D_H
StatesW$Figure3D_G

StatesW$Figure2D_S
StatesW$Figure2D_A
StatesW$Figure2D_H
StatesW$Figure2D_G

StatesW$Figure

htmlwidgets::saveWidget(as_widget(StatesW$Figure3D_S), "Taller 2/Taller_2_Entropía_3D.html")
htmlwidgets::saveWidget(as_widget(StatesW$Figure3D_A), "Taller 2/Taller_2_Helmholtz_3D.html")
htmlwidgets::saveWidget(as_widget(StatesW$Figure3D_H), "Taller 2/Taller_2_Entalpia_3D.html")
htmlwidgets::saveWidget(as_widget(StatesW$Figure3D_G), "Taller 2/Taller_2_Gibbs_3D.html")

htmlwidgets::saveWidget(as_widget(StatesW$Figure2D_S), "Taller 2/Taller_2_Entropía_2D.html")
htmlwidgets::saveWidget(as_widget(StatesW$Figure2D_A), "Taller 2/Taller_2_Helmholtz_2D.html")
htmlwidgets::saveWidget(as_widget(StatesW$Figure2D_H), "Taller 2/Taller_2_Entalpia_2D.html")
htmlwidgets::saveWidget(as_widget(StatesW$Figure2D_G), "Taller 2/Taller_2_Gibbs_2D.html")

htmlwidgets::saveWidget(as_widget(StatesW$Figure), "Taller 2/Taller_2_FuncionesEstado.html")
