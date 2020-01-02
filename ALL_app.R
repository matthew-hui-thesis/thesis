list.of.packages <- c("shiny", "shinyjs", "dplyr", "data.table", "ggplot2", "tictoc")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(shiny)
library(shinyjs)
library(dplyr)
library(data.table)
library(ggplot2)
library(tictoc)

#Input criteria
max.Ninf <- 4 #(minimum is 1)
max.Nobs <- 12 #(minimum is 2)

#Functions for EBE estimation
##Model parameters
TVCL=7.7293915062049097
TVV1=19.000870864793672
TVQ=0.28290322491080838
TVV2=6.6282022544476051
B_BSA_V1=0.98469240611826103
B_BSA_CL=0.72130221741781342
B_CLCR_CL=0.25634021485234787
B_AGE_Q=0.27798497399395100

SG=0.0910229819834100
OMCL=2.05154942269959763E-002
OMV1=0
OMQ=0
OMV2=0.11954530193107613
OMIOV=2.22935054783364207E-002

OM = diag(c(OMCL, OMV1, OMQ, OMV2))
SG = matrix(SG)

NsC1 = function(s, Rate, CL, V1, Q, V2, C10, C20) {
  V1*V2**2*C10*s**3 +
    (V1*V2*C10*Q+Rate*V2**2+V2**2*C20*Q+V1*V2*C10*Q)*s**2 +
    (2*Rate*V2*Q+V2*C20*Q**2+V1*C10*Q**2)*s +
    Rate*Q**2
}
comp2InfC1 <- function(t, Rate, popCL, etaCL, etaIOV, popV1, etaV1, popQ, etaQ, popV2, etaV2, C10, C20) {
  CL = popCL*exp(etaCL+etaIOV)
  V1 = popV1*exp(etaV1)
  Q = popQ*exp(etaQ)
  V2 = popV2*exp(etaV2)
  a = (Q/V2+Q/V1+CL/V1-sqrt((Q/V2+Q/V1+CL/V1)**2-4*CL*Q/V1/V2))/2
  b = (Q/V2+Q/V1+CL/V1+sqrt((Q/V2+Q/V1+CL/V1)**2-4*CL*Q/V1/V2))/2
  (Rate*Q/a/b-NsC1(-Q/V2, Rate, CL, V1, Q, V2, C10, C20)/(Q/V2)/(a-Q/V2)/(b-Q/V2)*exp(-Q/V2*t)+
      NsC1(-a, Rate, CL, V1, Q, V2, C10, C20)/a/(V2*a-Q)/(b-a)*exp(-a*t)+
      NsC1(-b, Rate, CL, V1, Q, V2, C10, C20)/b/(V2*b-Q)/(a-b)*exp(-b*t))/V1/V2
}
NsC2 = function(s, Rate, CL, V1, Q, V2, C10, C20) {
  V1*V2*C20*s**2 +
    (V2*Q*C20+V2*CL*C20+V1*Q*C10)*s +
    Rate*Q
}
comp2InfC2 <- function(t, Rate, popCL, etaCL, etaIOV, popV1, etaV1, popQ, etaQ, popV2, etaV2, C10, C20) {
  CL = popCL*exp(etaCL+etaIOV)
  V1 = popV1*exp(etaV1)
  Q = popQ*exp(etaQ)
  V2 = popV2*exp(etaV2)
  a = (Q/V2+Q/V1+CL/V1-sqrt((Q/V2+Q/V1+CL/V1)**2-4*CL*Q/V1/V2))/2
  b = (Q/V2+Q/V1+CL/V1+sqrt((Q/V2+Q/V1+CL/V1)**2-4*CL*Q/V1/V2))/2
  (Rate*Q/a/b+NsC2(-a, Rate, CL, V1, Q, V2, C10, C20)/a/(a-b)*exp(-a*t)+
      NsC2(-b, Rate, CL, V1, Q, V2, C10, C20)/b/(b-a)*exp(-b*t))/V1/V2
}

getComp2conc <- function(times, timeInf, Rate, arrPopCL, arrPopV1, arrPopQ, arrPopV2, etaCL, etaV1, etaQ, etaV2) {
  #Next observation carried backward is assumed
  if (!any(times %in% timeInf)&&timeInf<max(times)) {
    added <- T
    iID <- (sort(c(times, timeInf))==timeInf) %>% which
    cTimes <- sort(c(times, timeInf))
    cRate <- c(rep(Rate, iID), rep(0, length(cTimes)-iID))
    cPopCL <- c(arrPopCL[1:(iID-1)], arrPopCL[iID], arrPopCL[iID:length(arrPopCL)])
    cPopV1 <- c(arrPopV1[1:(iID-1)], arrPopV1[iID], arrPopV1[iID:length(arrPopV1)])
    cPopQ <- c(arrPopQ[1:(iID-1)], arrPopQ[iID], arrPopQ[iID:length(arrPopQ)])
    cPopV2 <- c(arrPopV2[1:(iID-1)], arrPopV2[iID], arrPopV2[iID:length(arrPopV2)])
  } else {
    added <- F
    cTimes <- times
    cRate <- sapply(cTimes, function(t) {if (t<=timeInf) {Rate} else {0}})
    cPopCL <- arrPopCL
    cPopV1 <- arrPopV1
    cPopQ <- arrPopQ
    cPopV2 <- arrPopV2
  }
  arrC1 <- 0
  arrC2 <- 0
  for (i in 2:length(cTimes)) {
    lst_args <- list(t = cTimes[i]-cTimes[i-1], Rate = cRate[i],
                     popCL = cPopCL[i], etaCL = etaCL, etaIOV = 0,
                     popV1 = cPopV1[i], etaV1 = etaV1,
                     popQ = cPopQ[i], etaQ = etaQ,
                     popV2 = cPopV2[i], etaV2 = etaV2,
                     arrC1[i-1], arrC2[i-1])
    arrC1 <- append(arrC1, do.call("comp2InfC1", lst_args))
    arrC2 <- append(arrC2, do.call("comp2InfC2", lst_args))
  }
  if (added) {arrC1 <- arrC1[-iID]; arrC2 <- arrC2[-iID]}
  return(matrix(c(arrC1, arrC1), byrow = F, ncol = 2))
}

EBE_optim <- function(ETA, DATA) {
  with(DATA, {
    etaCL = if (diag(OM)[1]>0) {ETA[1]} else {0}
    etaV1 = if (diag(OM)[2]>0) {ETA[2]} else {0}
    etaQ = if (diag(OM)[3]>0) {ETA[3]} else {0}
    etaV2 = if (diag(OM)[4]>0) {ETA[4]} else {0}
    lapply(seq(Ninf), function(i) {
      lst_args = list(
        times = TAD[[i]],
        timeInf = infTime[i],
        Rate = dosePerBSA[i]*BSAdose[i]*2200.51/infTime[i],
        arrPopCL = TVCL*(BSA[[i]]/0.735)**B_BSA_CL*(CLCR[[i]]/192)**B_CLCR_CL,
        arrPopV1 = TVV1*(BSA[[i]]/0.735)**B_BSA_V1,
        arrPopQ = TVQ*(AGE[[i]]/63.5)**B_AGE_Q,
        arrPopV2 = rep(TVV2, length(TAD[[i]])),
        etaCL = etaCL, etaV1 = etaV1, etaQ = etaQ, etaV2 = etaV2
      )
      return(do.call(getComp2conc, lst_args)[-1,])
    }) %>% do.call(rbind, .) -> FGH
    F = FGH[,1] %>% matrix(., ncol = 1)
    H = FGH[,2] %>% matrix(., ncol = 1)
    COV <- ((H) %*% SG %*% t(H)) %>% diag %>% diag(., ncol=length(F), nrow=length(F))
    Y <- do.call(c, lapply(seq(Ninf), function(i) {DV[[i]][-1]}))
    RES <- Y-F
    diagOM_on <- sapply(diag(OM), function(x) {x>0})
    nETA <- ETA[diagOM_on]
    nOM <- (diag(OM)[diagOM_on]) %>% diag(., ncol = length(which(diagOM_on)), nrow = length(which(diagOM_on)))
    OBJ <- log(det(COV)) + t(RES) %*% solve(COV) %*% RES + t(nETA) %*% solve(nOM) %*% nETA
    return(OBJ %>% as.vector)
  })
}

MTX_UI <- fluidPage(
  shinyjs::useShinyjs(),
  #Hidden panel for example button use
    conditionalPanel("input.Ninf > 999.99", fluidRow(disabled(numericInput("egDummy", NULL, 0)))),
  #Hidden panel for result panel use
    conditionalPanel("input.Ninf > 999.99", fluidRow(disabled(numericInput("resDummy", NULL, 0)))),
  #Title
    fluidRow(column(6, titlePanel("Methotrexate dose calculator for Acute Lymphoblastic Leukaemia"))),
    hr(),
  #Tab formatting
    tags$style(type = 'text/css', HTML(".container-fluid > .nav > li > a[data-value='tabOpt'] {background-color: #339966; color:white;}
                                       .container-fluid > .nav > li > a:hover[data-value='tabOpt'] {background-color: #228855; color:white;}
                                       .container-fluid > .nav > li > a:focus[data-value='tabOpt'] {background-color: #117744; color:white;}")),
  #tabs
  navbarPage(title = "", id = "selTab",
   tabPanel("Step 1 - Getting individual parameters", icon = icon("user", "fa-2x"), value = "tab1", wellPanel(
     #Getting individual parameters
     ##Number of infusions
     fluidRow(numericInput("Ninf", "Total number of infusions to input:", 1,
                           min = 1, max = max.Ninf, step = 1) %>% column(., width = 3),
              column(3, actionButton("exampleButton", "Example values (ID=1)"))),
     ##Use renderUI to produce the UI for the max.Ndose doses
     tags$div(id = "InfPanel")
   )),
   tabPanel("Step 2 - Setting therapeutic target(s)", icon = icon("bullseye", "fa-2x"), value = "tab2", wellPanel(
     ##Adjust therapeutic target
     fluidRow(column(12, textOutput("crteriaDetail"))),
     fluidRow(titlePanel(" ")),
     fluidRow(column(12, textOutput("tut_str_def_Cpss"))),
     fluidRow(column(3, numericInput("def_Cpss_LB", NULL, 16, min = 0, step = 1)))
   )),
   tabPanel("Step 3 - Current demographics", icon = icon("sliders", "fa-2x"), value = "tab3", wellPanel(
     #Info for this dose
     fluidRow(column(6,titlePanel("Current demographics"))),
     wellPanel(
       fluidRow(column(2, textInput("nowBSA", "BSA (m2):", 1)),
                column(2, textInput("nowCLCR", "CLCR (mL/min/1.73m2):", 1)),
                column(2, textInput("nowAGE", "AGE (months):", 1)))
     )
   )),
   tabPanel("DOSE OPTIMIZATION", icon = icon("cogs", "fa-2x"), value = "tabOpt", wellPanel(
     #Calculation button
     fluidRow(column(actionButton("calculate", "Calculate", style="color: #ffffff; background-color: #228855; border-color: #11b45d; font-size: 18px"), width = 2)),
     fluidRow(column(6, textOutput("msgCal"))),
     br(),
     br(),
     #Result part
     conditionalPanel("input.resDummy == 1",
                      wellPanel(fluidRow(column(6, titlePanel("Results"))),
                                fluidRow(column(12, textOutput("optDPA"))),
                                fluidRow(column(12, textOutput("optDose"))),
                                hr(),
                                fluidRow(column(6, plotOutput("plot_1")), column(6, plotOutput("plot_2"))),
                                fluidRow(column(6, textOutput("Legend_1")))
                      )
     ),
     #Styling
     tags$head(tags$style("#optDPA{font-size: 16px}")),
     tags$head(tags$style("#optDose{font-size: 16px}")),
     tags$head(tags$style("#crteriaDetail{font-size: 18px}")),
     tags$head(tags$style("#tut_str_def_Cpss{font-size: 16px}"))
   )))
)

MTX_Server <- function(input, output, session) {
#Initiation
  withProgress(
    {incProgress(0, message = "Loading... just a few seconds...")
    },
    message = "Calculator setup"
  )
  output$crteriaDetail <- renderText(paste0("Minimum dose is defined as the dose at which there would be a 95% chance of achieving steady-state plasma drug level of at least ", input$def_Cpss_LB, " uM."))
  output$tut_str_def_Cpss <- renderText("Or you may define the minimum steady-state plasma drug level target yourself:")
#UI creation
  change_Ninf <- reactive({input[["Ninf"]]})
  ExistNinf <<- 0
  
  #Creating inf panels
    observeEvent(change_Ninf(), {
      Ninf <- change_Ninf()
      if (Ninf>ExistNinf) {    
        #Decide list to insert
          if (ExistNinf==0) {lstToIns <- seq(Ninf)} else {lstToIns <- seq(Ninf)[-seq(ExistNinf)]}
        #Update list inserted
          ExistNinf <<- Ninf
        #Insert UI & reactive component
          lapply(lstToIns,
                 function(InfN) {
            #Insert UI
              insertUI(
                selector = '#InfPanel',
                ui = conditionalPanel(paste0("input.Ninf >= ", InfN),
                                      hr(),
                                      wellPanel(
                                        style = if (InfN%%2==0) {"background-color: #ccf5fe; border-color: #44535c"} else {"background-color: #aed7f1; border-color: #44535c"},
                                        fluidRow(column(3, titlePanel(paste0("Dose #", InfN)))),
                                        fluidRow(column(3, numericInput(paste0("Nobs_", InfN),
                                                                        "Number of observations (including start of infusion):", 2, min = 2, max = max.Nobs, step = 1)),
                                                 column(3, textInput(paste0("dosePerBSA_", InfN), "Total dose per BSA (g/m2):", 5)),
                                                 column(3, textInput(paste0("BSA_", InfN), "BSA (m2):", 0.47)),
                                                 column(3, textInput(paste0("infTime_", InfN), "Infusion time (hr):", 24))
                                        ),
                                        fluidRow(column(2, "Time after infusion started (hr):"),
                                                 column(2, "Conc (umol/L):"),
                                                 column(2, "BSA (m2):"),
                                                 column(2, "CLCR (mL/min/1.73m2):"),
                                                 column(2, "AGE (months):")
                                        ),
                                        fluidRow(
                                          column(2, offset = 4, checkboxInput(paste0("fix_BSA_", InfN), "Same", FALSE)),
                                          column(2, checkboxInput(paste0("fix_CLCR_", InfN), "Same", FALSE)),
                                          column(2, checkboxInput(paste0("fix_AGE_", InfN), "Same", FALSE))
                                        ),
                                        tags$div(id = paste0("ObsPanel_", InfN))
                                      ))
              )
            #Update Obs reactive components
              assign(paste0("change_Nobs_", InfN), reactive({input[[paste0("Nobs_", InfN)]]}), inherits = T)
            #Creating Obs components
              if (!exists(paste0("ExistNobs_", InfN))) {assign(paste0("ExistNobs_", InfN), 0)}
              observeEvent(get(paste0("change_Nobs_", InfN))(), {
                Nobs <- get(paste0("change_Nobs_", InfN))()
                ExistNobs <- get(paste0("ExistNobs_", InfN))
                lstToIns_obs <- numeric()
                if (Nobs>ExistNobs) {
                  #Decide list to insert
                    if (ExistNobs==0) {lstToIns_obs <- seq(Nobs)} else {lstToIns_obs <- seq(Nobs)[-seq(ExistNobs)]}
                  #Update list inserted
                    assign(paste0("ExistNobs_", InfN), Nobs, inherits = T)
                  #Insert UI
                    lapply(lstToIns_obs, function(ObsN) {
                      if (ObsN==1) {
                        insertUI(
                          selector = paste0("#ObsPanel_", InfN),
                          ui = conditionalPanel(paste0("input.Nobs_", InfN, " >= ", 1),
                                                fluidRow(
                                                  column(2, disabled(textInput(paste0("TAD_", InfN, "_", 1), NULL, 0))),
                                                  column(2, disabled(textInput(paste0("DV_", InfN, "_", 1), NULL, 0))),
                                                  column(2, textInput(paste0("BSA_", InfN, "_", 1), NULL, 1)),
                                                  column(2, textInput(paste0("CLCR_", InfN, "_", 1), NULL, 1)),
                                                  column(2, textInput(paste0("AGE_", InfN, "_", 1), NULL, 1))
                                                )),
                          session = session
                        )
                      } else {
                        insertUI(
                          selector = paste0("#ObsPanel_", InfN),
                          ui = conditionalPanel(paste0("input.Nobs_", InfN, " >= ", ObsN),
                                                fluidRow(
                                                  column(2, textInput(paste0("TAD_", InfN, "_", ObsN), NULL, 24*(ObsN-1))),
                                                  column(2, textInput(paste0("DV_", InfN, "_", ObsN), NULL, 1)),
                                                  column(2, {if (input[[paste0("fix_BSA_", InfN)]]) {
                                                    disabled(textInput(paste0("BSA_", InfN, "_", ObsN), NULL, 1))
                                                  } else {
                                                    textInput(paste0("BSA_", InfN, "_", ObsN), NULL, 1)
                                                  }}),
                                                  column(2, {if (input[[paste0("fix_CLCR_", InfN)]]) {
                                                    disabled(textInput(paste0("CLCR_", InfN, "_", ObsN), NULL, 1))
                                                  } else {
                                                    textInput(paste0("CLCR_", InfN, "_", ObsN), NULL, 1)
                                                  }}),
                                                  column(2, {if (input[[paste0("fix_AGE_", InfN)]]) {
                                                    disabled(textInput(paste0("AGE_", InfN, "_", ObsN), NULL, 1))
                                                  } else {
                                                    textInput(paste0("AGE_", InfN, "_", ObsN), NULL, 1)
                                                  }})
                                                ))
                        )
                      }
                    })
                    lapply(
                      lstToIns_obs,
                      function(ObsN) {
                        if (ObsN!=1) {
                          assign(paste0("react_fix_BSA_", InfN, "_", ObsN), reactive({input[[paste0("fix_BSA_", InfN)]]}), inherits = T)
                          observeEvent(get(paste0("react_fix_BSA_", InfN, "_", ObsN))(), {
                            shinyjs::toggleState(paste0("BSA_", InfN, "_", ObsN), !(get(paste0("react_fix_BSA_", InfN, "_", ObsN))()))
                          })
                          assign(paste0("react_fix_CLCR_", InfN, "_", ObsN), reactive({input[[paste0("fix_CLCR_", InfN)]]}), inherits = T)
                          observeEvent(get(paste0("react_fix_CLCR_", InfN, "_", ObsN))(), {
                            shinyjs::toggleState(paste0("CLCR_", InfN, "_", ObsN), !(get(paste0("react_fix_CLCR_", InfN, "_", ObsN))()))
                          })
                          assign(paste0("react_fix_AGE_", InfN, "_", ObsN), reactive({input[[paste0("fix_AGE_", InfN)]]}), inherits = T)
                          observeEvent(get(paste0("react_fix_AGE_", InfN, "_", ObsN))(), {
                            shinyjs::toggleState(paste0("AGE_", InfN, "_", ObsN), !(get(paste0("react_fix_AGE_", InfN, "_", ObsN))()))
                          })
                        }
                      }
                    )
                }
              })
          })
      }
  })
    
    
    
#Example button
  observeEvent(input$exampleButton, {
      withProgress({
        incProgress(1/7, message = "Inputing example")
        updateTextInput(session, "Ninf", value = 4)
        updateTextInput(session, "nowBSA", value = 1.36)
        updateTextInput(session, "nowCLCR", value = 200)
        updateTextInput(session, "nowAGE", value = 141)
        updateNumericInput(session, "egDummy", value = 1)
      })
  })
  ch_egDummy <- reactive({input$egDummy})
  observeEvent(ch_egDummy(), {
    if (ch_egDummy()==1) {
        withProgress({
          incProgress(2/7, message = "Inputing example")
          updateTextInput(session, "Nobs_1", value = 4)
          updateNumericInput(session, "egDummy", value = 2)
        })
    } else if (ch_egDummy()==2) {
        withProgress({
          incProgress(3/7, message = "Inputing example")
          updateTextInput(session, "Nobs_2", value = 3)
          updateNumericInput(session, "egDummy", value = 3)
        })
    } else if (ch_egDummy()==3) {
        withProgress({
          incProgress(4/7, message = "Inputing example")
          updateTextInput(session, "Nobs_3", value = 3)
          updateNumericInput(session, "egDummy", value = 4)
        })
    } else if (ch_egDummy()==4) {
        withProgress({
          incProgress(5/7, message = "Inputing example")
          updateTextInput(session, "Nobs_4", value = 3)
          updateNumericInput(session, "egDummy", value = 5)
        })
    } else if (ch_egDummy()==5) {
        withProgress({
          incProgress(6/7, message = "Inputing example")
          updateCheckboxInput(session, "fix_BSA_1", value = TRUE)
          updateCheckboxInput(session, "fix_CLCR_1", value = FALSE)
          updateCheckboxInput(session, "fix_AGE_1", value = TRUE)
          updateCheckboxInput(session, "fix_BSA_2", value = TRUE)
          updateCheckboxInput(session, "fix_CLCR_2", value = FALSE)
          updateCheckboxInput(session, "fix_AGE_2", value = TRUE)
          updateCheckboxInput(session, "fix_BSA_3", value = TRUE)
          updateCheckboxInput(session, "fix_CLCR_3", value = FALSE)
          updateCheckboxInput(session, "fix_AGE_3", value = TRUE)
          updateCheckboxInput(session, "fix_BSA_4", value = TRUE)
          updateCheckboxInput(session, "fix_CLCR_4", value = FALSE)
          updateCheckboxInput(session, "fix_AGE_4", value = TRUE)
          updateNumericInput(session, "egDummy", value = 6)
        })
    } else if (ch_egDummy()==6) { 
        withProgress({
          incProgress(7/7, message = "Inputing example")
          updateTextInput(session, "dosePerBSA_1", value = 5)
          updateTextInput(session, "BSA_1", value = 1.36)
          updateTextInput(session, "infTime_1", value = 24)
          
          updateTextInput(session, "TAD_1_2", value = 24)
          updateTextInput(session, "TAD_1_3", value = 48)
          updateTextInput(session, "TAD_1_4", value = 72)
          updateTextInput(session, "DV_1_2", value = 46)
          updateTextInput(session, "DV_1_3", value = 0.3)
          updateTextInput(session, "DV_1_4", value = 0.07)
          updateTextInput(session, "BSA_1_1", value = 1.36)
          updateTextInput(session, "CLCR_1_1", value = 173.6429)
          updateTextInput(session, "CLCR_1_2", value = 191.9211)
          updateTextInput(session, "CLCR_1_3", value = 208.3714)
          updateTextInput(session, "CLCR_1_4", value = 208.3714)
          updateTextInput(session, "AGE_1_1", value = 139)
          
          updateTextInput(session, "dosePerBSA_2", value = 5)
          updateTextInput(session, "BSA_2", value = 1.36)
          updateTextInput(session, "infTime_2", value = 24)
          
          updateTextInput(session, "TAD_2_2", value = 24)
          updateTextInput(session, "TAD_2_3", value = 48)
          updateTextInput(session, "DV_2_2", value = 53)
          updateTextInput(session, "DV_2_3", value = 0.23)
          updateTextInput(session, "BSA_2_1", value = 1.36)
          updateTextInput(session, "CLCR_2_1", value = 169.6047)
          updateTextInput(session, "CLCR_2_2", value = 214.5000)
          updateTextInput(session, "CLCR_2_3", value = 197.1081)
          updateTextInput(session, "AGE_2_1", value = 139)
          
          updateTextInput(session, "dosePerBSA_3", value = 5)
          updateTextInput(session, "BSA_3", value = 1.36)
          updateTextInput(session, "infTime_3", value = 24)
          
          updateTextInput(session, "TAD_3_2", value = 24)
          updateTextInput(session, "TAD_3_3", value = 48)
          updateTextInput(session, "DV_3_2", value = 55)
          updateTextInput(session, "DV_3_3", value = 0.24)
          updateTextInput(session, "BSA_3_1", value = 1.36)
          updateTextInput(session, "CLCR_3_1", value = 173.6429)
          updateTextInput(session, "CLCR_3_2", value = 187.0000)
          updateTextInput(session, "CLCR_3_3", value = 191.9211)
          updateTextInput(session, "AGE_3_1", value = 140)
          
          updateTextInput(session, "dosePerBSA_4", value = 5)
          updateTextInput(session, "BSA_4", value = 1.36)
          updateTextInput(session, "infTime_4", value = 24)
          
          updateTextInput(session, "TAD_4_2", value = 24)
          updateTextInput(session, "TAD_4_3", value = 48)
          updateTextInput(session, "DV_4_2", value = 44)
          updateTextInput(session, "DV_4_3", value = 0.23)
          updateTextInput(session, "BSA_4_1", value = 1.36)
          updateTextInput(session, "CLCR_4_1", value = 155.1702)
          updateTextInput(session, "CLCR_4_2", value = 182.3250)
          updateTextInput(session, "CLCR_4_3", value = 208.3714)
          updateTextInput(session, "AGE_4_1", value = 140)
          
          updateNumericInput(session, "egDummy", value = 0)
        })
    }
  })

#Calculate button
  observeEvent(input$calculate, {
    updateNumericInput(session, "resDummy", value = 0)
    ErrFree <- TRUE
    withProgress({
      tryCatch({
  #Reading input values
          tic("Reading input values")
          incProgress(1/4, message = "Reading input values")
          as.object <- function(name) {
            eval(parse(text = paste0("input$", name)))
          }
        #Cpss
          Cpss_LB <- as.object("def_Cpss_LB") %>% as.numeric
          if (is.numeric(Cpss_LB)&&(Cpss_LB>0)) {} else {stop()}
        #Ninf
          if ((input$Ninf%%1==0)&&(input$Ninf>0)) {Ninf <- input$Ninf} else {stop()}
        #Nobs, dosePerBSA, BSA, infTime
          Nobs <- integer(Ninf)
          dosePerBSA <- numeric(Ninf)
          BSAdose <- numeric(Ninf)
          infTime <- numeric(Ninf)
          TAD <- list()
          DV <- list()
          BSA <- list()
          CLCR <- list()
          AGE <- list()
          for (i in 1:Ninf) {
            Nobs[i] <- as.object(paste0("Nobs_", i))
            dosePerBSA[i] <- as.object(paste0("dosePerBSA_", i)) %>% as.numeric
            BSAdose[i] <- as.object(paste0("BSA_", i)) %>% as.numeric
            infTime[i] <- as.object(paste0("infTime_", i)) %>% as.numeric
            if ((Nobs[i]%%1==0)&&(Nobs[i]>1)) {} else {stop()}
            if (is.numeric(dosePerBSA[i])&&(dosePerBSA[i]>0)) {} else {stop()}
            if (is.numeric(BSAdose[i])&&(BSAdose[i]>0)) {} else {stop()}
            if (is.numeric(infTime[i])&&(infTime[i]>0)) {} else {stop()}
        #TAD, DV, BSA, CLCR, AGE
            curTAD <- numeric(Nobs[i])
            curDV <- numeric(Nobs[i])
            curBSA <- numeric(Nobs[i])
            curCLCR <- numeric(Nobs[i])
            curAGE <- numeric(Nobs[i])
            for (j in 1:Nobs[i]) {
              curTAD[j] <- as.object(paste0("TAD_", i, "_", j)) %>% as.numeric
              curDV[j] <- as.object(paste0("DV_", i, "_", j)) %>% as.numeric
              curBSA[j] <- as.object(paste0("BSA_", i, "_", j)) %>% as.numeric
              curCLCR[j] <- as.object(paste0("CLCR_", i, "_", j)) %>% as.numeric
              curAGE[j] <- as.object(paste0("AGE_", i, "_", j)) %>% as.numeric
              if (as.object(paste0("fix_BSA_", i))) {curBSA[j] <- curBSA[1]}
              if (as.object(paste0("fix_CLCR_", i))) {curCLCR[j] <- curCLCR[1]}
              if (as.object(paste0("fix_AGE_", i))) {curAGE[j] <- curAGE[1]}
              if (is.numeric(curTAD[j])&&(((j==1)&&(curTAD[j]==0))||((j!=1)&&(curTAD[j]>0)))) {} else {stop()}
              if (is.numeric(curDV[j])&&(((j==1)&&(curDV[j]==0))||((j!=1)&&(curDV[j]>0)))) {} else {stop()}
              if (is.numeric(curBSA[j])&&(curBSA[j]>0)) {} else {stop()}
              if (is.numeric(curCLCR[j])&&(curCLCR[j]>0)) {} else {stop()}
              if (is.numeric(curAGE[j])&&(curAGE[j]>0)) {} else {stop()}
            }
            if (identical(curTAD, sort(curTAD))&&((curTAD %>% unique %>% length)==(curTAD %>% length))) {} else {stop()}
            TAD <- c(TAD, list(curTAD))
            DV <- c(DV, list(curDV))
            BSA <- c(BSA, list(curBSA))
            CLCR <- c(CLCR, list(curCLCR))
            AGE <- c(AGE, list(curAGE))
          }
        #BSA, CLCR, AGE (for this dose)
          nowBSA <- as.object("nowBSA") %>% as.numeric
          nowCLCR <- as.object("nowCLCR") %>% as.numeric
          nowAGE <- as.object("nowAGE") %>% as.numeric
          if (is.numeric(nowBSA)&&(nowBSA>0)) {} else {stop()}
          if (is.numeric(nowCLCR)&&(nowCLCR>0)) {} else {stop()}
          if (is.numeric(nowAGE)&&(nowAGE>0)) {} else {stop()}
        #ending message
          output$msgCal <- renderText("Inputs values okay!")
    toc()
  #Getting individual parameters
    tic("Estimating individual parameters")
    incProgress(1/4, message = "Estimating individual parameters")
    resOptim <- optim(par = c(0,0,0,0), fn = EBE_optim, method = "L-BFGS-B",
                      DATA = list(
                        Ninf = Ninf, TAD = TAD, infTime = infTime,
                        dosePerBSA = dosePerBSA, BSAdose = BSAdose,
                        BSA = BSA, CLCR = CLCR, AGE = AGE,
                        DV = DV
                      ))
    assign("inspect_Optim", resOptim, envir = .GlobalEnv)
    #Get etas
      ETA_CL <- resOptim$par[1]
      ETA_V1 <- resOptim$par[2]
      ETA_Q <- resOptim$par[3]
      ETA_V2 <- resOptim$par[4]
      ICL_IND <- TVCL*(nowBSA/0.735)**B_BSA_CL*(nowCLCR/192)**B_CLCR_CL*exp(ETA_CL)
      IV1 <- TVV1*(nowBSA/0.735)**B_BSA_V1*exp(ETA_V1)
      IQ <- TVQ*(nowAGE/63.5)**B_AGE_Q*exp(ETA_Q)
      IV2 <- TVV2*exp(ETA_V2)
    toc()
  #Optimizing dose
    tic("Optimizing dose")
    incProgress(1/4, message = "Optimizing dose")
      #Defining Cp function
        CpFcn <- function(t) {
          if (t<=24) {
            (dose*2200.51*nowBSA/24*(IA/Ialpha*(1-exp(-Ialpha*t))+IB/Ibeta*(1-exp(-Ibeta*t)))) %>% return()
          } else {
            (dose*2200.51*nowBSA/24*(IA/Ialpha*(1-exp(-Ialpha*24))*exp(-Ialpha*(t-24))+IB/Ibeta*(1-exp(-Ibeta*24))*exp(-Ibeta*(t-24)))) %>% return()
          }
        }
      #Dose and Cpss
        #Calculating the 95% confidence dose
          final_P95_dose <- (Cpss_LB/exp(qnorm(0.05,0,1)*sqrt(OMIOV+SG))*24*ICL_IND/2200.51/nowBSA) %>% as.vector
      #Getting expected PK profile at optimal dose (by simulation)
        timePt <- seq(0.1, 96, by = 0.1)
        #Predicted expected Cp
          ICL <- ICL_IND
          Ik12 <- IQ/IV1
          Ik21 <- IQ/IV2
          Ik10 <- ICL/IV1
          Ibeta <- (Ik12+Ik21+Ik10-sqrt((Ik12+Ik21+Ik10)**2-4*Ik21*Ik10))/2
          Ialpha <- Ik21*Ik10/Ibeta
          IA <- (Ialpha-Ik21)/IV1/(Ialpha-Ibeta)
          IB <- (Ibeta-Ik21)/IV1/(Ibeta-Ialpha)
          dose <- final_P95_dose
          ExpectedConc <- sapply(timePt, CpFcn)
          PredProf <- data.frame(timePt, ExpectedConc)
          sampleSize <- 5000
          ETA_IOV <- rnorm(sampleSize, 0, sqrt(OMIOV))
          EPS <- rnorm(sampleSize*10, 0, sqrt(SG))
          PredConc <- sapply(ETA_IOV, function(eta_IOV) {
            ICL <- ICL_IND*exp(eta_IOV)
            Ik12 <- IQ/IV1
            Ik21 <- IQ/IV2
            Ik10 <- ICL/IV1
            Ibeta <- (Ik12+Ik21+Ik10-sqrt((Ik12+Ik21+Ik10)**2-4*Ik21*Ik10))/2
            Ialpha <- Ik21*Ik10/Ibeta
            IA <- (Ialpha-Ik21)/IV1/(Ialpha-Ibeta)
            IB <- (Ibeta-Ik21)/IV1/(Ibeta-Ialpha)
            dose <- final_P95_dose
            Conc_24h <- CpFcn(24)*exp(sample(EPS, 1))
            Conc_48h <- CpFcn(48)*exp(sample(EPS, 1))
            Conc_72h <- CpFcn(72)*exp(sample(EPS, 1))
            return(c(Conc_24h, Conc_48h, Conc_72h))
          }) %>% t %>% as.array %>% c(rep(24,sampleSize), rep(48,sampleSize), rep(72, sampleSize), .) %>% matrix(., ncol = 2) %>% as.data.frame
          colnames(PredConc) <- c("timePt", "Conc")
      #Generating graphs
        #Cpss-dose graph
          dose_LB <- 0.1
          dose_UB <- ceiling(final_P95_dose + 3)
          dose_seq <- seq(dose_LB, dose_UB, 0.1)
          lst_Cpss_P025 <- exp(qnorm(0.025, log(dose_seq*2200.51*nowBSA/24/ICL_IND), sqrt(OMIOV+SG)))
          min025 <- min(lst_Cpss_P025)
          lst_Cpss_P050 <- exp(qnorm(0.050, log(dose_seq*2200.51*nowBSA/24/ICL_IND), sqrt(OMIOV+SG)))
          lst_Cpss_P500 <- exp(qnorm(0.500, log(dose_seq*2200.51*nowBSA/24/ICL_IND), sqrt(OMIOV+SG)))
          lst_Cpss_P950 <- exp(qnorm(0.950, log(dose_seq*2200.51*nowBSA/24/ICL_IND), sqrt(OMIOV+SG)))
          lst_Cpss_P975 <- exp(qnorm(0.975, log(dose_seq*2200.51*nowBSA/24/ICL_IND), sqrt(OMIOV+SG)))
          max975 <- max(lst_Cpss_P975)
          assign("test_df", data.frame(dose_seq, lst_Cpss_P025, lst_Cpss_P050, lst_Cpss_P500, lst_Cpss_P950, lst_Cpss_P975), envir = .GlobalEnv)
          plot_1 <- ggplot(data.frame(dose_seq, lst_Cpss_P025, lst_Cpss_P050, lst_Cpss_P500, lst_Cpss_P950, lst_Cpss_P975), aes(x=dose_seq, y=lst_Cpss_P500)) +
                      geom_line(aes(y=lst_Cpss_P500), linetype = "solid", size = 1) +
                      geom_line(aes(y=lst_Cpss_P050), linetype = "dashed", size = 1) +
                      geom_ribbon(aes(ymin=lst_Cpss_P025, ymax=lst_Cpss_P975), fill="blue", alpha=0.15) +
                      geom_ribbon(aes(ymin=lst_Cpss_P050, ymax=lst_Cpss_P950), fill="green", alpha=0.15) +
                      geom_hline(yintercept=Cpss_LB, linetype = "dashed", size = 1) +
                      geom_vline(xintercept=final_P95_dose, linetype = "dotted", size = 1) +
                      scale_x_continuous(name = "Dose (g/m2)",
                                         limits = c(dose_LB, dose_UB),
                                         breaks = c(round(dose_LB,1), round(final_P95_dose,1), round(dose_UB,1))) +
                      scale_y_continuous(name = "Steady-state plasma level (uM)",
                                         limits = c(min025, max975),
                                         breaks = c(1,10,Cpss_LB,100,1000),
                                         trans = "log10") +
                      ggtitle("Predicted steady-state level against dose per BSA") +
                      theme(axis.text=element_text(size=18),
                            axis.title=element_text(size=21,face="bold"),
                            title = element_text(size = 22, vjust = 0.5))
          plot_2 <- ggplot(PredProf, aes(x=timePt, y=ExpectedConc)) +
            stat_boxplot(data = PredConc, aes(group = factor(timePt), y = Conc), geom = "errorbar", width = 10) +
            geom_boxplot(data = PredConc, aes(group = factor(timePt), y = Conc), width = 10) +
            geom_line(size = 1) +  
            scale_x_continuous(name = "Time after dose (hours)", breaks = c(24,48,72,96)) +
            scale_y_continuous(name = "Predicted plasma concentration (uM/L)", trans = "log10", breaks = c(0.1,1,10,100)) + 
            ggtitle("Predicted plasma concentration profile") +
            theme(axis.text=element_text(size=12),
                  axis.title=element_text(size=15,face="bold"),
                  title = element_text(size=16, vjust = 0.5))
        
      toc()
  #Generating output
    tic("Generating output")
    incProgress(1/4, message = "Generating output")
      output$optDPA <- renderText(paste0("Minimum dose per BSA (rounded to 1 d.p.) = ", round(final_P95_dose,1), " g/m2"))
      output$optDose <- renderText(paste0("Optimal dose (rounded to 2 d.p.) = ", round(final_P95_dose*nowBSA,2), " g"))
      output$plot_1 <- renderPlot({plot_1})
      output$plot_2 <- renderPlot({plot_2})
      output$Legend_1 <- renderText("Blue region shows the 95% confidence interval. Green region shows the 90% confidence interval")
    toc()
    output$msgCal <- renderText("Done!")
      }, warning = function(w) {
      }, error = function(e) {
        output$msgCal <- renderText("Problems with input values! Please check!")
      }, finally = {
      })
    })
    updateNumericInput(session, "resDummy", value = 1)
  })
}

shinyApp(MTX_UI, MTX_Server)