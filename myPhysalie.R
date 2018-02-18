#####################################
### Script Shiny myPhysalie ###
#####################################

# Master 2 Bioinformatique
# Par Maira BARCA, Mauriac RAZAFIMAFONTY et Rokhaya BA

#--
#indiquer le chemin vers le repertoire courant
#setwd ("") 

#lancement de l'application shiny en cliquant sur le bouton "Run App" 

#--
#instalation des packages 
#une fois les packages installés une 1ere fois (ou si deja installes sur l'appareil), la commande suivante peut etre mise en commentaire
#install.packages (c ("shiny", "ggmap", "maptools", "maps", "ggplot2", "FactoMineR", "ade4", "PMCMR"))

#chargement de toutes les librairies
library(shiny)
library("ggmap")
library(maptools)
library(maps)
library(ggplot2) 
library(FactoMineR)
library(ade4)
library("PMCMR")

#chargement des fichiers
coord<-read.csv("base1_myPhysalie.csv",sep = ';',dec=',',header=TRUE)
coord <- na.omit(coord)

data<-read.csv("base2_myPhysalie.csv", sep = ';',dec=',',header=TRUE)
data = data[,-1]
colnames(data)[9]="RTP"
colnames(data)[63]="ratio1"
colnames(data)[64]="ratio2"
colnames(data)[65]="ratio3"
colnames(data)[66]="ratio4"
colnames(data)[67]="ratio5"
colnames(data)[68]="ratio6"
colnames(data)[69]="ratio7"
colnames(data)[70]="ratio8"
data <-na.omit(data)

#separation en plusieurs sous tables
geo<-data[,c(1,5:8)]
reseau<-data[,c(1,9:13)]
ratio<-data[,c(1,63:73)]
activite<-data[,c(74:100)]
rtp<-data[,c(10:13)]
ratios<-ratio[,-1]

#definition des variables de l'onglet 1
#systeme du reseau
sys <- c("BayofBiscay", "Mediterranee_Baie_Bizerte", "Marais_desseches_Charente_Maritime", "Polynesie_Atoll_Takapoto", "Bassin_Arcachon_interne", 
         "Sargasse_center_gyres", "Sargasse_edge_gyres", "Sargasse_Bermuda", "Arabian_Sea",
         "Pacifique_equatorial", "Mer_Floride", "Bassin_Marennes_Oleron")

#definition des variables de l'onglet 3
#selection des variables d'interet du script d'emma pour l'analyse factorielle multiple
orgstab<-data[,c(75,76,78,89,90,91,93,94,99,100)]
activitetab<-data[,c(79,82,84,85,86,103)]
#concatenation des 3 matrices de donnees
tab<-data.frame(rtp,activitetab,orgstab); names(tab) 
gr<-c(ncol(rtp),ncol(activitetab),ncol(orgstab))
#realisation de l'analyse multifactorielle
mfa<-MFA(tab, group=gr, type=c("c","s","c"), ncp=25, name.group= c("RTP","Act", "Org"))
#variables a utiliser dans l'UI
choices <- c("Visualisation des variables", "Visualisation des stations")
stations <- c("Dim1vsDim2", "Dim1vsDim3","Dim2vsDim3")

# Definition de l'interface UI
ui <- fluidPage( # Application title
  titlePanel("ShinyApp."),
  tabsetPanel(
  #onglet1
  
    tabPanel("Localisation des systèmes dans l'espace",
             sidebarLayout(
               sidebarPanel(
                 h3("Localisation des systemes en fonction de leurs coordonnees geographiques", align="center"),
                 selectInput("syst", "Choix du système", choices = sys, selected = "BayofBiscay")), 
               
               mainPanel(
                 plotOutput(outputId="distPlot", height = "800px", width = "1000px")))
    ),
  
  #onglet2
    tabPanel("Analyse de variance des rations en fonction des RTP",
             titlePanel("Deuxième onglet"),
             sidebarLayout(
               sidebarPanel(
                 selectInput("ratioss", "Ratios:", choices = colnames(ratios)),
                 checkboxInput(inputId = "postHoc", label = strong("Post-hoc de Nemenyi"), value = FALSE)),
               
               mainPanel(
                 plotOutput(outputId="boxPlot"),
                 verbatimTextOutput("printKW"),
                 verbatimTextOutput("postHoc")))
             ),
  
  #onglet3
    tabPanel("Analyse multivariée (AFM) des variables d'interet",
             titlePanel("Troisieme onglet"),
             
             sidebarPanel(
               selectInput("type", "Type de graphique:", choices = choices, selected = "Visualisation des variables"),
               selectInput("station", "Choix des stations:", choices = stations, selected = "Dim1vsDim2"),
               p(),
               h4("Informations sur les variables :"),
               p("AoDC : Ascendance relative"),
               p("AMI : Information mutuelle moyenne"),
               p("CCI : Indice du cycle de Fynn"),
               p("DH : Détrivorie/Herbivorie"),
               p("HR : Entropie, diversité de flux"),
               p("APL : Longueur moyenne du trajet d'une unité de C"),
               p("SOI : Indice d'omnivorie moyen"),
               p("RoDC : Redondance relative"),
               p("PB : Production/Biomasse"),
               p("TST : Somme des flux entrants ou sommes des flux sortants"),
               p("totRES : Production totale"),
               p("totNPP : Production nette totale"),
               p("ProdoRes : Productivité du réseau"),
               p("T.. : Somme totale des flux du systeme")),
             
             mainPanel(
               plotOutput(outputId="drawCircle",height = "850px", width = "850px"))
    )
  
  ))


#--
# Definition du serveur
server <- function(input, output) {
  
  #fonction pour l'onglet1
  #creation d'une carte du monde qui va etre rempli avec les coordonnees du reseau
  output$distPlot <- renderPlot({
    mp <- NULL
    mapWorld <- borders("world", colour="gray50", fill="gray50") # create a layer of borders
    mp <- ggplot() +   mapWorld
    
    #en fonction du systeme choisie, on fait une sous table avec l'ensemble des sites de ce systeme
    #calcul de la moyenne de la latitude et longitude pour n'avoir qu'une valeur correspondante au systeme demande
    df1<-subset(coord, Système==input$syst)
    lat = mean(df1$Latitude)
    long = mean(df1$Longitude)
    
    #toutes les coordonnees du reseau sont inserees et seront visibles en bleu
    mp <- mp+ (geom_point(aes(x=coord$Longitude, y=coord$Latitude) ,color="blue", size=3))
    #les coordonnees choisie par l'utilisateur vont etre visible en rouge
    mp <- mp+ (geom_point(aes(x=long, y=lat) ,color="red", size=3))  
    mp
  })
  
  #--
  #fonctions pour l'onglet2
  output$boxPlot <- renderPlot({
    #creation des boxplot des ratios en fonction des RTP
    pl <- ggplot(data = data,
                 aes_string(fill="RTP", x="RTP",
                            y=input$ratioss))
    # draw the boxplot for the specified trait
    pl + geom_boxplot() 
  })
  
  output$printKW <- renderPrint({
    #test de Kruskal-wallis en fonction du ration demande par l'utilisateur
    if(kruskal.test(get(input$ratioss)~RTP, data)$p.value<0.05){paste("Test de Kruskall Wallis significatif")}
      })
  
  output$postHoc <- renderPrint({
    #test post Hoc de nemenyi en fonction du ration demande par l'utilisateur
    if (input$postHoc) {
      post <- posthoc.kruskal.nemenyi.test(data=data, get(input$ratioss)~RTP) 
      summary(post)}
  })
  
  #-- 
  #fonction pour l'onglet 3
  output$drawCircle <- renderPlot({
    #methode pour visualisation des resultats de l'afm
    #graphe va etre modifie en fonction des choix de l'utilisateur
    if(input$type=="Visualisation des variables"){
      if(input$station=="Dim1vsDim2"){
        plot(mfa,choix="var",axes = c(1,2),cex=1)}
      else if(input$station=="Dim1vsDim3"){
        plot(mfa,choix="var",axes = c(1,3),cex=1)}
      else if(input$station=="Dim2vsDim3"){
        plot(mfa,choix="var",axes = c(2,3),cex=1) }
    }
    
    if(input$type=="Visualisation des stations"){
      if(input$station=="Dim1vsDim2"){
        s.class(mfa$ind$coord[,c(1,2)] , fac=data$Systeme)}
      else if(input$station=="Dim1vsDim3"){
        s.class(mfa$ind$coord[,c(1,3)] , fac=data$Systeme)}
      else if(input$station=="Dim2vsDim3"){
        s.class(mfa$ind$coord[,c(2,3)] , fac=data$Systeme) }
    }
  })
  
}


#--
# Creation du Shiny
shinyApp(ui = ui, server = server)