This is a Dockerfile repository to build a shinyproxy image for the IOTC_SS3 shiny app
the IOTC_SS3 app is hosted at 
https://github.com/aenieblas/IOTC_SS3_Shiny

The docker image hosted at 
https://hub.docker.com/r/nathanvaughan/shinyproxy-iotc_ss3/

it is automaticaly built from the shiny app repository hosted at  
https://github.com/nathanvaughan1/IOTC_SS3_Shiny

and will be updated if this repository is edited.

Additional shiny proxy information
    
name: SS3 Diagnostic Plots

display-name: SS3 Diagnostic Plots

docker-cmd: ["R", "-e rmarkdown::run('/root/IOTC_SS3/ss3_dashboard_final.Rmd')"]

docker-image: nathanvaughan/shinyproxy-iotc_ss3
