This is a Dockerfile repository to build a shinyproxy image for the IOTC_SS3 shiny app
the IOTC_SS3 app is hosted at 
https://github.com/aenieblas/IOTC_SS3_Shiny

The docker image hosted at 
https://hub.docker.com/r/nathanvaughan/shinyproxy-IOTC_SS3/

it is automaticaly built from the shiny app repository hosted at  
https://github.com/nathanvaughan1/IOTC_SS3_Shiny

and will be updated if this repository is edited.

Additional shiny proxy information
    
name: IOTC_SS3

display-name: IOTC_SS3 shiny app

docker-cmd: ["R", "-e rmarkdown::render('/root/IOTC_SS3/ss3_dashboard_final.rmd')"]

docker-image: nathanvaughan/shinyproxy-IOTC_SS3
