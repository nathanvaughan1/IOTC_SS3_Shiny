This is a Dockerfile repository to build a shinyproxy image for the BFT-E shiny app
the BFT-E app is hosted at 
https://github.com/TahaImz/shinyapp_BFT-E_retrosTable

The docker image hosted at 
https://hub.docker.com/r/nathanvaughan/shinyproxy-bft-e/

it is automaticaly built from the shiny app repository hosted at  
https://github.com/nathanvaughan1/shinyproxy-BFT-E

and will be updated if this repository is edited.

Additional shiny proxy information
    
name: BFT-E

display-name: BFT-E shiny app

docker-cmd: ["R", "-e shiny::runApp('/root/BFT')"]

docker-image: nathanvaughan/shinyproxy-BFT-E
