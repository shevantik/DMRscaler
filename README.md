# DMRscaler
Readme under construction:
To install DMRscaler, first download the package repository from https://github.com/leroybondhus/DMRscaler
Then make sure you have "devtools" installed, this can be done with the command:

>  if(!require("devtools", quietly = TRUE )){
>   install.packages("devtools")
>   library("devtools")
>  }

After this, from R run the command:

> devtools::install(path_to_DMRscaler_directory)

replacing the variable _path_to_DMRscaler_directory_ with the path to the downloaded DMRscaler directory. 
If running R from the downloaed DMRscaler directory, this can be done by setting _path_to_DMRscaler_directory = "../DMRscaler" 


To get set up quickly with DMRscaler, follow the instructions in the DMRscaler_User_Guide_simplified.html,
a more detailed user manual that includes some real data analysis that highlights the utility of DMRscaler's functionality is
included in the User manual available at https://leroybondhus.github.io/DMRscaler/
