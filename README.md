# Hawkes process modeling of COVID-19 with mobility leading indicators and spatial covariates
## Overall Framework on Hawkess processes component
<img src="./figure/Frame_work.png" width="800px">

## Dependencies 
Matlab 2018b

## Quick run
### Use MLE estimate shape and scale parameters for Weibull distribution
#### HawkPR('./input_data/NYT_Dconfirmed.csv', './input_data/GoogleMobi_Dconfirmed.csv', './input_data/Demo_Dconfirmed.csv', 14, '', '', 200, 7, 100, './output/mdl.mat', './output/pred.csv')
### Specify shape and scale parameters for Weibull distribution
#### HawkPR('./input_data/NYT_Dconfirmed.csv', './input_data/GoogleMobi_Dconfirmed.csv', './input_data/Demo_Dconfirmed.csv', 14, 8, 4, 200, 7, 100, './output/mdl.mat', './output/pred.csv')


