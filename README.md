# Basic code to solve differential equations in [Herd immunity thresholds for SARS-CoV-2 estimated from unfolding epidemics](https://www.medrxiv.org/content/10.1101/2020.07.23.20160762v2)

Core models of COVID19 accounting for individual variation in susceptibility or exposure to infection implemented in MATLAB.

These models were analysed and used to analyse the #SARSCoV2 pandemic in our 2 medRxiv preprints. 

First initially posted in April:  
[M. Gabriela M. Gomes, Rodrigo M. Corder, Jessica G. King, Kate E. Langwig, Caetano Souto-Maior, Jorge Carneiro, Guilherme Goncalves, Carlos Penha-Goncalves, Marcelo U. Ferreira, Ricardo Aguas, Individual variation in susceptibility or exposure to SARS-CoV-2 lowers the herd immunity threshold](https://www.medrxiv.org/content/10.1101/2020.04.27.20081893v3)

Second in July:  
[R. Aguas, R. M. Corder, J. G. King, G. Gonçalves, M. U. Ferreira, M. G. M. Gomes, Herd
immunity thresholds for SARS-CoV-2 estimated from unfolding epidemics](https://www.medrxiv.org/content/10.1101/2020.07.23.20160762v2)

Neither has completed peer review yet. We will continue to update and add code to this repository and post preprints with new analyses. 

Happy to take comments and suggestions. Happy to collaborate. Keep well!

## How to generate plots in Figures 1, 2, and Extended Data Figure 5:

Run `Epidemic.m` and answer questions as follows:

```
*** Choose a model: ***
1. Variable susceptibility
2. Variable connectivity
3. Variable connectivity (reducing CV during social distancing)
```

_Enter:_
- _1 for Figure 1;_ 
- _2 for Figure 2;_ 
- _3 for Extended Data Figure 5._

```
*** Insert parameters: ***
Initial R0 =
```

_Enter:_ <img src="https://render.githubusercontent.com/render/math?math=R_0"> _displayed in plot you wish to reproduce (also in Extended Data Tables 1, 2, 3, for
models 1, 2, 3, respectively)._

```
Coefficient of variation in susceptibility (> or = 0): CV =
```

_Enter: 𝐶𝑉 displayed in plot you wish to reproduce (also in Extended Data Tables 1, 2, 3, for
models 1, 2, 3, respectively)._

```
Social distancing (0 - 1) =
```

_Enter:_ <img src="https://render.githubusercontent.com/render/math?math=d_{max}"> _displayed in plot you wish to reproduce (also in Extended Data Tables 1, 2, 3, for
models 1, 2, 3, respectively)._

```
Time (in days) to initial social distancing measures =
```

_Enter:_ <img src="https://render.githubusercontent.com/render/math?math=round(t_0^d)-t_0">, _where_ <img src="https://render.githubusercontent.com/render/math?math=t_0^d"> _is as displayed in Extended Data Tables 1, 2, 3, for models
1, 2, 3, respectively, and_ <img src="https://render.githubusercontent.com/render/math?math=t_0"> _is model-independent: Belgium (1 day); England (29 days);
Portugal (3 days); Spain (8 days)._

Note: Resulting plots will be approximations due to inevitable rounding errors in reported
parameter estimates.
