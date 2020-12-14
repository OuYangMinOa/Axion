## single_bin.jl 
    Simulation the signal out from the  cavity and after the  Amplifier and the FFT
## How to use ?
    Download Julia in https://julialang.org/downloads/
open julia and type
```
] add Distributed PyPlot CurveFit CSV DataFrames Distribution
include("single_bin.jl")
```
and you can start the code
### functions 
* single_bin( plotting = true , number=1 , SNR=2)
    ```
    Simulate the single bin signal
    ```
*  six_bin(plotting = true , number=1 , SNR=2 , v = 0.00375 , M_num=6)
    ```
    Simulate the wide bin signal
    ```
* rescan_six_bin_no_plot(change=0, SNR=1.645, change2=-1, v=0.0041, M_num=4)
    ```
    Simulate the wide bin signal and rescan onetime without potting
    ```
* rescan_six_bin_plot(change=0, SNR=1.645, change2=-1, v=0.0041, M_num=4)
    ```
    Simulate the wide bin signal and rescan onetime with plotting
    ```

