
# Revenue Optimization using Dynamic Programming/ ADP.(On the Value Iteration algorithm of Dynamic Programming) 

This MATLAB script solves a revenue optimization problem using Dynamic Programming. It calculates the optimal revenue and decision-making strategy over a given time horizon, considering the constraints on prices and resources.

## Requirements
- MATLAB (R2018a or above)

## Usage
1. Open the MATLAB software.
2. Set the current directory to the location where the script is saved.
3. Run the script.
4. Follow the instructions prompted by the script:
   - Enter the number of prices (m).
   - Enter the number of resources (N).
   - Provide any additional inputs if required (such as loading matrices).

## Main .m files

**Exact DP**: To run the DP algorithm for the MDP you should run `Exact DP.m`!

**Approximate DP**: To run the ADP algorithm for the MDP you shoudl run `ADP_folder/LSTD_max_norm_version_3.m`

**Exact DP  Periodic**: To run th periodical DP algorithm for the MDP with Sinusoidal, square wave and sawtooth wave you should run `Exact_DP_Periodic.m`!

- sinusoidal transition probabilities:

    - $\lambda=[0.5+0.1*cos(k*pi/10);0.3+.1*sin(k*pi/10);0.2+0.1*cos(k*pi/10); 0.15-0.05*cos(k*pi/10)]$;

    - $mu=[0.3-0.1*cos(k*pi/10);0.4-0.1*sin(k*pi/10);0.5-0.1*cos(k*pi/10); 0.6+0.05*cos(k*pi/10)]$;
    


## Script Explanation
The script consists of the following main sections:

1. **Input**: The user is prompted to enter the number of prices (m) and the number of resources (N). Additional inputs may be required, such as loading matrices containing price vectors, time horizons, and state spaces.

2. **State Space Construction**: The script constructs a state space matrix (S) containing all possible combinations of at most N resources on m places. The base conversion routine is used to generate combinations, and any combinations exceeding the resource limit are excluded.

3. **Terminal Period**: The value of each state at the terminal time step is calculated and saved in the `Terminal_revenue` variable. This represents the revenue at the end of the time horizon.

4. **Dynamic Programming**: The script applies Dynamic Programming to evaluate the optimal decision-making strategy and revenue for each state and time step. It iterates over the time steps in reverse order, evaluating each state and action based on transition probabilities and previously calculated values.

5. **Frequency Distribution of Different Decisions**: The script calculates the frequency distribution of decisions (actions) before reservation for each time step. The results are saved in the `frequency_decision_BD` variable.

6. **Plotting**: The script includes plotting commands to visualize the frequency distribution of decisions over time. The results can be displayed in multiple subplots.

## Additional Notes
- Make sure to provide all required inputs, such as loading matrices containing price vectors, time horizons, and state spaces, before running the script.
- The script might require modifications or additional code to handle specific problem instances or constraints.
- If you encounter any issues or have questions, please feel free to contact the author.

### Examples:

1. **Example 1**: Applying DP algorithm for a resource allocation problem with `m=4`, and `N=4`, with $\alpha=0.9$, $c= [0.8,\ 0.9,\ 1,\ 1.1]$. For this configuration the set of birth probabilities is $\lambda= [0.090596,\ 0.048632,\ 0.015657,\ 0.005088],$ and the death probabilities is 
$\mu = [0.483723,\ 0.444019,\ 0.024843,\ 0.335103]$.

2. **Example 2**: Applying random policy on the MDP in example 1

3. **Example 3**: Applying off-policy LSTD approach on the MDP for large scale resource allocation problem  `LSTD_max_norm_version_3.m`

- We might use $\eta$ or $\sigma$ in place of $\lambda$ and $\mu$!

4. **Example 7**: Sinusoidal perturbation on MDP.

5. **Example 8**: Square wave perturbation on MDP.

6. **Example 9**: Sawtooth wave perturbation on MDP.

7. **Appendix/floquet_vs_nofloquet.m**: Floquet Transvormation.

### ADP Folder 

This directory includes all the necessary files to run the LSTD approach. The file `LSTD_max_norm_version_3.m` should be run and you required to insert `m` and `N` as the input parameters! 
- Note that you need to specify $\lambda$ and $\mu$, here we use $\lambda= [0.090596,\ 0.048632,\ 0.015657,\ 0.005088],$ and the death probabilities is 
$\mu = [0.483723,\ 0.444019,\ 0.024843,\ 0.335103]$.


## Articles

1. Forootani, Ali, et al. "Modelling and solving resource allocation problems via a dynamic programming approach." International Journal of Control 94.6 (2021): 1544-1555.
2. Forootani, Ali, et al. "Approximate dynamic programming for stochastic resource allocation problems." IEEE/CAA Journal of Automatica Sinica 7.4 (2020): 975-990.
3. Forootani, Ali, et al. "A least-squares temporal difference based method for solving resource allocation problems." IFAC Journal of Systems and Control 13 (2020): 100106.
4. Forootani, Ali, et al. "Allocating resources via price management systems: a dynamic programming-based approach." International Journal of Control 94.8 (2021): 2123-2143.
5. Forootani, Ali, et al. "Off-Policy Temporal Difference Learning for Perturbed Markov Decision Processes." IEEE Control Systems Letters, 8, (2024): 3488-3493.

you can contact me via: aliforootani@ieee.org/aliforootani@gmail.com

## License
This script is provided under the [MIT License](https://opensource.org/licenses/MIT).
