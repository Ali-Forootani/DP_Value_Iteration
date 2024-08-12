# Dynamic-Programming-for-MDP
Sure! Here's a README file explaining the MATLAB script for revenue optimization using Dynamic Programming:

# Revenue Optimization using Dynamic Programming

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

## Script Explanation
The script consists of the following main sections:

1. **Input**: The user is prompted to enter the number of prices (m) and the number of resources (N). Additional inputs may be required, such as loading matrices containing price vectors, time horizons, and state spaces.

2. **State Space Construction**: The script constructs a state space matrix (S) containing all possible combinations of at most N resources on m places. The base conversion routine is used to generate combinations, and any combinations exceeding the resource limit are excluded.

3. **Terminal Period**: The value of each state at the terminal time step is calculated and saved in the `Terminal_revenue` variable. This represents the revenue at the end of the time horizon.

4. **Dynamic Programming**: The script applies Dynamic Programming to evaluate the optimal decision-making strategy and revenue for each state and time step. It iterates over the time steps in reverse order, evaluating each state and action based on transition probabilities and previously calculated values.

5. **Frequency Distribution Before Reservation**: The script calculates the frequency distribution of decisions (actions) before reservation for each time step. The results are saved in the `frequency_decision_BD` variable.

6. **Plotting**: The script includes plotting commands to visualize the frequency distribution of decisions over time. The results can be displayed in multiple subplots.

## Additional Notes
- Make sure to provide all required inputs, such as loading matrices containing price vectors, time horizons, and state spaces, before running the script.
- The script might require modifications or additional code to handle specific problem instances or constraints.
- If you encounter any issues or have questions, please feel free to contact the author.

## Article

1. Forootani, Ali, et al. "Modelling and solving resource allocation problems via a dynamic programming approach." International Journal of Control 94.6 (2021): 1544-1555.

you can contact me via: aliforootani@ieee.org

## License
This script is provided under the [MIT License](https://opensource.org/licenses/MIT).
