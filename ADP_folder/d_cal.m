function [d_k summation_dk] = d_cal(ph_k,k,summation,current_state,g_state)

summation_dk=summation+(ph_k*g_state);

d_k=(1/(k+1))*(summation_dk);