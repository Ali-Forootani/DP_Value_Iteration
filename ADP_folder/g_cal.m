function [g_state] = g_cal(current_state,price_vector)
g_state=dot(price_vector,current_state);