function MSE = get_mse(in_ref,in_test)

    MSE = mean(abs((in_ref-in_test)*(in_ref-in_test)'));