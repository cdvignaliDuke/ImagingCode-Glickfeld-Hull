function Ca_ker = Calcium_kernel (tau, frame_length)

tau_in_frame = tau / frame_length;
Ca_ker_size = ceil(tau_in_frame*3);
Ca_ker=exp(-[1:Ca_ker_size]/tau_in_frame);
Ca_ker=Ca_ker/sum(Ca_ker);
