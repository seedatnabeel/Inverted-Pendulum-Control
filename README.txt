ELEN 4016: CONTROL 2 PROJECT

SIMULINK AND MATLAB SCRIPT README

Nabeel Seedat (719484)
Irfaan Mohamed (547347)

Compile and Run the .slx models on MATLAB 2015A with Symbolic and Fuzzy Logic Toolbox

Make sure the MATLAB working directory is in the same folder as the Simulink and .m files
i.e. the MATLAB working directory must be the folder on the USB.


Simulink files
 
N.B. The .m files don't need to be run: they automatically run when the Simulink model is run
 - Provided the working directory is correct

(1) single_pend_linear.slx : Simulates the control via state feedback of the linear single inverted pendulum on a cart

Associated .m file: single_pend.m

(2) single_pend_non_linear.slx: Simulates the control via state feedback of the non-linear single inverted pendulum on a cart

Associated .m file:single_pend.m

(3) fuzzy_controller.slx: Simulates the Fuzzy Logic Controller for the Double Inverted pendulum on a cart

Associated .m file: fuzzy_logic.m

(4) double_feedback.slx: Simulates the state feedback (pole placement and LQR) control of a Double Inverted pendulum on a cart

Associated .m file: double_pend_feedback.m

(5) double_luenberger.slx: Simulates the Luenberger observer with state feedback control of a Double Inverted pendulum on a cart
                           It also shows graph of mic measurement vs real, the scaling factor and measurement error.

Associated .m file: double_pend_feedback.m 