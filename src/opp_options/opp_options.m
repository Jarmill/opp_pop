classdef opp_options < handle
    %OPP_OPTIONS Options data structure for the optimal pulse pattern task
    %   Detailed explanation goes here

    properties
        %device parameters
        f0(1, 1)double{mustBeNonnegative} = 50; %frequency in Hz
        Ts(1, 1)double{mustBeNonnegative} = 1e-4; %minimal inter-switch sample time (sec)
        Z_load(1, 1)double = 0+0j; %impedance of the load (given frequency f0)
        L = [-1, 0, 1]; %levels of the inverter
        

        %harmonics constraints
        harmonics = [];
        harmonics_load = [];


        %TODO:
        %power budget (and allocation of switching device parameters)
        Power_max(1, 1)double = 1000; 
        %power loss model
        %grid-side filter definitions (includes DC tank capacitor with RMS current constraint)

        %switching parameters
        k(1, 1)int32 = 4; %number of switches in the sequence
        Symmetry(1, 1)opp_symmetry = opp_symmetry.Full %symmetry type
        three_phase(1, 1)opp_three_phase = opp_three_phase.Ignore; %how to deal with three-phase considerations
        Unipolar(1, 1)logical = false %unipolar (in case of half-wave or quarter-wave)
        early_stop(1, 1)logical = true %stop at even pulse numbers (instead of k)
        start_level(1, 1)int32 = 0; %start at an initial level (if not zero)


        %polynomial optimization parameters
        partition(1, 1)int32 = 4; %number of partitions of the disc c^2+s^2=1
        solver = 'mosek';
        TIME_INDEP(1, 1)logical = false; %include time as an explicit state

    end

end

