classdef subsystem_base < subsystem_interface
    %SUBSYSTEM A subsystem x'=f(t, x) of dynamical system in measure 
    %analysis where there is no uncertainty
    
    properties
        
        
        %additional measures for box uncertainty
        meas_box = {};  %box occupation measures
        meas_comp = {}; %box-complement occupation measures                
        
        varnames = {'t', 'x'};
        
        f_box = {};     %affine decomposition of dynamics 
                        %{no input, input 1, input 2, ...}                                                       
    end

    
    methods
        %% Constructor
        function obj = subsystem_base(loc_supp, f, sys_id, loc_id)
            %SUBSYSTEM Construct a continuous (possibly uncertain) 
            %subsystem, fill in information
            
            %process input
            if nargin < 3
                sys_id = 1;
            end
            
            if nargin < 4
                loc_id = [];
            end
            
            %superclass constructor
            obj@subsystem_interface(loc_supp, f, sys_id, loc_id);
                                                                                                 
        end
               

        %% Constraints        
        
       
        function Ay = cons_liou(obj, d)
            %CONS_LIOU Liouville Equation includes an affine combination of
            %Lie derivatives (continuous systems only)
            
            %no box inputs, simple to perform
             Ay = obj.meas_occ.mom_lie(d, obj.vars, obj.f);
           
            
        end       
       
        

        
        
        %% Dual Recovery 
        function obj = dual_process(obj, v, zeta)
            %DUAL_PROCESS store dual functions and compute nonnegative
            %functions for this subsystem
            
            %auxiliary function v            
            obj.dual.v = v;
            

            obj.dual.Lv = diff(v, obj.vars.x)*obj.f;
            
            if ~isempty(obj.vars.t)
                obj.dual.Lv = obj.dual.Lv + diff(v, obj.vars.t);
            end
            
            
            obj.nn_ = obj.dual.nn;
        end        
        
        %% Function Evaluation (for sampling)
        %TODO: Implement this
        %nonnegativity evaluation: include supports        
        

    end
end

