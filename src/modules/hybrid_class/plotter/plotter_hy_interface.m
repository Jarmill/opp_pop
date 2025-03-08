classdef plotter_hy_interface < handle
    %PLOTTER_HY_INTERFACE an interface to plotting functions for
    %trajectories of hybrid systems
    
    properties
        FS_axis = 14;       %font size of axis
        FS_title = 16;      %font size of title        
    end
    
    properties(Access=protected)
        %OSM: out_sim_multi:    a cell array of trajectories
        %OSD: out_sim_deal:     trajectories split into locations and
        %                       guards
        %the two output of sampler_hy.sample_traj_multi
        %
        %hybrid trajectories are 'dealt' into their residing locations and
        %guards. osd is a struct with fields 'locations' and 'guards'
        osm; 
        osd;
    end
    
    methods
        function obj = plotter_hy_interface(osm, osd)
            %PLOTTER_HY_INTERFACE Construct an instance of this class
            obj.osm = osm;
            obj.osd = osd;
        end
              
        
        %% trajectory setter and getter
        function obj = set_osd(obj, osd_new)
            obj.osd = osd_new;
        end
        
        function obj = set_osm(obj, osm_new)
            obj.osm = osm_new;
        end
                        
        function [osm_out, osd_out] = get_traj(obj)
            osm_out = obj.osm;
            osd_out = obj.osd;
        end
        
        
        %% plot nonnegative functions along trajectories
        function F = nonneg_jump(obj)
            %NONNEG_JUMP plot nonnegative transitions from the jump            
            F = figure(60);
            clf
            hold on
            title('Jump Nonnegativity', 'FontSize', obj.FS_title)
            ylabel("$v_i(x) - v(R_{i \rightarrow i'} x)$", 'interpreter', 'latex', 'FontSize', obj.FS_axis);
            xlabel('time', 'FontSize', obj.FS_axis)
            for i = 1:length(obj.osd.guards)
                os_guards = obj.osd.guards{i};
                for j = 1:length(os_guards)
%                 osn = obj.osm{i};
%                 for j = 1:length(osn.jump)
                    j_curr = os_guards{j};            
                    stem(j_curr.t, j_curr.nonneg, 'c') 
                end
            end
        end 
        
                
        function F = nonneg_loc(obj)
            %NONNEG_LOC plot nonnegative functions in the locations
            
            %TODO: update this with the new nonnegative evaluators
            %also where there are uncertainty processes, and therefore a
            %varying number of nonnegative functions
            F = figure(61);
            clf
            nonneg_title = {'Initial Value', 'Cost Proxy', 'Decrease in Value', };
            ax_loc = {'$\gamma - v(x)$', '$v(x) - p(x)$', '$-L_{f} v(x)$'};
            for k = 1:3
                subplot(3, 1, k)
                hold on
                xlabel('time', 'FontSize', obj.FS_axis)
                ylabel(ax_loc{k}, 'interpreter', 'latex', 'FontSize', obj.FS_axis);
                title(nonneg_title{k}, 'FontSize', obj.FS_title);   
                
                %TODO: write code determining which nonneg should be blank
                
                for i = 1:length(obj.osd.locations)
                    os_loc = obj.osd.locations{i};
                    for j = 1:length(os_loc)
                        os_curr = os_loc{j};
                        plot(os_curr.t, os_curr.nonneg(k, :)', 'c');
                    end
                end
            end

        end
        
        
         function F = aux_plot(obj)
            %aux_plot plot the auxiliary function v along trajectories (in
            %all locations)
            
            %TODO: update this with the new nonnegative evaluators
            %also where there are uncertainty processes, and therefore a
            %varying number of nonnegative functions
            F = figure(62);
            clf
%             nonneg_title = {'Initial Value', 'Decrease in Value', 'Cost Proxy'};
%             ax_loc = {'$\gamma - v(x)$', '$v(x) - p(x)$', '$-L_{f} v(x)$'};

                hold on
                xlabel('time', 'FontSize', obj.FS_axis)
                ylabel('$v_i(x)$', 'interpreter', 'latex', 'FontSize', obj.FS_axis);
                title('Auxiliary Functions v', 'FontSize', obj.FS_title);   
                
                for i = 1:length(obj.osd.locations)
                    os_loc = obj.osd.locations{i};
                    for j = 1:length(os_loc)
                        os_curr = os_loc{j};
                        plot(os_curr.t, os_curr.v', 'c');
                    end
                end

        end
        
        %% plot trajectories and properties
        
                
        function F = objective_plot(obj, p_est)
            %OBJECTIVE_PLOT
            F = figure(65);
            clf
            hold on
            for i = 1:length(obj.osm)
                osn = obj.osm{i};
                for j = 1:length(osn.sim)
                    traj_curr = osn.sim{j};
                    if ~isempty(traj_curr.objective)
                        plot(traj_curr.t, traj_curr.objective, 'c')
                    end
%                     scatter(traj_curr.t(end), traj_curr.x(end, 2).^2, 30, 'ok')
                end
            end
                
            if nargin == 2
                plot(xlim, p_est*[1;1], 'r--', 'LineWidth', 2)
            end
        %     xlabel('time', 'FontSize', FS_axis);
            xlabel('$t$', 'interpreter', 'latex', 'FontSize', obj.FS_axis);
            ylabel('$p(x)$', 'interpreter', 'latex', 'FontSize', obj.FS_axis);
%             zlabel('$\omega$', 'interpreter', 'latex', 'FontSize', obj.FS_axis);
           title('Objective', 'FontSize', obj.FS_title);   
        end
        
    end
end

