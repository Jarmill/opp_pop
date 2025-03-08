classdef opp_three_phase < uint32
   enumeration
      Ignore (0)    %only deal with single-phase considerations
      Balanced (1)  %impose that the current is symmetric and balanced (i(t)+i(t+2pi/3) + i(t+4pi/3) = 0)
      Floating (2)  %penalize only non-triplen harmonics in the energy definition
   end
end