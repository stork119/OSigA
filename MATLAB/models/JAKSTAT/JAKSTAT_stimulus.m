function result = JAKSTAT_stimulus(t,i)
 

% 
%  x=[-1,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,20,30,40,50,60,80,100];
%  y=[0,0,0,0,0,0,0,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0];
     
 x=[0,1,2,3,4,5];
 y=[0,1,1,1,1,0];
     
 % if t<10 z =0;
 % elseif t > 15 z = 0;
 % else z = interp1(x,y,t,'linear');
      
   
    
 stim={};
    stim{1}=(t>0&t<5)*interp1(x,y,t,'pchip');
    result=stim{i};    
   % result = 1;
end



% tmesh_stm = 0:0.01:5;
% stm_res = [];
% for ti = 1:max(size(tmesh_stm))
%     t = tmesh_stm(ti);
%     stm_res(ti) = JAKSTAT_stimulus(t,1);
% end
%  csvwrite([pathoutputfolder{1},'/interpolation.csv'],stm_res);
