function [AIneq, BIneq] = constraint(P, Ju, Jv, x, deltaF, disp_min, disp_max)

% y = Px + Ju deltaU + Jv deltaF <= ymax
%=> Ju deltaU <= ymax - Px -Jv deltaF
% where, y = [  displacement  ;  velocity  ;  control_force/a  ];

[Pn, Pm] = size(P);

BMaxIneqTotal1 = disp_max(1,1)*ones(Pn,1) - P*x - Jv*deltaF;    % displacement constraint
BMaxIneq1 = BMaxIneqTotal1(1:3:end,1);
AMaxIneq1 = Ju(1:3:end,:);

BMaxIneqTotal2 = disp_max(2,1)*ones(Pn,1) - P*x - Jv*deltaF;    % velocity constraint
BMaxIneq2 = BMaxIneqTotal2(2:3:end,1);
AMaxIneq2 = Ju(2:3:end,:);

BMaxIneqTotal3 = disp_max(3,1)*ones(Pn,1) - P*x - Jv*deltaF;    % control force constraint
BMaxIneq3 = BMaxIneqTotal3(3:3:end,1);
AMaxIneq3 = Ju(3:3:end,:);

% y = Px + Ju deltaU + Jv deltaF >= ymin
%=> -Ju deltaU <= -ymin + Px  + Jv deltaF
BMinIneqTotal1 = -disp_min(1,1)*ones(Pn,1) + P*x + Jv*deltaF;   % displacement constraint
BMinIneq1 = BMinIneqTotal1(1:3:end,1);
AMinIneq1 = -Ju(1:3:end,:);

BMinIneqTotal2 = -disp_min(2,1)*ones(Pn,1) + P*x + Jv*deltaF;   % velocity constraint
BMinIneq2 = BMinIneqTotal2(2:3:end,1);
AMinIneq2 = -Ju(2:3:end,:);

BMinIneqTotal3 = -disp_min(3,1)*ones(Pn,1) + P*x + Jv*deltaF;   % control force constraint
BMinIneq3 = BMinIneqTotal3(3:3:end,1);
AMinIneq3 = -Ju(3:3:end,:);

AIneq = [];
BIneq = [];

if (disp_min(1,1) ~= 0 || disp_max(1,1) ~= 0)
    AIneq = [  AMinIneq1  ;  AMaxIneq1  ];
    BIneq = [  BMinIneq1  ;  BMaxIneq1  ];
end
if (disp_min(2,1) ~= 0 || disp_max(2,1) ~= 0)
    AIneq = [  AIneq  ;  AMinIneq2  ;  AMaxIneq2  ];
    BIneq = [  BIneq  ;  BMinIneq2  ;  BMaxIneq2  ];
end
if (disp_min(3,1) ~= 0 || disp_max(3,1) ~= 0)
    AIneq = [  AIneq  ;  AMinIneq3  ;  AMaxIneq3  ];
    BIneq = [  BIneq  ;  BMinIneq3  ;  BMaxIneq3  ];
end