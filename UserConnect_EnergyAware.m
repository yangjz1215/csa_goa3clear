% UserConnect_EnergyAware.m - 基于任务优先级的连接算法
% 参考 Code_main0\UserConnect.m
% UAV_type: 0=普通, 1=增强型(eUAV)
% RRH_type: 0=普通, 1=增强型(eRRH)
function Connect = UserConnect_EnergyAware(User, N_User, RRH, N_RRH, RRH_type, UAV, N_UAV, UAV_type, priorities, E_remaining, E_max, params)
    H = 50;
    omega0 = 1e6;
    alpha0 = 1.42e-4;
    sigma2 = 3.98e-12;
    
    D_RRH = 150;
    D_UAV = 150;
    
    if nargin >= 12 && ~isempty(params)
        if isfield(params, 'cover_radius')
            D_UAV = params.cover_radius;
        end
        if isfield(params, 'RRH_radius')
            D_RRH = params.RRH_radius;
        end
    end
    
    Ptx = 1;
    PtxR = 10;
    PtxU = 10;
    PtxEU = 5;
    
    N_CR = 10;
    N_CU = 8;
    N_CEU = 6;
    N_ERT = 10;
    N_EUT = 5;
    
    Dis_UserRRH = zeros(N_User, N_RRH);
    Dis_UserUAV = zeros(N_User, N_UAV);
    for i = 1:N_User
        Dis_UserRRH(i, :) = sqrt(sum((User(i, :) - RRH) .^ 2, 2));
        Dis_UserUAV(i, :) = sqrt(sum((User(i, :) - UAV) .^ 2, 2) + H^2);
    end
    
    Dis_RRHBBU = zeros(N_RRH, 1);
    Dis_RRHBBU = sqrt(sum(RRH.^2, 2));
    Dis_UAVBBU = zeros(N_UAV, 1);
    Dis_UAVBBU = sqrt(sum(UAV.^2, 2) + H^2);
    
    Rate_UserRRH = zeros(N_User, N_RRH);
    Rate_UserRRH = omega0*log2(1 + alpha0*Ptx./(sigma2*(Dis_UserRRH.^2)));
    
    Rate_UserUAV = zeros(N_User, N_UAV);
    Rate_UserUAV = omega0*log2(1 + alpha0*Ptx./(sigma2*(Dis_UserUAV.^2)));
    
    Rate_RRHBBU = zeros(N_RRH, 1);
    Rate_RRHBBU = omega0*log2(1 + alpha0*PtxR./(sigma2*(Dis_RRHBBU.^2)));
    
    Rate_UAVBBU = zeros(N_UAV, 1);
    for i = 1:N_UAV
        if UAV_type(i) == 1
            Rate_UAVBBU(i) = omega0*log2(1 + alpha0*PtxEU/(sigma2*(Dis_UAVBBU(i)^2)));
        else
            Rate_UAVBBU(i) = omega0*log2(1 + alpha0*PtxU/(sigma2*(Dis_UAVBBU(i)^2)));
        end
    end
    
    Connect = zeros(1, N_User);
    Range = [D_RRH*ones(1,N_RRH), D_UAV*ones(1,N_UAV)];
    
    [~, index] = sort(priorities, 'descend');
    Dis = [Dis_UserRRH, Dis_UserUAV];
    
    UAV_type_converted = zeros(size(UAV_type));
    UAV_type_converted(UAV_type == 1) = 1;
    types = [RRH_type(:)', UAV_type_converted(:)'];
    
    for i = 1:N_User
        user = index(i);
        ES_index = find(types == 1);
        valid_index = ES_index(Dis(user, ES_index) <= Range(ES_index));
        
        if ~isempty(valid_index)
            tmpDist = Dis(user, valid_index);
            [~, tindex] = sort(tmpDist);
            
            for kk = 1:length(valid_index)
                j = valid_index(tindex(kk));
                if j <= N_RRH
                    if sum(Connect == j) < N_CR && sum(Connect == j) < N_ERT
                        Connect(user) = j;
                        break;
                    end
                else
                    uav_index = j - N_RRH;
                    if sum(Connect == j) < N_CEU && sum(Connect == j) < N_EUT
                        Connect(user) = j;
                        break;
                    end
                end
            end
        end
        
        if Connect(user) == 0
            tmp_dis = Dis(user, :);
            tmp_dis(types == 1) = Inf;
            [~, tmp_index] = sort(tmp_dis);
            
            for k = 1:N_RRH + N_UAV
                j = tmp_index(k);
                if tmp_dis(j) <= Range(j)
                    if j <= N_RRH
                        if sum(Connect == j) < N_CR
                            Connect(user) = j;
                            break;
                        end
                    else
                        uav_index = j - N_RRH;
                        if sum(Connect == j) < N_CU
                            Connect(user) = j;
                            break;
                        end
                    end
                end
            end
        end
    end
end
