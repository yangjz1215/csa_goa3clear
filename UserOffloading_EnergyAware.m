% UserOffloading_EnergyAware.m - 基于任务优先级的卸载算法
% 参考 Code_main0\UserOffloading.m
% 注意：Connect数组中，1到N_RRH是RRH索引，N_RRH+1到N_RRH+N_UAV是UAV索引
% UAV_type: 0=普通, 1=增强型(eUAV)
% RRH_type: 0=普通, 1=增强型(eRRH)
function off_decision = UserOffloading_EnergyAware(User, N_User, N_RRH, RRH, RRH_type, UAV, N_UAV, UAV_type, ...
    D, C, priorities, Connect, E_remaining, E_max, params)
    
    H = 50;
    omega0 = 1e6;
    alpha0 = 1.42e-4;
    sigma2 = 3.98e-12;
    
    PtxU = 10;
    PtxEU = 5;
    ki = 1e-27;
    
    N_ERT = 10;
    N_EUT = 5;
    
    f_eUAV = 2e9;
    Pho = 100;
    
    if isscalar(E_max)
        E_max_array = E_max * ones(N_UAV, 1);
    else
        E_max_array = E_max(:);
    end
    
    Dis_UAVBBU = zeros(N_UAV, 1);
    Dis_UAVBBU = sqrt(sum(UAV.^2,2) + H^2);
    
    Rate_UAVBBU = zeros(N_UAV, 1);
    for i = 1:N_UAV
        if UAV_type(i) == 1
            Rate_UAVBBU(i) = omega0*log2(1+alpha0*PtxEU/(sigma2*(Dis_UAVBBU(i)^2)));
        else
            Rate_UAVBBU(i) = omega0*log2(1+alpha0*PtxU/(sigma2*(Dis_UAVBBU(i)^2)));
        end
    end
    
    off_decision = zeros(1, N_User);
    UAV_energy = E_remaining(:)';
    UAV_Time = zeros(1, N_UAV);
    
    [~, index] = sort(priorities, 'descend');
    
    for i = 1:N_User
        user = index(i);
        if Connect(user) == 0
            continue;
        else
            device_idx = Connect(user);
            
            if device_idx <= N_RRH
                rrh_idx = device_idx;
                if RRH_type(rrh_idx) == 1
                    if sum(off_decision == rrh_idx) < N_ERT
                        off_decision(user) = rrh_idx;
                    else
                        off_decision(user) = N_RRH + N_UAV + 1;
                    end
                else
                    off_decision(user) = N_RRH + N_UAV + 1;
                end
            else
                uav_idx = device_idx - N_RRH;
                
                if UAV_type(uav_idx) == 1
                    if isvector(C) && length(C) > 1
                        C_user = C(user);
                    else
                        C_user = C;
                    end
                    Energy = ki*(f_eUAV^2)*C_user;
                    
                    h0 = alpha0/sum((UAV(uav_idx,:)-User(user,:)).^2 + H^2);
                    r = omega0*log2(1+(PtxEU*h0)/sigma2);
                    
                    if isvector(D) && length(D) > 1
                        D_user = D(user);
                    else
                        D_user = D;
                    end
                    t_tx = D_user/r;
                    t_ex = C_user/f_eUAV;
                    
                    if sum(off_decision == device_idx) < N_EUT && ...
                       (Pho*max(UAV_Time(uav_idx),t_tx+t_ex)+Energy) <= UAV_energy(uav_idx)
                        off_decision(user) = device_idx;
                        UAV_energy(uav_idx) = max(0, UAV_energy(uav_idx) - Energy);
                        UAV_Time(uav_idx) = max(UAV_Time(uav_idx), t_tx+t_ex);
                    else
                        Tx = D_user/Rate_UAVBBU(uav_idx);
                        Energy_tx = PtxEU*Tx;
                        if (Pho*max(UAV_Time(uav_idx),t_tx+Tx)+Energy_tx) <= UAV_energy(uav_idx)
                            off_decision(user) = N_RRH + N_UAV + 1;
                            UAV_energy(uav_idx) = max(0, UAV_energy(uav_idx) - Energy_tx);
                            UAV_Time(uav_idx) = max(UAV_Time(uav_idx), t_tx+Tx);
                        end
                    end
                else
                    h0 = alpha0/sum((UAV(uav_idx,:)-User(user,:)).^2 + H^2);
                    r = omega0*log2(1+(PtxU*h0)/sigma2);
                    
                    if isvector(D) && length(D) > 1
                        D_user = D(user);
                    else
                        D_user = D;
                    end
                    t_tx = D_user/r;
                    Tx = D_user/Rate_UAVBBU(uav_idx);
                    Energy_tx = PtxU*Tx;
                    
                    if (Pho*max(UAV_Time(uav_idx),t_tx+Tx)+Energy_tx) <= UAV_energy(uav_idx)
                        off_decision(user) = N_RRH + N_UAV + 1;
                        UAV_energy(uav_idx) = max(0, UAV_energy(uav_idx) - Energy_tx);
                        UAV_Time(uav_idx) = max(UAV_Time(uav_idx), t_tx+Tx);
                    end
                end
            end
        end
    end
end
