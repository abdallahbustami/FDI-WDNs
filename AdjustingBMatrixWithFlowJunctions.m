    %%This function updates B matrix for chlorine according to the booster
    %%station location
    function[BAo]=AdjustingBMatrixWithFlowJunctions(JJ,BAo,q_B,Location_B,qout)
    [~,nB] = size(Location_B);
    for bb= 1:nB
        if JJ== Location_B(bb)
            %   BAo(JJ,bb)=q_B(JJ)/(qout+q_B(JJ));
            %         BAo(JJ,bb)=q_B(JJ)/(qout);
            BAo(JJ,bb)=q_B(JJ);
        end
    end
    end