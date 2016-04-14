%function Decoder = Viterbi(ConCode,Length)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%该函数实现维特比译码
%%%%ConCode为待译的卷积码
%%%%Length为原始序列的长度
%%%%Decoder为译码后所得码字
%%%%S0->00,S1->10,S2->01,S3->11
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    RecSeq = ConCode;
    Decoder = zeros(1,Length);      %存储维特比译码后的码字
    D = zeros(1,4);                 %记录到达某一状态时的累计汉明距离
    %DTemp = D;
    RouteS0 = zeros(1,2*Length);    %记录到达S0状态所经过的路径
    RouteS1 = zeros(1,2*Length);    %记录到达S1状态所经过的路径
    RouteS2 = zeros(1,2*Length);    %记录到达S2状态所经过的路径
    RouteS3 = zeros(1,2*Length);    %记录到达S3状态所经过的路径                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            = zeros(1,2*Length);
  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%首先得到最优路径
    for k = 1 : Length
        if k == 1                   %从S0状态开始译码，第一次只需记录从S0状态到
            D(1) = RecSeq(1) + RecSeq(2);%S1状态和S2状态时的汉明距离
            D(2) = abs(RecSeq(1)-1) + abs(RecSeq(2)-1);
        elseif k == 2                %在第二时刻到达各个状态的路径分别为S0->S0
            DTemp = D;               %S0->S1,S1->S2,S1->S3
            DS0_S0 = RecSeq(3) + RecSeq(4);
            DS0_S1 = abs(RecSeq(3)-1) + abs(RecSeq(4)-1);
            DS1_S2 = abs(RecSeq(3)-1) + RecSeq(4);
            DS1_S3 = RecSeq(3) + abs(RecSeq(4)-1);
            D(1) = DTemp(1) + DS0_S0;
            D(2) = DTemp(1) + DS0_S1;
            D(3) = DTemp(2) + DS1_S2;
            D(4) = DTemp(2) + DS1_S3;
            RouteS0(1:4) = [0 0 0 0];
            RouteS1(1:4) = [0 0 1 0];
            RouteS2(1:4) = [1 0 0 1];
            RouteS3(1:4) = [1 0 1 1];
        else
            %%%%从第3时刻开始出现路径的选择
            %%计算从上一时刻到该时刻的一步汉明距离，到达每一状态都有两条路径
            DS0_S0 = RecSeq(2*k-1) + RecSeq(2*k);  
            DS2_S0 = abs(RecSeq(2*k-1)-1) + abs(RecSeq(2*k)-1);
            DS0_S1 = DS2_S0;
            DS2_S1 = DS0_S0;
            DS1_S2 = abs(RecSeq(2*k-1)-1) + RecSeq(2*k);
            DS3_S2 = RecSeq(2*k-1) + abs(RecSeq(2*k)-1);
            DS1_S3 = DS3_S2;
            DS3_S3 = DS1_S2;
            
            %%保存原来的路径及累计汉明距离
            RS0Temp = RouteS0(1:2*k-2);
            RS1Temp = RouteS1(1:2*k-2);
            RS2Temp = RouteS2(1:2*k-2);
            RS3Temp = RouteS3(1:2*k-2);
            DTemp = D;
            
            %%通过比较来得到新的路径和累计汉明距离
            if (DTemp+DS0_S0 <= DTemp(3)+DS2_S0)
                D(1) = DTemp(1) + DS0_S0;
                RouteS0(1:2*k) = [RS0Temp 0 0];
            else
                D(1) = DTemp(3) + DS2_S0;
                RouteS0(1:2*k) = [RS2Temp 0 0];
            end
            
            if (DTemp(1)+DS0_S1 <= DTemp(3)+DS2_S1)
                D(2) = DTemp(1) + DS0_S1;
                RouteS1(1:2*k) = [RS0Temp 1 0];
            else
                D(2) = DTemp(3) + DS2_S1;
                RouteS1(1:2*k) = [RS2Temp 1 0];
            end
            
            if (DTemp(2)+DS1_S2 <= DTemp(4)+DS3_S2)
                D(3) = DTemp(2) + DS1_S2;
                RouteS2(1:2*k) = [RS1Temp 0 1];
            else
                D(3) = DTemp(4) + DS3_S2;
                RouteS2(1:2*k) = [RS3Temp 0 1];
            end
            
            if (DTemp(2)+DS1_S3 <= DTemp(4)+DS3_S3)
                D(4) = DTemp(2) + DS1_S3;
                RouteS3(1:2*k) = [RS1Temp 1 1];
            else
                D(4) = DTemp(4) + DS3_S3;
                RouteS3(1:2*k) = [RS3Temp 1 1];
            end
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%根据最优路径及状态转移图来译码得到码字
    
   a00 = [0 0];a01 = [0 1];
   a11 = [1 1];a10 = [1 0];
    for m = 1 : Length
        if (m == 1)
            if (RouteS0(1) == 0)
                    Decoder(1) = 0;
            else
                    Decoder(1) = 1;
            end
        elseif (m == 2)
                if (RouteS0(3) == 0)
                    Decoder(2) = 0;
                else
                    Decoder(2) = 1;
                end
        else
            L1 = 2*m - 3;
            R1 = 2*m - 2;
            L2 = 2*m - 1;
            R2 = 2*m;
            if ((RouteS0(L1:R1)==a00||RouteS0(L1:R1)==a11) && RouteS0(L2:R2)==a00)
                Decoder(m) = 0;
            elseif ((RouteS0(L1:R1)==a10||RouteS0(L1:R1)==a11) && RouteS0(L2:R2)==a01)
                Decoder(m) = 0;
            else
                Decoder(m) = 1;
            end
        end
    end