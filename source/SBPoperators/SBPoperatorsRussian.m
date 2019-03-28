function [D1,HI,H,DISS,Dp,Dm,HIfull] = SBPoperatorsRussian(N,h,order)

switch order
    case 4
        [D1,HI,H,DISS,Dp,Dm] = SBP4_Upwind(N,h);
    case 6
        [D1,HI,H,DISS,Dp,Dm] = SBP6_Upwind(N,h);
    case 8
        [D1,HI,H,DISS,Dp,Dm] = SBP8_Upwind(N,h);
end

HIfull = HI;
HI = HI(1,1);

