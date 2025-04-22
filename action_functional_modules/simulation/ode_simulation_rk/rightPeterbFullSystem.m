% Для численного интегрирования модели объекта управления замкнутого
% обратной связью
% x' = Ax + Bu + Gw; (w - noise)

function dx = rightPeterbFullSystem(t, x, A,B,G, cntrl,ptb)
% dx = A * x + B*cntrl + G*ptb;

dx = A * x + cntrl + ptb;
end

