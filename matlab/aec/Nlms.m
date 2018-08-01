function [cleSpec] = Nlms(micSpec, spkSpec, nTaps, mu)

threshold = 0.032;

[nBands, nFrames] = size(micSpec);
nFrames = min(nFrames, size(spkSpec, 2));
w = zeros(nTaps, nBands);
cleSpec = zeros(nBands, nFrames);
power = zeros(nBands, nFrames);
erle = zeros(nFrames);

for iFrames = 1 + nTaps:nFrames
    for iBands = 1:nBands
        x = spkSpec(iBands, iFrames + (-nTaps:-1));
        x = x.';
        d = micSpec(iBands, iFrames);
        d = d.';
        e = d - w(:, iBands)' * x;
        
        power(iBands, iFrames) = x' *  x;
        
        if (power(iBands, iFrames)/nTaps > threshold^2)
            w(:, iBands) = w(:, iBands) + mu * conj(e) * x./power(iBands, iFrames);
        end;
        
        cleSpec(iBands, iFrames) = e;
    end;
    if(mod(iFrames, 100) == 0)
        fprintf('%d\n', iFrames);
    end;
    erle(iFrames) = db(sqrt(sum(micSpec(:,iFrames).*conj(micSpec(:,iFrames)))...
        /sum(cleSpec(:,iFrames).*conj(cleSpec(:,iFrames)))) + 1e-10);
end;

if(0)
close all;
figure(1);
plot(erle);
end;

end