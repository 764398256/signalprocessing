function [micOut, erle] = Apa(micSpec, spkSpec, nTaps, mu)

L = 2;
N = nTaps;
r = 0.000001;
threshold = 0.032;

assert(size(micSpec,1) == size(spkSpec,1));
nBands = size(micSpec,1);
nFrames = min(size(micSpec, 2), size(spkSpec,2));

% pre-allocate memory
micOut = complex(zeros(nBands, nFrames));
w = complex(zeros(N, nBands));   % filter
Xap = complex(zeros(N, L));

spkHistory = complex(zeros(L + N -1, nBands));
Dap = complex(zeros(L, nBands));
Eap = complex(zeros(L, nBands));

erle = zeros(1, nFrames);

for (ii = 1:nFrames)
    spkHistory(2:end,:) = spkHistory(1:end-1,:);
    spkHistory(1, :) = spkSpec(:, ii).';
    Dap(2:end,:) = Dap(1:end-1,:);
    Dap(1, :) = micSpec(:, ii).';
    
    for i = 1:nBands
        for j = 1:N;
            Xap(j,:) = spkHistory(j:j+L-1, i);
        end;
        Eap(:,i) = conj(conj(Dap(:, i)) - Xap'*w(:,i));
        if (sum(spkHistory(1:N,i).*conj(spkHistory(1:N,i)))/N > threshold^2)
            temp =  Xap * inv(Xap' * Xap + r.*eye(L)) * (conj(Eap(:,i)));
            w(:,i) = w(:,i) + mu * temp;
        end;
    end;
    
    micOut(:, ii) = Eap(1,:).';
    erle(ii) = db(sum(micSpec(:, ii).*conj(micSpec(:, ii)))/...
        (sum(micOut(:,1).*conj(micOut(:,1))) + 1.0E-10))/2;
    if mod(ii, 100) == 0
        fprintf('%d %d\n', ii, erle(ii));
    end;
end;