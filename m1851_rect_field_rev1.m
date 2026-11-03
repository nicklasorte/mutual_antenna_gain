function F = m1851_rect_field_rev1(n, mu)
%M1851_RECT_FIELD_REV1 Normalized field pattern of a cos^n aperture taper.
%   F = m1851_rect_field_rev1(n, mu) returns F(mu) normalized to F(0) = 1
%   for the rectangular aperture field distribution f(x) = cos^n(pi*x/2),
%   -1 <= x <= 1, without pedestal (Rec. ITU-R M.1851-2 Table 4), using
%   the finite Fourier transform of M.1851-2 Eq. (1).
%
%   The transform is evaluated exactly by expanding cos^n into complex
%   exponentials: cos^n(a) = 2^-n * sum_k C(n,k) exp(j(n-2k)a), so
%       F(mu) = 2^-n * sum_k C(n,k) * S(mu + (n-2k)*pi/2),  S(b) = sin(b)/b.
%   This reproduces the Table 4 closed forms for n = 0..3 and avoids a
%   transcription error in the published n = 4 expression. Built-in sinc
%   (sin(pi*x)/(pi*x)) supplies S(b) = sinc(b/pi).
%
%   n   taper power, integer 0..4
%   mu  pattern variable (rad), any size

arguments
    n (1,1) double {mustBeInteger, mustBeInRange(n, 0, 4)}
    mu double
end

F = zeros(size(mu));
for k = 0:n
    F = F + nchoosek(n, k) * sinc((mu + (n - 2*k)*pi/2) / pi);
end
F = F / 2^n;
F0 = sum(arrayfun(@(k) nchoosek(n, k) * sinc((n - 2*k)/2), 0:n)) / 2^n;
F = F / F0;
end
