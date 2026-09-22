function [q, qd] = renormalize_quaternion(q, qd)
%RENORMALIZE_QUATERNION  Unit-length correction of Euler parameters after a step.
%   Position: q <- q/||q||, Terze et al., Multibody Syst Dyn 38 (2016), App. C, (C.7).
%   If qd is requested, also qdot <- qdot - (qdot'*q)*q so that q'*qdot = 0, (C.8).

    q = q / norm(q);
    if nargout > 1
        qd = qd - (qd.' * q) * q;
    end
end
