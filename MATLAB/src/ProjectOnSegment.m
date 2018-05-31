

function cpq = ProjectOnSegment(c1, p1, q1)


cmp = c1 - p1
qmp = q1 - p1

lambda = (cmp * qmp') / (sum(qmp.^2));
lambda = max(0, min(lambda, 1));

cpq = p1 + lambda * qmp;




end