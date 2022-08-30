function Q_r = Q_right(q)

Q_r = [q(1), -q(2), -q(3), -q(4);...
    q(2), q(1), q(4), -q(3);...
    q(3), -q(4), q(1), q(2);...
    q(4), q(3), -q(2), q(1)];

end
