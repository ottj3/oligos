package edu.tcnj.oligos.data;

import com.google.common.collect.Lists;

import java.util.List;

public enum Base {
    A,
    C,
    T,
    G,

    W(A, T),
    S(C, G),
    M(A, C),
    K(G, T),
    R(A, G),
    Y(C, T),

    B(C, G, T, S, K, Y),
    D(A, G, T, R, W, K),
    H(A, C, T, M, W, Y),
    V(A, C, G, M, R, S),

    N(A, C, T, G, W, S, M, K, R, Y, B, D, H, V);

    List<Base> included;

    Base(Base... included) {
        if (included == null) {
            this.included = Lists.newArrayList(this);
        } else {
            this.included = Lists.newArrayList(included);
        }
    }

    public boolean matches(Base other) {
        if (this == other || this == N || other == N) {
            return true;
        } else if (this.included.contains(other) || other.included.contains(this)) {
            return true;
        }
        return false;
    }
}