import numpy as np

class Civ:
    def __init__(self, alpha, beta, gamma, d=2, capital=None, color=np.array([1,0,0])):
        self.alpha = alpha
        self.beta = beta
        self.gamma = gamma
        self.d = d
        self.color = color
        self.capital = capital

    def c(self, P, R):
        return self.alpha*P

    def p(self, P, R):
        if self.alpha*P == 0:
            return 0
        return self.beta*R

    def __eq__(self, other):
        return (self.color == other.color).all()

    def __ne__(self, other):
        return (self.color != other.color).any()
