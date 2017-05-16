"""
This module computes the probabilities and functions defined in the article, for binary trees only.
In order to do that, we need to import these modules.
Variable `a` represents the alpha in Ford's Alpha model.
"""
from __future__ import division
import Shape, PhyloTree, newick
from sage.all import *
from sympy import *

a = var('a')

def Pi(T):
    """
    Compute the product in Proposition 2.
    :param T: A given `Shape` object (it depends only on the topology of T).
    :return: A symbolic expression depending on `a`.
    """
    if T.is_leaf or T.children == []:
        return 1
    elif len(T.children) > 1:
        b = T.children[0].count_leaves()
        c = T.children[1].count_leaves()
        if not b*c == 0:
            return phi(b, c)*prod(Pi(x) for x in T.children)
        else:
            return 1
    else:
        return 1

def phi(n, m):
    """
    Compute phi as defined in the article.
    :param n, m: Given integers.
    :return: A symbolic expression depending on `a`.
    """
    return a/Integer(2)*binomial(n+m, n) + (1-2*a)*binomial(n+m-2, n-1)

def Gamma_alpha(n):
    """
    Compute the Gamma-Alpha function of n.
    :param n: Given integer.
    :return: A symbolic expression depending on `a`.
    """
    if n == 2:
        return 1 - a
    elif n > 2:
        return (n - 1 - a)*Gamma_alpha(n-1)
    else:
        return 0

def prob_shape(T):
    """
    Compute the probability of generating the shape of T under the Alpha model by Ford, assuming T is an instance of `Shape` --i.e.,
    that it is only the shape of a tree.
    :param T: The instance of `Shape` of which we want to compute the probability.
    :return: A symbolic expression depending on `a`.
    """
    if not bool(T.children):
        return 1
    T1 = T.shape()
    n  = T1.count_leaves()
    Gamma = Gamma_alpha(n)
    if Gamma == 0:
        return "ERROR: division by zero"
    else:
        k = T1.count_symmetries()
        eq = 2 ** (n - k - 1) / Gamma * Pi(T)
        return eq.factor()

def prob_tree(T):
    """
    Compute the probability of generating T under the Alpha model by Ford, assuming T is an instance of `PhyloTree`. 
    Given an instance of `Shape`, it computes the probability of a phylogenetic tree with its shape.
    :param T: The instance of `PhyloTree` of which we want to compute the probability.
    :return: A symbolic expression depending on `a`.
    """
    if not bool(T.children):
        return 1
    n = T.count_leaves()
    Gamma = Gamma_alpha(n)
    if Gamma == 0:
        return "ERROR: division by zero"
    else:
        eq = 2 ** (n - 1) / (factorial(n)*Gamma) * Pi(T)
        return eq.factor()

def prob(T):
    """
    Compute the probability of generating T under the Alpha model by Ford. 
    Given an instance of `Shape`, it computes the probability of generating the given shape.
    Given an instance of `PhyloTree` that is a phylogenetic tree --i.e., such that does not have repeated labels--, it 
    computes the probability of generating the given tree.
    :param T: The instance of `Shape` of which we want to compute the probability.
    :return: A symbolic expression depending on `a`.
    """
    if not bool(T.children):
        return 1
    n = T.count_leaves()
    Gamma = Gamma_alpha(n)
    if Gamma == 0:
        return "ERROR: division by zero"
    else:
        k = T.count_symmetries()
        nlabels = len(T.labels())
        eq = 2**(n-k-1) / (factorial(nlabels)*Gamma) * Pi(T)
        return eq.factor()

def shape_probs_from_list(X):
    """
    Return the probabilities of all the instances of `Shape` in a list (in string) of Newick codes.
    :param X: a string containing several Newick codes.
    :return: A list of symbolic expressions depending on `a`.
    """
    l = Shape.from_newick_list(X)
    return [prob_shape(t) for t in l]

def tree_probs_from_list(X):
    """
    Return the probabilities of all the instances of `PhyloTree` in a list (in string) of Newick codes.
    :param X: a string containing several Newick codes.
    :return: A list of symbolic expressions depending on `a`.
    """
    l = PhyloTree.from_newick_list(X)
    return [prob_tree(t) for t in l]

def shape_probs_from_file(fname, encoding='utf8', strip_comments=False, **kw):
    """
    Load a list of instances of `Shape` from a Newick formatted file and return their probabilities in a list, assuming
    they are instances of `Shape`.
    :param fname: file path.
    :param strip_comments: Flag signaling whether to strip comments enclosed in square \
    brackets.
    :param kw: Keyword arguments are passed through to `Node.read`.
    :return: A list of symbolic expressions depending on `a`.
    """
    l = newick.read(fname, encoding, strip_comments, **kw)
    return [prob_shape(Shape.newick_node_to_shape(t)) for t in l]

def tree_probs_from_file(fname, encoding='utf8', strip_comments=False, **kw):
    """
    Load a list of instances of `PhyloTree` from a Newick formatted file and return their probabilities in a list, assuming
    they are instances of `PhyloTree`.
    :param fname: file path.
    :param strip_comments: Flag signaling whether to strip comments enclosed in square \
    brackets.
    :param kw: Keyword arguments are passed through to `Node.read`.
    :return: A list of symbolic expressions depending on `a`.
    """
    l = newick.read(fname, encoding, strip_comments, **kw)
    return [prob_tree(PhyloTree.newick_node_to_tree(t)) for t in l]