"""
A `Shape` represents a topological tree. The data structure implemented here is of recursive type: a `Shape` can be either
a leaf or a list of `Shape` objects. Leaves are not distinguishable, but we know that they are leaves.
We choose a sorted shape to be the class representant of all shapes isomorphic to it.
In order to read Newick codes we import the newick module from https://github.com/glottobank/python-newick.
"""

import newick


class Shape(object):
    """
    A `Shape` instance is either a leaf or a list of `Shape` instances that hang from a root.
    """
    def __init__(self, children):
        """
        Create a new `Shape` object. 
        The boolean is_leaf is True if the object is a leaf; it is False otherwise.
        :param children: if not is_leaf, the `Shape` objects which are descendants of the considered object; 
        otherwise, `None`.
        :return: `Shape` instance.
        """
        self.is_leaf  = children is None
        self.children = children
        assert self.is_leaf or len(children) > 0

    def sort(self):
        """
        Sorts self using graduate lexicographical order.
        """
        if not self.is_leaf:
            children.sort()

    def compare(self, T2):
        """
        Compare self with another `Shape` object. We use lexicographical order in order to compare two `Shape` instances. 
        Leaves in this case are indistinguishable. It returns anint c, which is 0 if self and T2 are equal, < 0 if self < T2, 
        and > 0 if self > T2.
        :param T2: the `Shape` object against which we compare self.
        :return: int instance.
        """
        if self.is_leaf and T2.is_leaf:
            return 0
        elif self.is_leaf:
            return -1
        elif T2.is_leaf:
            return 1
        else:
            c = len(self.children) - len(T2.children)
            if c != 0:
                return c

            for i in range(0, len(self.children)):
                c = self.children[i].compare(T2.children[i])
                if c != 0:
                    return c
            return 0

    def __lt__(self, T2):
        """
        Uses the comparing method above to decide if self is less than T2.
        :param T2: the `Shape` object against which we compare self.
        :return: bool instance.
        """
        return self.compare(T2) < 0

    def __le__(self, T2):
        """
        Uses the comparing method above to decide if self is less or equal than T2.
        :param T2: the `Shape` object against which we compare self.
        :return: bool instance.
        """
        return self.compare(T2) <= 0

    def __eq__(self, T2):
        """
        Uses the comparing method above to decide if self is equal to T2.
        :param T2: the `Shape` object against which we compare self.
        :return: bool instance.
        """
        return self.compare(T2) == 0

    def __ne__(self, T2):
        """
        Uses the comparing method above to decide if self is not equal to T2.
        :param T2: the `Shape` object against which we compare self.
        :return: bool instance.
        """
        return self.compare(T2) != 0

    def __ge__(self, T2):
        """
        Uses the comparing method above to decide if self is greater or equal than T2.
        :param T2: the `Shape` object against which we compare self.
        :return: bool instance.
        """
        return self.compare(T2) >= 0

    def __gt__(self, T2):
        """
        Uses the comparing method above to decide if self is greater than T2.
        :param T2: the `Shape` object against which we compare self.
        :return: bool instance.
        """
        return self.compare(T2) > 0

    def iso(self, T2):
        """
        Since our `Shape` objects are sorted, to know whether two trees are isomorphic or not it suffices to know whether
        they are equal or not.
        :param T2: the `Shape` object against which we compare self.
        :return: bool instance.
        """
        return self == T2

    def to_newick_tuple(self):
        """
        Returns a tuple representing the simplified Newick code of self. Leaves are marked as 1's, since we do not distinguish them.
        :return: tuple instance.
        """
        if self.is_leaf:
            return 1
        else:
            return tuple(x.to_newick_tuple() for x in self.children)


    def to_newick(self):
        """
        Returns a string representing the simplified Newick code of self.
        :return: string instance.
        """
        return str(self.to_newick_tuple()) + ";"


    def is_symmetric(self):
        """
        Returns True if the root of self is a symmetric node, and False otherwise. If self is a leaf, it returns True:
        ex falso quodlibet.
        :return: int instance.
        """
        if not self.is_leaf:
            return all(self.children[0].iso(x) for x in self.children)
        else:
            return True

    def count_symmetries(self):
        """
        Returns the number of symmetric interior nodes in self.
        :return: int instance.
        """
        if self.is_leaf:
            return 0
        elif all(self.children[0].iso(x) for x in self.children):
            return 1 + sum(x.count_symmetries() for x in self.children)
        else:
            return sum(x.count_symmetries() for x in self.children)

    def count_leaves(self):
        """
        Returns the number of leaves in self.
        :return: int instance.
        """

        if self.is_leaf:
            return 1
        else:
            return sum(x.count_leaves() for x in self.children)

    def shape(self):
        """
        Returns the `Shape` associated to self. Namely, it "forgets" the labels of the leafs.
        :return: `Shape` instance.
        """
        return self

    def labels(self):
        """
        Returns a list with the labels that appear in self, sorted in lexicographical order.
        Since we only use 1's as flags for leaves, the output is [1].
        :return: list instance.
        """
        return [1]



def from_newick(X):
    """
    Create a `Shape` object from a Newick code entered as a string.
    :param X: a string representing a Newick code.
    :return: `Shape` instance.
    """
    return newick_node_to_shape(newick.loads(X)[0])

def from_newick_list(X):
    """
    Create a list of `Shape` objects from a list of Newick codes entered as a string.
    :param X: a string representing a list of Newick codes.
    :return: [`Shape`] instance.
    """
    return[newick_node_to_shape(n) for n in newick.loads(X)]

def newick_node_to_shape(N):
    """
    Create a `Shape` object from a `Node` object.
    :param N: a `Node`.
    :return: `Shape` instance.
    """
    if not bool(N.descendants):
        return Shape(None)
    return Shape(sorted([newick_node_to_shape(x) for x in N.descendants]))


def shapes_from_file(fname, encoding='utf8', strip_comments=False, **kw):
    """
    Load a list of shapes from a Newick formatted file.
    :param fname: file path.
    :param strip_comments: Flag signaling whether to strip comments enclosed in square \
    brackets.
    :param kw: Keyword arguments are passed through to `Node.read`.
    :return: [`Shape`] instance.
    """
    l = newick.read(fname, encoding, strip_comments, **kw)
    return [newick_node_to_shape(n) for n in l]