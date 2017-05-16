"""
For us, a Phylogenetic Tree (a `PhyloTree` instance) is a special case of Shape in which labels of leaves can be distinguished.
We import modules Shape.py and newick.py; the latter will be used for reading Newick code in string format and turning it
into trees; it can be found in https://github.com/glottobank/python-newick.
"""
from Shape import *
import newick



class PhyloTree(Shape):
    """
     Thus, the class `PhyloTree` must extend the class `Shape`, overriding those methods in `Shape` that do not take into
     account the distinguishability of leaves.
    """

    def __init__(self, leaf, children):

        """
        Create a new `PhyloTree` object. 
        The boolean is_leaf is True if the object is a leaf; it is False otherwise.
        :param leaf: if is_leaf, its label; otherwise, `None`.
        :param children: if not is_leaf, the `PhyloTree` objects which are descendants of the considered object; 
        otherwise, `None`.
        :return: `PhyloTree` instance.
        """
        super(PhyloTree, self).__init__(children)
        self.leaf = leaf
        assert not self.is_leaf or (leaf is not None)


    def compare(self, T2):
        """
        Compare self with another `PhyloTree` object. We use lexicographical order in order to compare two `PhyloTree` 
        instances. This method overrides that of `Shape`, since any two `PhyloTree` leaves are comparable. It returns an
        int c, which is 0 if self and T2 are equal, < 0 if self < T2, and > 0 if self > T2.
        :param T2: the `PhyloTree` object against which we compare self.
        :return: int instance.
        """
        if self.is_leaf and T2.is_leaf:
            if self.leaf   > T2.leaf:
                return 1
            elif self.leaf < T2.leaf:
                return -1
            else:
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


    def shape(self):
        """
        Returns the `Shape` associated to self. Namely, it "forgets" the labels of the leafs.
        :return: `Shape` instance.
        """
        if self.is_leaf:
            return Shape(None)
        else:
            return Shape([x.shape() for x in self.children])


    def to_newick_tuple(self):
        """
        Returns a tuple representing the simplified Newick code of self.
        :return: tuple instance.
        """
        if self.is_leaf:
            return self.leaf
        else:
            return tuple(x.to_newick_tuple() for x in self.children)

    def leaves(self):
        """
        Yields the (labels of the) leaves of self.
        :return: `PhyloTree` instance.
        """
        if self.is_leaf:
            yield self.leaf
        else:
            for x in self.children:
                for l in x.leaves():
                    yield l


    def labels(self):
        """
        Returns a list with the labels that appear in self, sorted in lexicographical order.
        Repetitions may arise if the user enters trees which are not phylogenetic.
        :return: list instance.
        """
        return sorted(list(self.leaves()))


    def is_phylo(self):
        """
        Returns True if self is phylogenetic (namely, if it has no repeated leaves). Returns False otherwise.
        :return: bool instance.
        """
        L = self.labels()
        for i in range(1, len(L)):
            if L[i] == L[i-1]:
                return False
        return True


def from_newick(X):
    """
    Create a `PhyloTree` object from a Newick code entered as a string.
    :param X: a string representing a Newick code.
    :return: `PhyloTree` instance.
    """
    return newick_node_to_tree(newick.loads(X)[0])

def from_newick_list(X):
    """
    Create a list of `PhyloTree` objects from a list of Newick codes entered as a string.
    :param X: a string representing a list of Newick codes.
    :return: [`PhyloTree`] instance.
    """
    return[newick_node_to_tree(n) for n in newick.loads(X)]

def newick_node_to_tree(N):
    """
    Create a `PhyloTree` object from a `Node` object.
    :param N: a Node.
    :return: `PhyloTree` instance.
    """
    if not bool(N.descendants):
        return PhyloTree(N.name, None)
    return PhyloTree(None, sorted([newick_node_to_tree(x) for x in N.descendants]))

def trees_from_file(fname, encoding='utf8', strip_comments=False, **kw):
    """
    Load a list of trees from a Newick formatted file.
    :param fname: file path.
    :param strip_comments: Flag signaling whether to strip comments enclosed in square \
    brackets.
    :param kw: Keyword arguments are passed through to `Node.read`.
    :return: [`PhyloTree`] instance.
    """
    l = newick.read(fname, encoding, strip_comments, **kw)
    return [newick_node_to_tree(n) for n in l]