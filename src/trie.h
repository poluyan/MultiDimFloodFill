#ifndef TRIE_H_INCLUDED
#define TRIE_H_INCLUDED

#include <vector>

namespace trie_c
{

struct TrieNode
{
    TrieNode **children;
    bool isEndOfWord;
};

TrieNode *getNode(int asize);
void insert(TrieNode *root, std::vector<int> key, int asize);
bool search(TrieNode *root, std::vector<int> key);

}

namespace trie_cpp
{

template <template <typename> class T, typename I>
class TrieNode
{
public:
    I index;
    std::vector<T<I>*> children;
    bool isEndOfWord;
    TrieNode() : index(0), isEndOfWord(0) { }
    TrieNode(I ind) : index(ind), isEndOfWord(0) { }
};

template <typename I>
class Node: public TrieNode<Node, I>
{
public:
    Node() : TrieNode<Node, I>() {}
    Node(I ind) : TrieNode<Node, I>(ind) {}
};

template <typename I>
class NodeCount: public TrieNode<NodeCount, I>
{
public:
    size_t count;
    NodeCount() : TrieNode<NodeCount, I>(), count(0) {}
    NodeCount(I ind) : TrieNode<NodeCount, I>(ind), count(0) {}
};

template <typename T, typename I>
class Trie
{
public:
    T *root;

    Trie()
    {
        root = new T();
    }
    ~Trie()
    {
        delete root;
    }
};

}

#endif
