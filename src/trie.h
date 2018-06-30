#ifndef TRIE_H_INCLUDED
#define TRIE_H_INCLUDED

#include <vector>

namespace trie
{

struct TrieNode
{
    TrieNode **children;
    bool isEndOfWord;
};

TrieNode *getNode(int asize);
void insert(TrieNode *root, std::vector<int> key, int asize);
bool search(TrieNode *root, std::vector<int> key);

class Trie
{
private:
    TrieNode *root;
    Trie()
    {
        root = new TrieNode();
    }
};

}

#endif
