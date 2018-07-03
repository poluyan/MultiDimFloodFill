#include "trie.h"

namespace trie
{
   
TrieNode *getNode(int asize)
{
    TrieNode *p =  new TrieNode;
    p->isEndOfWord = false;
    p->children = new TrieNode* [asize];
    for(int i = 0; i < asize; i++)
        p->children[i] = 0;
    return p;
}

void insert(TrieNode *root, std::vector<int> key, int asize)
{
    TrieNode *p = root;
    for(size_t i = 0; i < key.size(); i++)
    {
        int index = key[i];
        if(!p->children[index])
            p->children[index] = getNode(asize);
        p = p->children[index];
    }
    p->isEndOfWord = true;
}

bool search(TrieNode *root, std::vector<int> key)
{
    TrieNode *p = root;
    for(size_t i = 0; i < key.size(); i++)
    {
        int index = key[i];
        if(!p->children[index])
            return false;
        p = p->children[index];
    }
    return (p != 0 && p->isEndOfWord);
}

}