#ifndef TRIE_H
#define TRIE_H

#include <vector>
#include <algorithm>

namespace trie_cpp
{
    template <template <typename> class T, typename I>
    struct TrieNode
    {
        I index;
        std::vector<T<I>*> children;
        bool isEndOfWord;
        TrieNode() : index(0), isEndOfWord(0) { }
        TrieNode(I ind) : index(ind), isEndOfWord(0) { }
    };

    template <typename I>
    struct Node: public TrieNode<Node, I>
    {
        Node() : TrieNode<Node, I>() {}
        Node(I ind) : TrieNode<Node, I>(ind) {}
    };

    template <typename I>
    struct NodeCount: public TrieNode<NodeCount, I>
    {
        size_t count;
        NodeCount() : TrieNode<NodeCount, I>(), count(0) {}
        NodeCount(I ind) : TrieNode<NodeCount, I>(ind), count(0) {}
    };

    template <typename T, typename I>
    class Trie
    {
    public:
        T *root;
        Trie() noexcept;
        ~Trie();
        bool empty() const;
        void insert(std::vector<I> key);
        bool search(std::vector<I> key) const;
        size_t get_total_count() const;
        void fill_tree_count();
        void remove_tree();
        std::vector<I> get_and_remove_last();
        size_t get_number_of_nodes() const;
    private:
        void number_of_nodes(T *p, size_t &count) const;
        void get_number(T *p, size_t &count) const;
        void fill_tree_count(T *p);
        void remove_tree(T *p);
        bool is_all_end_of_word_true(T *p) const;
        void delete_last();
    };

    template <typename T, typename I>
    Trie<T,I>::Trie() noexcept
    {
        root = new T();
    }
    template <typename T, typename I>
    Trie<T,I>::~Trie()
    {
        delete root;
    }

    template <typename T, typename I>
    bool Trie<T,I>::empty() const
    {
        return root->children.empty();
    }

    template <typename T, typename I>
    void Trie<T,I>::insert(std::vector<I> key)
    {
        T *p = root;
        for(auto &i : key)
        {
            auto it = std::find_if(p->children.begin(), p->children.end(), [&i](const T* obj)
            {
                return obj->index == i;
            });
            if(it == p->children.end())
            {
                T *temp = new T(i);
                p->children.emplace_back(std::move(temp));
                p = p->children.back();
            }
            else
            {
                p = p->children[std::distance(p->children.begin(), it)];
            }
        }
        p->isEndOfWord = true;
    }

    template <typename T, typename I>
    bool Trie<T,I>::search(std::vector<I> key) const
    {
        T *p = root;
        for(auto &i : key)
        {
            auto it = std::find_if(p->children.begin(), p->children.end(), [&i](const T* obj)
            {
                return obj->index == i;
            });
            if(it == p->children.end())
            {
                return false;
            }
            else
            {
                p = p->children[std::distance(p->children.begin(), it)];
            }
        }
        return p->isEndOfWord;
    }

    template <typename T, typename I>
    size_t Trie<T,I>::get_total_count() const
    {
        size_t count = 0;
        get_number(root, count);
        return count;
    }

    template <typename T, typename I>
    void Trie<T,I>::fill_tree_count()
    {
        fill_tree_count(root);
        //root->count = get_total_count();
    }

    template <typename T, typename I>
    void Trie<T,I>::remove_tree()
    {
        T *p = root;
        while(!is_all_end_of_word_true(p))
        {
            remove_tree(p);
        }
        for(auto &i : root->children)
            delete i;
        p->children.clear();
    }

    template <typename T, typename I>
    std::vector<I> Trie<T,I>::get_and_remove_last()
    {
        std::vector<I> sample;
        T *p = root;
        if(p->children.empty())
            return sample;

        while(!p->children.back()->isEndOfWord)
        {
            sample.push_back(p->children.back()->index);
            p = p->children.back();
        }
        sample.push_back(p->children.back()->index);

        delete_last();

        return sample;
    }

    template <typename T, typename I>
    size_t Trie<T,I>::get_number_of_nodes() const
    {
        size_t count = 0;
        number_of_nodes(root, count);
        return count;
    }

    template <typename T, typename I>
    void Trie<T,I>::number_of_nodes(T *p, size_t &count) const
    {
        count += p->children.size();
        for(auto &i : p->children)
        {
            number_of_nodes(i, count);
        }
    }
    template <typename T, typename I>
    void Trie<T,I>::get_number(T *p, size_t &count) const
    {
        for(auto &i : p->children)
        {
            if(i->isEndOfWord)
                ++count;
            get_number(i, count);
        }
    }
    template <typename T, typename I>
    void Trie<T,I>::fill_tree_count(T *p)
    {
        for(auto &i : p->children)
        {
            size_t count = 0;
            get_number(i, count);
            i->count = count > 0 ? count : 1;
            fill_tree_count(i);
        }
    }
    template <typename T, typename I>
    void Trie<T,I>::remove_tree(T *p)
    {
        bool is_end = false;
        for(auto &i : p->children)
        {
            if(i->isEndOfWord)
            {
                is_end = true;
                break;
            }
            remove_tree(i);
        }
        if(is_end)
        {
            for(auto &i : p->children)
                delete i;
            p->children.clear();
            p->isEndOfWord = true;
        }
    }
    template <typename T, typename I>
    bool Trie<T,I>::is_all_end_of_word_true(T *p) const
    {
        for(auto &i : p->children)
        {
            if(!i->isEndOfWord)
                return false;
        }
        return true;
    }
    template <typename T, typename I>
    void Trie<T,I>::delete_last()
    {
        T *p = root;
        if(p->children.empty())
            return;

        while(!p->children.back()->isEndOfWord)
        {
            p = p->children.back();
        }

        delete p->children.back();
        p->children.pop_back();

        if(p->children.empty())
        {
            p->isEndOfWord = true;
            delete_last();
        }
    }
}

#endif
