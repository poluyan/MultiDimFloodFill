#ifndef TRIE_H
#define TRIE_H

#include <vector>
#include <algorithm>
#include <memory>

namespace trie_cpp
{
    template <template <typename> class T, typename I>
    struct TrieNode
    {
        I index;
        std::vector<std::unique_ptr<T<I>>> children;
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
        bool all_sorted;
    public:
        std::unique_ptr<T> root;
        Trie();
        ~Trie();
        bool empty() const;
        void insert(const std::vector<I> &key);
        //void binary_insert(const std::vector<I> &key);
        bool search(const std::vector<I> &key) const;
        bool binary_search(const std::vector<I> &key) const;
        size_t get_total_count() const;
        void fill_tree_count();
        void remove_tree();
        std::vector<I> get_and_remove_last();
        size_t get_number_of_nodes() const;
        void sort();
    private:
        void number_of_nodes(T *p, size_t &count) const;
        void get_number(T *p, size_t &count) const;
        void fill_tree_count(T *p);
        void remove_tree(T *p);
        bool is_all_end_of_word_true(T *p) const;
        void delete_last();
        void sort(T *p);
    };

    template <typename T, typename I>
    Trie<T,I>::Trie() : all_sorted(false), root(new T()) { }
    template <typename T, typename I>
    Trie<T,I>::~Trie()
    {
    }

    template <typename T, typename I>
    bool Trie<T,I>::empty() const
    {
        return root->children.empty();
    }

    template <typename T, typename I>
    void Trie<T,I>::insert(const std::vector<I> &key)
    {
        all_sorted = false;
        auto p = root.get();
        for(auto &i : key)
        {
            auto it = std::find_if(p->children.begin(), p->children.end(), [&i](const std::unique_ptr<T> &obj)
            {
                return obj->index == i;
            });
            if(it == p->children.end())
            {
                p->children.emplace_back(new T(i));
                p = p->children.back().get();
            }
            else
            {
                p = p->children[std::distance(p->children.begin(), it)].get();
            }
        }
        p->isEndOfWord = true;
    }

    template <typename T, typename I>
    bool Trie<T,I>::search(const std::vector<I> &key) const
    {
        auto p = root.get();
        for(auto &i : key)
        {
            auto it = std::find_if(p->children.begin(), p->children.end(), [&i](const std::unique_ptr<T> &obj)
            {
                return obj->index == i;
            });
            if(it == p->children.end())
            {
                return false;
            }
            else
            {
                p = p->children[std::distance(p->children.begin(), it)].get();
            }
        }
        return p->isEndOfWord;
        //return (p != 0 && p->isEndOfWord);
    }

    template <typename T, typename I>
    bool Trie<T,I>::binary_search(const std::vector<I> &key) const
    {
        auto p = root.get();
        std::unique_ptr<T> val(new T());
        for(auto &i : key)
        {
            val->index = i;
            auto it = std::lower_bound(p->children.begin(), p->children.end(), val, [](const std::unique_ptr<T> &l, const std::unique_ptr<T> &r)
            {
                return l->index < r->index;
            });
            if(it == p->children.end())
            {
                return false;
            }
            else
            {
                p = p->children[std::distance(p->children.begin(), it)].get();
            }
        }
        return p->isEndOfWord;
    }

    template <typename T, typename I>
    size_t Trie<T,I>::get_total_count() const
    {
        size_t count = 0;
        get_number(root.get(), count);
        return count;
    }

    template <typename T, typename I>
    void Trie<T,I>::fill_tree_count()
    {
        fill_tree_count(root);
    }

    template <typename T, typename I>
    void Trie<T,I>::remove_tree()
    {
        auto p = root.get();
        while(!is_all_end_of_word_true(p))
        {
            remove_tree(p);
        }
        p->children.clear();
    }

    template <typename T, typename I>
    void Trie<T,I>::sort()
    {
        sort(root.get());
        all_sorted = true;
    }

    template <typename T, typename I>
    void Trie<T,I>::sort(T *p)
    {
        std::sort(p->children.begin(),p->children.end(),[](const std::unique_ptr<T> &l, const std::unique_ptr<T> &r)
        {
            return l->index < r->index;
        });
        for(auto &i : p->children)
        {
            sort(i.get());
        }
    }

    template <typename T, typename I>
    std::vector<I> Trie<T,I>::get_and_remove_last()
    {
        std::vector<I> sample;
        auto p = root.get();
        if(p->children.empty())
            return sample;

        while(!p->children.back()->isEndOfWord)
        {
            sample.push_back(p->children.back()->index);
            p = p->children.back().get();
        }
        sample.push_back(p->children.back()->index);

        delete_last();

        return sample;
    }

    template <typename T, typename I>
    size_t Trie<T,I>::get_number_of_nodes() const
    {
        size_t count = 0;
        number_of_nodes(root.get(), count);
        return count;
    }

    template <typename T, typename I>
    void Trie<T,I>::number_of_nodes(T *p, size_t &count) const
    {
        count += p->children.size();
        for(auto &i : p->children)
        {
            number_of_nodes(i.get(), count);
        }
    }
    template <typename T, typename I>
    void Trie<T,I>::get_number(T* p, size_t &count) const
    {
        for(auto &i : p->children)
        {
            if(i->isEndOfWord)
                ++count;
            get_number(i.get(), count);
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
        for(const auto &i : p->children) // double ref
        {
            if(i->isEndOfWord)
            {
                is_end = true;
                break;
            }
            remove_tree(i.get());
        }
        if(is_end)
        {
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
        auto p = root.get();
        if(p->children.empty())
            return;

        while(!p->children.back()->isEndOfWord)
        {
            p = p->children.back().get();
        }
        p->children.pop_back();

        if(p->children.empty())
        {
            p->isEndOfWord = true;
            delete_last();
        }
    }
}

#endif
