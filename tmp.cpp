#include<iostream>
#include<algorithm>
#include<vector>
#include<string>
#include<map>
#include<set>
using namespace std;

class Solution {
public:
    map<int, int> table;
    vector<vector<int>> res;
    map<int, int> upper;
    
    void dfs(vector<int>& can, int target, int s) {    
        if(target == 0) {
            vector<int> r;
            for(auto& c : can) {
                if(table.count(c) > 0) {
                    for(int i = 0; i < table[c]; i++) {
                        r.push_back(c);
                    }
                }
            }
            res.push_back(r);
            return;
        }
        if(s == can.size())return;
        
        table[can[s]] = 0;
        do {
            dfs(can, target - table[can[s]] * can[s], s + 1);
            table[can[s]]++;
        } while(table[can[s]] * can[s] <= target && table[can[s]] <= upper[can[s]]);
        // cout << can[s] << ":" << table[can[s]] << endl;
        table[can[s]] = 0;
        return;
    }
    
    vector<vector<int>> combinationSum2(vector<int>& candidates, int target) {
        sort(candidates.begin(), candidates.end());
        for(auto& c : candidates) {
            if(upper.count(c) == 0)upper[c] = 0;
            upper[c]++;
        }
        // set<int> s(candidates.begin(), candidates.end());
        // candidates = vector<int>(s.begin(), s.end());
        
        dfs(candidates, target, 0);
        return res;
    }
};

int main() {
    Solution s;
    vector<int> v({2, 2});
    auto res = s.combinationSum2(v, 4);
    for(auto& r1 : res) {
        for(auto& r2 : r1) {
            cout << r2 << ",";
        }
        cout << endl;
    }
    return 0;
}
