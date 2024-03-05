#include <iostream>
#include <vector>


class Solution {
public:
    static int compress(std::vector<char>& chars) {
        std::vector <char> s;
        for (int i = 0; i < chars.size(); i++) {
            for (int j = i; j <= chars.size(); j++) {
                if (j == chars.size() || chars[j] != chars[i]) {
                    if (j - i == 1) {
                        s.push_back(chars[i]);
                    }
                    else {
                        s.push_back(chars[i]);
                        std::vector <char> num;
                        int col = j - i;
                        while(col) {
                            num.push_back((col % 10) + '0');
                            col /= 10;
                        }
                        for (int k = num.size() - 1; k >= 0; k--) {
                            s.push_back(num[k]);
                        }
                    }
                    i = j - 1;
                    break;
                }
            }
        }
        chars = s;
        return chars.size();
    }
};

int main() {
    std::vector <char> s = {'a', 'a', 'b', 'b', 'b', 'b', 'b', 'b', 'b', 'b', 'b', 'b'};
    std::cout << Solution::compress(s) << std::endl;
    for (auto letter : s) {
        std::cout << letter << " ";
    }
}