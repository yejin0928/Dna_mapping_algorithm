//2019112007_�ǿ���
#define _CRT_SECURE_NO_WARNINGS

#define MAX 20000//��ü reference ����
#define N_SIZE 1000//short read ����
#define K_SIZE 70//short read ����
//needleman_wunsch�˰��� ����ϴ� ���� ü��
#define MATCH 1//��ġ
#define MISMATCH - 1//����ġ
#define GAP -2//��

#include <iostream>
#include <iomanip>
#include <utility>
#include<stdio.h>
#include <string>
#include <vector>
#include <map>
#include <math.h>
#include <sstream>
#include <algorithm>
#include<cstdlib>
#include<time.h>
#include<fstream>
#include <string.h>
#include<random>
#include <queue>
#include <list>
#include <cstring>
#include <cstdio>
using namespace std;

char* restruct_m_DNA = new char[MAX + 1];//������ my dna ���� ������ ����
void make_rDNA(string& r);//reference dna�� ����� �Լ�
void make_my_DNA(string& m, string& r);//my dna�� ����� �Լ�
void make_myDNA_s(string& m, string& s, int k, int n);//short read�� ����� �Լ�

void make_rDNA(string& r) {
    /*�޸��� Ʈ������ ���� �����⸦ �̿��� 0-3�������� �յ��ϰ� ��еǾ� 
    �߻��� ���ڿ� 0�϶�, A, 1�϶�, C, 2�϶� G, 3�϶� T���� �� ���ڿ� ���ڸ� �����Ͽ�
    DNA�� ������. */
    random_device rn;
    mt19937_64 rand(rn());
    uniform_int_distribution<int> range(0, 3);
    //srand(6783);
    int cnt = MAX;
    char* r_DNA = new char[MAX + 1];

    char num_array[5] = "ACGT";

    for (int i = 0; i < MAX; i++) {
        int random = range(rand);
        r_DNA[i] = num_array[random];
    }
    r_DNA[MAX] = '\0';

    string str(r_DNA);
    r = str;
    ofstream fout("Reference_DNA.txt");
    fout << r_DNA;
    fout.close();
    delete r_DNA;

}

void make_my_DNA(string& m, string& r) {
    /*MY DNA�� �ռ� ������ reference DNA���� ���̸� 50%�� ������, 
    for���ȿ��� MAX �� ���� �ݺ����� �����ϴµ�, 
    (reference DNA �ε����� ����+1)�� 2�� ����� reference DNA������ ����,
    while���ȿ��� ���� referenceDNA �ε����� ���ڿ� �ٸ� ���� ���� ������ 
    �޸��� Ʈ������ ���� �����⸦ �̿��� �߻��� 0-3������ ������ ���� 
    0�϶�, A, 1�϶�, C, 2�϶� G, 3�϶� T ���ڷ� �ٲ㼭 reference DNA�� �ٸ��� �����. */
    random_device rn;
    mt19937_64 rand(rn());
    uniform_int_distribution<int> range(0, 3);
    char* my_DNA = new char[MAX + 1];
    int diff = 50;
    for (int i = 0; i < MAX; i++) {
        my_DNA[i] = r[i];
        if ((i != 0) && ((i + 1) % diff == 0))
        {
            while (my_DNA[i] == r[i]) {
                int dif = range(rand);
                switch (dif) {
                case 0:
                    my_DNA[i] = 'A';
                    break;
                case 1:
                    my_DNA[i] = 'T';
                    break;
                case 2:
                    my_DNA[i] = 'G';
                    break;
                case 3:
                    my_DNA[i] = 'C';
                    break;

                }
            }

        }
    }
    my_DNA[MAX] = '\0';
    string str(my_DNA);
    m = str;

    ofstream fo("my_DNA.txt");
    fo << my_DNA;
    fo.close();
    delete my_DNA;
}

void make_myDNA_s(string& m, string& s, int k, int n) {
    /*Short read�� �޸��� Ʈ������ ���� �����⸦ �̿��� ���� 0-MAX-k-1 �������� �յ��ϰ� ��еǾ� 
    �߻��� ���� ���� short read�� �ڸ� ��ġ�� �����Ͽ� for���ȿ��� �ݺ� Ƚ�� N�� ���� �ش� ��ġ���� 
    ���̰� k�� short read�� �����Ѵ�. */
    string buffer = "";

    ofstream o("mydna_short.txt");
    random_device rn1;
    mt19937_64 rand(rn1());
    uniform_int_distribution<int> range1(0, MAX - k - 1);

    for (int i = 0; i < n; i++) {
        int i_rand = range1(rand);
        buffer = "";
        for (int j = i_rand; j < i_rand + k; j++) {
            buffer += m[j];
            //s += buffer;
            //s += "\0";
        }
        o << buffer << endl;
    }
    o.close();


}

class score {//needleman wunsch�˰��򿡼� ���̴� Ŭ����
public:
    int num;//����
    string prev_dir;//����
    score() :num(0), prev_dir("None") {}
};
map <string, unsigned short> string_to_short;



void result(vector <string> seq, string w, int f_position, int l_position, int m_index, int k, int n) {
    /*��� ������ ���� ������ short read�� sequence�� ���� ������ for���ȿ��� üũ�ؼ� ó������ ACGT �� 
    �ϳ��� ���� ��Ÿ���� ���� �ε����� ���������� ��Ÿ���� ���� �ε����� ���� first�� last������ �����Ѵ�. 
    �������� ��� sequence���� string ������ w�� �����ؼ�, while���ȿ��� first���� last������ �ε��� ����, 
    ����ġ�� ���ϴ� -�� �ƴ� ��(ACGT �� �ϳ��� ��)�� �߰ߵǾ��� �ÿ� �ش� �ε������� short read�� ���� K���� 
    while���ȿ��� ����ġ�� ���ϴ� -�� �ƴ� ��(ACGT�� �ϳ��� ��)�� ��Ÿ���� Ƚ���� �����Ѵ�. 
    ��Ÿ�� Ƚ���� max_count���� ũ�� max_count���� �ش� Ƚ���� ��ü�ϰ� �ε��� ���� max_index�� �����Ѵ�. 
    �� ������ while ���ȿ��� �ݺ��ϸ�, short�� reference DNA�� ���� ���� �κ� ���ڿ��� ��ġ�ϴ� ��ġ�� ã�� �� �ִ�. 
    ��, max_index�� reference DNA�� �κ� ���ڿ��� short read�� ��ü�� ���� �ε����̴�. 
    ���� ���������� max_index���� ���� K���� reference DNA�� ���ڿ��� short read�� ��ü�Ѵ�.
    */
    //cout << "result" << endl;
    map <string, unsigned short>::iterator iter;
    int count = 0;
    int first = 0;
    int last = 0;
    int len = 0;
    try {
        for (int i = 0; i < seq.size(); i++) {
            for (iter = string_to_short.begin(); iter != string_to_short.end(); iter++) {
                if (to_string(iter->second) == seq[i]) {
                    if (count == 0) {
                        first = i;
                    }
                    else if (count == (k - 1))
                        last = i;
                 
                    count++;
                  
                    w += iter->first;
                }
            }
            if (seq[i] == "___") {
               
                w += "-";
            }
        }

        f_position = first;
        l_position = last;
        w += "\0";

        int index = f_position;
        int max_count = 0;
        int max_index = 0;
        while (index != l_position) {
            if (w[index] != '-') {
                int count = 0;
                int dna_c = 0;
                while ((count < k) && (index + count) <= l_position) {
                    if (w[index + count] != '-')
                        dna_c++;
                    count++;
                }
                if (dna_c > max_count) {
                    max_count = dna_c;
                    max_index = index;
                }
                index++;
            }
            else
                index++;
        }
        m_index = max_index;

        if ((m_index + k) <= MAX) {
            for (int yes = m_index; yes < m_index + k; yes++)
                if(w[yes]!='-')
                 restruct_m_DNA[yes] = w[yes];
        }
    }
    catch (bad_alloc ex) {
        cout << ex.what() << endl;

    }

}


void back_track(vector<vector<score>> match_matrix, vector <unsigned short> a, vector <unsigned short> b, int k, int n) {
    /*�ڿ������� �Ž��� ���鼭 matrix�� ����� ������ ����ؼ� �� ���ڿ� ������ sequence���� string���� ��ȯ�ؼ� �����Ѵ�.
    ��ġ �� �ش� ���� ����, ����ġ�� ��-���� �����Ѵ�. �׸��� ���� ��� ������ short read�� sequence�� �Ű������� ������ 
    �Լ� result�� ȣ���Ѵ�.*/
   // cout << "back" << endl;
    vector <string> seq_a;
    vector <string> seq_b;

    try {
        int match_index = max(a.size(), b.size()) - 1;
        int m_i = match_index + 1;
        seq_a.resize(m_i);
        seq_b.resize(m_i);
        int i = a.size();
        int j = b.size();
        for (; match_index >= 0; match_index--) {
            if (match_matrix[i][j].prev_dir == "Left") {
                seq_a[match_index] = "___";
                seq_b[match_index] = to_string(b[j - 1]);
                j--;
            }
            else if (match_matrix[i][j].prev_dir == "Up") {
                seq_a[match_index] = to_string(a[i - 1]);
                seq_b[match_index] = "___";
                i--;
            }
            else if (match_matrix[i][j].prev_dir == "Left Up") {
                seq_a[match_index] = to_string(a[i - 1]);
                seq_b[match_index] = to_string(b[j - 1]);
                i--;
                j--;
            }
            else {
                if (match_matrix[i][j].prev_dir == "None" && i == 0 && j == 0) {
                    break;
                }
            }
        }

        string w_string = "";
        int l_p = 0, f_p = 0;

        int mat_index = 0;
        result(seq_b, w_string, f_p, l_p, mat_index,k,n);
    }
    catch (bad_alloc ex) {
        cout << ex.what() << endl;
    }
}


void forward(vector<vector<score>>& match_matrix, vector <unsigned short> a, vector <unsigned short> b, int k, int n) {
    /*�տ������� ����� �����Ѵ�. 
    ����� ���� ���ڷ� �����ϴµ�, ���� �ο� �ý����� �� ���ڿ��� ���� ��ġ �� 1, ����ġ�� -1, �� -2�̴�. 
    �������������� ����ؼ� gap�� ���� ������ ������ ������ �߰��Ѵ�. �׸��� �ְ����� �ִ� ������ ���Ѵ�.*/
   // cout << "forward" << endl;
    int match_score = MATCH;
    int mismatch_score = MISMATCH;
    int gap_score = GAP;
    int as = a.size() + 1;
    int bs = b.size() + 1;
    try {
        for (int i = 0; i < as; i++) {
            for (int j = 0; j < bs; j++) {
                if (i == 0 && j == 0) {
                    match_matrix[i][j].num = 0;
                }
                else if (i == 0) {
                    match_matrix[i][j].num = match_matrix[i][j - 1].num + gap_score;
                    match_matrix[i][j].prev_dir = "Left";
                }
                else if (j == 0) {
                    match_matrix[i][j].num = match_matrix[i - 1][j].num + gap_score;
                    match_matrix[i][j].prev_dir = "Up";
                }
                else {
                    int from_left_score = match_matrix[i][j - 1].num + gap_score;
                    int from_up_score = match_matrix[i - 1][j].num + gap_score;
                    int from_left_up_score = 0;
                    if (a[i - 1] == b[j - 1]) {
                        //��ġ
                        from_left_up_score = match_matrix[i - 1][j - 1].num + match_score;
                    }
                    else {
                        //����ġ
                        from_left_up_score = match_matrix[i - 1][j - 1].num + mismatch_score;
                    }

                    //�ִ밪 ���� ��
                    if (from_left_score > from_up_score) {
                        if (from_left_score >= from_left_up_score) {
                            match_matrix[i][j].num = from_left_score;
                            match_matrix[i][j].prev_dir = "Left";
                        }
                        else {
                            match_matrix[i][j].num = from_left_up_score;
                            match_matrix[i][j].prev_dir = "Left Up";
                        }
                    }
                    else {
                        if (from_up_score >= from_left_up_score) {
                            match_matrix[i][j].num = from_up_score;
                            match_matrix[i][j].prev_dir = "Up";
                        }
                        else {
                            match_matrix[i][j].num = from_left_up_score;
                            match_matrix[i][j].prev_dir = "Left Up";
                        }
                    }
                }
            }
        }
    }
    catch (bad_alloc ex) {
        cout << ex.what() << endl;
    }
}

void needleman_wunsch(vector <unsigned short> a, vector <unsigned short> b, int k, int n) {
    /*��:reference����+1, ��:short read�� ����+1�� matrix�� �����Ѵ�. 
    �׸��� forward�Լ��� back_track�Լ��� ȣ���Ѵ�.*/
   // cout << "needlman" << endl;
    vector<vector<score>> match_matrix;
    int as = a.size() + 1;
    int bs = b.size() + 1;
    try {
        match_matrix.resize(as);
        for (int i = 0; i < as; i++) {
            match_matrix[i].resize(bs);
        }

        forward(match_matrix, a, b,k,n);
        back_track(match_matrix, a, b, k,n);
    }
    catch (bad_alloc ex) {
        cout << ex.what() << endl;

    }
}


void convert(string a, string b, int k, int n) {
    /*�Է� ���� �� ���ڿ� reference DNA�� short read�� string_to_short�� ��ȯ�ϴ� ������ ��ģ��.
    ���⼭ string_to_short�� map()�� ����Ͽ� ���ڿ��� ������ ��ȯ�ϴ� ������ ���Ѵ�.
    �׸��� ���������� needleman_wunsch �Լ��� ȣ���Ѵ�.*/
   // cout << "convert" << endl;

    vector <unsigned short> seq_a;
    vector <unsigned short> seq_b;

    for (int i = 0; i < a.length(); i++) {
        stringstream ss;
        string s;
        ss << a[i];
        ss >> s;
        seq_a.push_back(string_to_short[s]);
    }
    for (int i = 0; i < b.length(); i++) {
        stringstream ss;
        string s;
        ss << b[i];
        ss >> s;
        seq_b.push_back(string_to_short[s]);
    }
    needleman_wunsch(seq_a, seq_b, k, n);
}

void input(string& r, int k, int n) {
    /*Needleman-Wunsch �˰����� �̿��� ���� �� ���ڿ��� �����Ѵ�. 
    Short read�� ����ִ� ������ ���� �� �� �� ����Ǿ� �ִ� short read�� 
    string������buffer�� �����ϰ� ���������� ������ �Լ� convert()�� ȣ���Ѵ�. 
    �׸��� buffer�� �� ���ڿ��� �ʱ�ȭ�Ѵ�. �� ������ while���ȿ��� ������ ������ �� �� �� ��� ���� ������, 
    ��, short read�� ������ N����ŭ �ݺ��Ѵ�.   */
    //cout << "input" << endl;
    clock_t start_b = clock();
    for (int i = 0; i < MAX; i++)
        restruct_m_DNA[i] = r[i];
    restruct_m_DNA[MAX] = '\0';
    ifstream f("mydna_short.txt");
    string buffer = "";


    double total = 0;
    while (getline(f, buffer)) {

        string a, b;
        a = r;
        b = buffer;

        convert(a, b,k,n);
        buffer = "";
    }

    f.close();
    ofstream foo("restruct_m_DNA.txt");
    foo << restruct_m_DNA;
    foo.close();
    string stra(restruct_m_DNA);
    r = stra;

    clock_t end_b = clock();
    cout << "�ɸ� �ð�:" << total + ((double)(end_b - start_b) / CLOCKS_PER_SEC) << "��" << endl;
    delete restruct_m_DNA;

}


int main(void)
{
    int n, k;
    string r_DNA = "";
    string m_DNA = "";
    string m_DNA_s = "";
    string restruct_m_DNA;

    make_rDNA(r_DNA);
    r_DNA += '\0';
    restruct_m_DNA = r_DNA;

    make_my_DNA(m_DNA, r_DNA);
    double dif = 0;
    for (int i = 0; i < MAX; i++)
        if (r_DNA[i] != m_DNA[i])
            dif++;
    cout << "������: " << 100 - (dif / MAX) * 100 << "%" << endl;//reference DNA�� my dna�� ����
    cout << "���� k�� �Է��ϼ���:";
    cin >> k;
    cout << "���� n�� �Է��ϼ���:";
    cin >> n;
    //k = K_SIZE;
   // n = N_SIZE;
    make_myDNA_s(m_DNA, m_DNA_s, k, n);
    m_DNA_s += '\0';

    //ACGT DNA ����� �������� ���ڰ��� ��������
    string_to_short.insert(pair<string, unsigned short>("G", 15));
    string_to_short.insert(pair<string, unsigned short>("C", 256));
    string_to_short.insert(pair<string, unsigned short>("A", 326));
    string_to_short.insert(pair<string, unsigned short>("T", 789));

    input(restruct_m_DNA, k, n);

    double diff = 0;
    for (int i = 0; i < MAX; i++)
        if (restruct_m_DNA[i] != m_DNA[i])
            diff++;
    cout << "������: " << 100 - (diff / MAX) * 100 << "%" << endl;//������ MY DNA�� ���� MY DNA�� ����

    system("pause");
}

