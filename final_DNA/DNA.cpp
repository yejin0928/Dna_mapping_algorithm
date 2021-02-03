//2019112007_권예진
#define _CRT_SECURE_NO_WARNINGS

#define MAX 20000//전체 reference 길이
#define N_SIZE 1000//short read 개수
#define K_SIZE 70//short read 길이
//needleman_wunsch알고리즘에 사용하는 점수 체계
#define MATCH 1//일치
#define MISMATCH - 1//불일치
#define GAP -2//갭

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

char* restruct_m_DNA = new char[MAX + 1];//복원한 my dna 값을 저장할 변수
void make_rDNA(string& r);//reference dna를 만드는 함수
void make_my_DNA(string& m, string& r);//my dna를 만드는 함수
void make_myDNA_s(string& m, string& s, int k, int n);//short read를 만드는 함수

void make_rDNA(string& r) {
    /*메르센 트위스터 난수 생성기를 이용해 0-3구간으로 균등하게 배분되어 
    발생된 숫자에 0일때, A, 1일때, C, 2일때 G, 3일때 T으로 각 문자에 숫자를 배정하여
    DNA를 만들어낸다. */
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
    /*MY DNA는 앞서 생성된 reference DNA와의 차이를 50%로 설정해, 
    for문안에서 MAX 번 동안 반복문을 수행하는데, 
    (reference DNA 인덱스의 숫자+1)이 2의 배수인 reference DNA문자의 값을,
    while문안에서 현재 referenceDNA 인덱스의 문자와 다른 값을 가질 때까지 
    메르센 트위스터 난수 생성기를 이용해 발생한 0-3구간의 난수에 따라 
    0일때, A, 1일때, C, 2일때 G, 3일때 T 문자로 바꿔서 reference DNA와 다르게 만든다. */
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
    /*Short read는 메르센 트위스터 난수 생성기를 이용해 숫자 0-MAX-k-1 구간에서 균등하게 배분되어 
    발생된 숫자 값을 short read를 자를 위치로 지정하여 for문안에서 반복 횟수 N번 동안 해당 위치부터 
    길이가 k인 short read를 생성한다. */
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

class score {//needleman wunsch알고리즘에서 쓰이는 클래스
public:
    int num;//점수
    string prev_dir;//방향
    score() :num(0), prev_dir("None") {}
};
map <string, unsigned short> string_to_short;



void result(vector <string> seq, string w, int f_position, int l_position, int m_index, int k, int n) {
    /*모든 과정을 거쳐 생성된 short read의 sequence가 끝날 때까지 for문안에서 체크해서 처음으로 ACGT 중 
    하나의 값이 나타났을 때의 인덱스와 마지막으로 나타났을 때의 인덱스를 각각 first와 last변수에 저장한다. 
    다음으로 모든 sequence값을 string 변수인 w에 저장해서, while문안에서 first부터 last사이의 인덱스 동안, 
    비일치를 뜻하는 -가 아닌 값(ACGT 중 하나의 값)이 발견되었을 시에 해당 인덱스부터 short read의 길이 K까지 
    while문안에서 비일치를 뜻하는 -가 아닌 값(ACGT중 하나의 값)이 나타나는 횟수를 저장한다. 
    나타난 횟수가 max_count보다 크면 max_count값을 해당 횟수로 교체하고 인덱스 값은 max_index에 저장한다. 
    이 과정을 while 문안에서 반복하면, short가 reference DNA와 가장 많은 부분 문자열이 일치하는 위치를 찾을 수 있다. 
    즉, max_index는 reference DNA의 부분 문자열을 short read로 교체할 시작 인덱스이다. 
    따라서 마지막으로 max_index부터 길이 K까지 reference DNA의 문자열을 short read로 교체한다.
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
    /*뒤에서부터 거슬러 가면서 matrix에 저장된 방향을 계산해서 두 문자열 각각의 sequence값을 string으로 변환해서 저장한다.
    일치 시 해당 문자 값을, 비일치시 “-“를 저장한다. 그리고 계산된 결과 생성된 short read의 sequence를 매개변수로 전달한 
    함수 result를 호출한다.*/
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
    /*앞에서부터 행렬을 구성한다. 
    행렬의 값은 숫자로 저장하는데, 숫자 부여 시스템은 두 문자열의 값이 일치 시 1, 불일치시 -1, 갭 -2이다. 
    시작점에서부터 출발해서 gap을 각각 끝까지 도달할 때까지 추가한다. 그리고 최고점이 있는 방향을 비교한다.*/
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
                        //일치
                        from_left_up_score = match_matrix[i - 1][j - 1].num + match_score;
                    }
                    else {
                        //불일치
                        from_left_up_score = match_matrix[i - 1][j - 1].num + mismatch_score;
                    }

                    //최대값 방향 비교
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
    /*행:reference길이+1, 열:short read의 길이+1로 matrix를 생성한다. 
    그리고 forward함수와 back_track함수를 호출한다.*/
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
    /*입력 받은 두 문자열 reference DNA와 short read를 string_to_short로 변환하는 과정을 거친다.
    여기서 string_to_short란 map()을 사용하여 문자열을 정수로 변환하는 과정을 말한다.
    그리고 마지막으로 needleman_wunsch 함수를 호출한다.*/
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
    /*Needleman-Wunsch 알고리즘을 이용해 비교할 두 문자열을 전달한다. 
    Short read가 담겨있는 파일을 열어 한 줄 씩 저장되어 있는 short read를 
    string변수인buffer에 저장하고 다음과정을 수행할 함수 convert()를 호출한다. 
    그리고 buffer를 빈 문자열로 초기화한다. 이 과정을 while문안에서 파일의 내용을 한 줄 씩 모두 읽을 때까지, 
    즉, short read의 개수인 N번만큼 반복한다.   */
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
    cout << "걸린 시간:" << total + ((double)(end_b - start_b) / CLOCKS_PER_SEC) << "초" << endl;
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
    cout << "복원률: " << 100 - (dif / MAX) * 100 << "%" << endl;//reference DNA와 my dna의 차이
    cout << "길이 k를 입력하세요:";
    cin >> k;
    cout << "개수 n을 입력하세요:";
    cin >> n;
    //k = K_SIZE;
   // n = N_SIZE;
    make_myDNA_s(m_DNA, m_DNA_s, k, n);
    m_DNA_s += '\0';

    //ACGT DNA 요소의 정수값과 문자값을 저장해줌
    string_to_short.insert(pair<string, unsigned short>("G", 15));
    string_to_short.insert(pair<string, unsigned short>("C", 256));
    string_to_short.insert(pair<string, unsigned short>("A", 326));
    string_to_short.insert(pair<string, unsigned short>("T", 789));

    input(restruct_m_DNA, k, n);

    double diff = 0;
    for (int i = 0; i < MAX; i++)
        if (restruct_m_DNA[i] != m_DNA[i])
            diff++;
    cout << "복원률: " << 100 - (diff / MAX) * 100 << "%" << endl;//복원한 MY DNA와 원래 MY DNA의 차이

    system("pause");
}

