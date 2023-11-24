#include <iostream>
#include <string>
#include <fstream>
#include <map>
#include <tuple>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "mt19937ar.h"
#include "MemoryOperation.h"
#include "include/stats.hpp"

using namespace std;

string EdgeFile;
int NodeNum;
double Eps;
string Eps_s;
double Mu;
string Mu_s;
char *EpsMu_s[2];

double EpsNsDeg;
double EpsAllbutNsDeg;
double Eps1st, Eps2ndTrSt;

int NSType;
double EClip;
double TClip;
string Clip_s;
int ItrNum;
int Alg;
double Balloc[2];
char *Balloc_s[2];

//struct Node
//{
//    int vi,vj;
//    double x;
//};

// Initialization of statslib
stats::rand_engine_t engine(1776);

FILE *FileOpen(string filename, const char *mode)
{
    FILE *fp;

    if ((fp = fopen(filename.c_str(), mode)) == NULL)
    {
        cout << "cannot open " << filename << endl;
        exit(-1);
    }
    return fp;
}

int compare_double(double *x, double *y)
{
    if (*x > *y)       return(1);   /* return positive integer */
    else if (*x == *y) return(0);   /* return zero     integer */
    else              return(-1);  /* return negative integer */
}

bool checkFileExistence(const std::string& str)
{
    std::ifstream ifs(str);
    return ifs.is_open();
}

// Randomly generate 0, 1, 2, ..., size-1, and store the first num values into rndperm
void MakeRndPerm(int *rndperm, int size, int num)
{
    int rnd;
    int *ordperm;
    int i, j;

    // 0, 1, 2, ..., size-1 --> ordperm
    ordperm = (int *)malloc(size * sizeof(int));
    for (i = 0; i < size; i++)
    {
        ordperm[i] = i;
    }

    for (i = 0; i < num; i++)
    {
        rnd = genrand_int32() % (size - i);
        rndperm[i] = ordperm[rnd];
        for (j = rnd + 1; j < size - i; j++)
        {
            ordperm[j - 1] = ordperm[j];
        }
    }

    free(ordperm);
}

// Read edges from the edge file
void ReadEdges(map<int, int> *a_mat, int *node_order)
{
    int node1, node2;
    int i;
    char s[1025];
    char *tok;
    FILE *fp;

    fp = FileOpen(EdgeFile, "r");
    for(i=0; i<3; i++) fgets(s, 1024, fp);
    while(fgets(s, 1024, fp) != NULL)
    {
        // 1st node --> node1
        tok = strtok(s, ",");
        node1 = atoi(tok);
        // 2nd node --> node2
        tok = strtok(NULL, ",");
        node2 = atoi(tok);
        if(node1 == node2) continue;
        // If both nodes exist, add the edge
        if(node_order[node1] < NodeNum && node_order[node2] < NodeNum)
        {
            a_mat[node_order[node1]][node_order[node2]] = 1;
            a_mat[node_order[node2]][node_order[node1]] = 1;
        }
    }
    fclose(fp);
}

void CalcBaseline(map<int, int> *a_mat, int *deg, string outfile, double &butt_num_ns, double &sen_tri)
{
    map<int, int> *a_mat_ns;			// noisy adjacency matrix
    map<int, int> *a_mat_del;			// deleted adjacency matrix after projection
    map<int, int>::iterator aitr;
    map<int, int>::iterator aitr2;
    map<int, int>::iterator aitr3;
    map<tuple<int, int>, int> a2_mat;
    map<tuple<int, int>, int>::iterator a2_itr;
    double p, q;
    double *butt_num_u, *st2_num_u, *thop_num_u, *butt2_num_u;
    int max_deg;
    int sen_st2_ij;
    int del_num;
    int *rndperm;
    double rnd;
    int i, j, k, x;
    FILE *fp;

    double *deg_ns;
    double sen;
    int deg_ns_floor;

    // Initialization
    a_mat_ns = new map<int, int>[NodeNum];
    a_mat_del = new map<int, int>[NodeNum];
    malloc1D(&butt_num_u, NodeNum);
    malloc1D(&st2_num_u, NodeNum);
    malloc1D(&thop_num_u, NodeNum);
    malloc1D(&butt2_num_u, NodeNum);
    malloc1D(&deg_ns, NodeNum);

    Eps1st = Eps;
    
    // Flip probability --> q
    q = 1.0 / (exp(Eps1st) + 1.0);
    p = 1.0 - q;

    double q1=0.9;

    // RR --> a_mat_ns
    // Count #noisy triangles and #noisy 2-stars for each user --> tri_num_u, st2_num_u
    for(i=0; i<NodeNum; i++)
    {
        for (aitr = a_mat[i].begin(); aitr != a_mat[i].end(); aitr++)
        {
            j = aitr->first;
            for (aitr2 = a_mat[i].begin(); aitr2 != a_mat[i].end(); aitr2++)
            {
                k = aitr2->first;
                // it must be j < k < i.
                if (j >= k||j==i || k==i) continue;

                st2_num_u[i] += 1.0;
                for (int kk=0;kk<NodeNum;++kk)
                {
                    //int kk = aitr3->first;
                    if (kk >= i||j==kk || k==kk) continue;
                    
                    rnd=genrand_real2();
                    
                    if (a_mat[k].count(kk)>=1)
                    {
                        if (rnd < q) a_mat_ns[k][kk]=0; else a_mat_ns[k][kk]=1;
                    }
                    else 
                    {
                        if (rnd < q) a_mat_ns[k][kk]=1; else a_mat_ns[k][kk]=0;
                    }
                    rnd=genrand_real2();

                    if (a_mat[j].count(kk)>=1)
                    {
                        if (rnd < q) a_mat_ns[j][kk]=0;else a_mat_ns[j][kk]=1;
                    }
                    else
                    {
                        if (rnd < q) a_mat_ns[j][kk]=1;else a_mat_ns[j][kk]=0;
                    }

                    if (a_mat[j][i]==1&&a_mat[k][i]==1)
                    {
                        if (a_mat_ns[k][kk] == 1 && a_mat_ns[j][kk] == 1) butt_num_u[i]+=1.0;
                        if (a_mat_ns[k][kk] == 1) thop_num_u[i]+=1.0;
                        if (a_mat_ns[j][kk] == 1) thop_num_u[i]+=1.0;
                    }
                }
            }
        }
        butt2_num_u[i] = (butt_num_u[i] - q*thop_num_u[i])/((1-2*q)*(1-2*q));
    }

    // Empirical estimate --> butt_num_ns
    butt_num_ns = 0;
    for(i=0; i<NodeNum; i++) butt_num_ns += butt2_num_u[i];

    delete[] a_mat_ns;
    delete[] a_mat_del;
    free1D(butt_num_u);
    free1D(st2_num_u);
    free1D(thop_num_u);
    free1D(butt2_num_u);
    free1D(deg_ns);
}

void CalcKstars(map<int, int> *a_mat, int *deg, string outfile, double &butt_num_ns, double &sen_tri)
{
    int *wedge_mat_ns;
    map<int, int>::iterator aitr;
    map<int, int>::iterator aitr2;
    map<int, int>::iterator aitr3;
    map<tuple<int, int>, int> a2_mat;
    map<tuple<int, int>, int>::iterator a2_itr;
    double p, q, q2, q3;
    double *butt_num_u, *st2_num_u, *thop_num_u, *butt2_num_u;
    int max_deg;
    int sen_st2_ij;
    int del_num;
    int *rndperm;
    double rnd;
    int i, j, k, x, kk;
    FILE *fp;

    double *deg_ns;
    double sen;
    int deg_ns_floor;

    // Initialization
    malloc1D(&butt_num_u, NodeNum);
    malloc1D(&st2_num_u, NodeNum);
    malloc1D(&thop_num_u, NodeNum);
    malloc1D(&butt2_num_u, NodeNum);
    malloc1D(&deg_ns, NodeNum);

    Eps1st = Eps;

    // Flip probability --> q
    q = 1.0 / (exp(Eps1st) + 1.0);
    p = 1.0 - q;
    q2 = 0.9;
    q3=0.9;

    // RR --> a_mat_ns
    for (k=0;k<NodeNum;++k)
    {
        for (j=0;j<k;++j)
        {
            malloc1D(&wedge_mat_ns,NodeNum);
            for (i=0;i<NodeNum;++i)
            {
                rnd=genrand_real2();
                if (a_mat[i].count(j)>0&&a_mat[i].count(k)>0)
                {
                    if (rnd<q) wedge_mat_ns[i]=0;else wedge_mat_ns[i]=1;
                }
                else
                {
                    if (rnd<q) wedge_mat_ns[i]=1;else wedge_mat_ns[i]=0;
                }
            }
            for (i=0; i<NodeNum; ++i)
            {
                if (i==j||i==k) continue;
                if (a_mat[i].count(j)<=0||a_mat[i].count(k)<=0) continue;
                st2_num_u[i] += 1.0;
                for (kk=i+1; kk<NodeNum; ++kk)
                {
                    if (j==kk || k==kk) continue;
                    if (wedge_mat_ns[kk]==1) butt_num_u[i]+=1.0;
                }
            }
            free1D(wedge_mat_ns);
        }
    }
    for (i=0;i<NodeNum;++i)
    {
        butt2_num_u[i] = (butt_num_u[i]-q*st2_num_u[i]) / (1.0-2.0*q);
    }

    // Empirical estimate --> butt_num_ns
    butt_num_ns = 0.0;
    for(i=0; i<NodeNum; i++) butt_num_ns += butt2_num_u[i];
    free1D(butt_num_u);
    free1D(st2_num_u);
    free1D(thop_num_u);
    free1D(butt2_num_u);
    free1D(deg_ns);
}

int main(int argc, char *argv[])
{
    int all_node_num;
    int triplet_num;
    int **node_order;
    map<int, int> *a_mat;			// adjacency matrix
    map<int, int>::iterator aitr;
    map<int, int>::iterator aitr2;
    map<int, int>::iterator aitr3;

    int *deg;									// degree
    int *deg_lower;								// degree in the lower-triangular part of a_max
    int max_deg;
    double Balloc_sum = 0.0;
    long long tot_edge_num;
    long long tri_num, st2_num, st3_num, ed2_num, ed1_num, non_num, butt_num;
    map<int, int> pa2_mat;			// 2-path matrix
    long long pa3_num, cy4_num, pa2_pow2;
    double clst;
    double tri_num_ns, sen_tri, butt_num_ns,butt_num_ns1,butt_num_ns2,butt_num_ns3, butt_num_ns4;
    double st2_num_ns, st3_num_ns, sen_st2, sen_st3;
    double st2_num_ns_2, sen_st_2;
    double clst_ns;
    double butt_re_ns, butt_l2_ns;
    double butt_re_ns_avg, butt_l2_ns_avg;
    double butt_re_ns1, butt_l2_ns1;
    double butt_re_ns_avg1, butt_l2_ns_avg1;
    double butt_re_ns2, butt_l2_ns2;
    double butt_re_ns_avg2, butt_l2_ns_avg2;
    double butt_re_ns3, butt_l2_ns3;
    double butt_re_ns_avg3, butt_l2_ns_avg3;
    double butt_re_ns4, butt_l2_ns4;
    double butt_re_ns_avg4, butt_l2_ns_avg4;
    double tri_re_ns, tri_l2_ns;
    double tri_re_ns_avg, tri_l2_ns_avg;
    double st2_re_ns, st2_l2_ns;
    double st2_re_ns_avg, st2_l2_ns_avg;
    double st3_re_ns, st3_l2_ns;
    double st3_re_ns_avg, st3_l2_ns_avg;
    double clst_re_ns, clst_l2_ns;
    double clst_re_ns_avg, clst_l2_ns_avg;
    double eclip_sum, tclip_sum;
    int eclip_num, tclip_num;
    int itr;
    int i, j, k, x;
    string outdir;
    string outfile;
    char s[1025], *str;
    char str_1[] = "1";
    char *tok;
    FILE *fp;

    int fix_perm;

    // Initialization of Mersennne Twister
    unsigned long init[4] = { 0x123, 0x234, 0x345, 0x456 }, length = 4;
    init_by_array(init, length);

    if (argc < 2)
    {
        printf("Usage: %s [EdgeFile (in)] ([#nodes (default: -1)] [epsilon-mu/mu2/mu3 (default: 1-1)] [NSType (default: -1)] [tclip-eclip (default: -1)] [#itr(-1) (default: 1)] [alg (default: 0)] [Balloc (default: 1-1)])\n\n", argv[0]);
        printf("[EdgeFile]: Edge file\n");
        printf("[#nodes]: Number of nodes (-1: all)\n");
        printf("[epsilon-mu/mu2/mu3]: Parameters epsilon and mu/mu2/mu3 (I/II/III) (mu/mu2/mu3 = exp(Eps1st)/(exp(Eps1st)+1) when it is set to 1)\n");
        printf("[NSType]: Noise type (-1: no noise, 0: Lap (max degree), 1: Lap (degree + clip), 2: Lap (noisy degree + clip), 3: Lap (noisy degree + clip + ignore users))\n");
        printf("[tclip(-eclip)]: Triangle and edge clipping parameters (-1: no clipping) (alg=2-4; set eclip when NSType=2)\n");
        printf("[#itr(-1)]: Number of iterations (set #itr-1 to fix the permutation of nodes)\n");
        printf("[alg]: Algorithm (1: interactive local, 2: efficient interactive local I, 3: efficient interactive local II, 4: efficient interactive local III, 5: non-interactive local (RR w/ emp), 6: non-interactive local (RR w/o emp), 7: [Ye+, T-KDE (mean)], 8: [Ye+, T-KDE (median)], 9: [Ye+, T-KDE (most frequent degree)], 10: non-interactive local (ARR w/ emp))\n");
        printf("[Balloc]: Privacy budget allocation (alg=1-3): Eps1st-Eps2ndTrSt\n");
        return -1;
    }

    EdgeFile = argv[1];

    NodeNum = -1;
    if (argc >= 3) NodeNum = atoi(argv[2]);

    Eps = 1.0;
    Eps_s = "1";
    Mu = 1.0;
    Mu_s = "1";

    if (argc >= 4)
    {
        if((EpsMu_s[0] = strtok(argv[3], "-")) == NULL)
        {
            printf("Error: incorrect [epsilon-mu]\n");
            exit(-1);
        }
        if((EpsMu_s[1] = strtok(NULL, "-")) == NULL)
        {
            printf("Error: incorrect [epsilon-mu]\n");
            exit(-1);
        }
        Eps = atof(EpsMu_s[0]);
        Mu = atof(EpsMu_s[1]);
        Eps_s = EpsMu_s[0];
        Mu_s = EpsMu_s[1];
    }

    NSType = -1;
    if (argc >= 5) NSType = atoi(argv[4]);

    TClip = -1;
    EClip = -1;
    Clip_s = "-1";
    if (argc >= 6)
    {
        Clip_s = argv[5];
        if(strcmp(argv[5], "-1") != 0)
        {
            if(argv[5][0] !=  '-')
            {
                if((tok  = strtok(argv[5], "-")) == NULL)
                {
                    printf("Error: incorrect [tclip(-eclip)]\n");
                    exit(-1);
                }
                TClip = atof(tok);
                if((tok  = strtok(NULL, "-")) != NULL) EClip = atof(tok);
            }
            else
            {
                if((tok  = strtok(argv[5], "-")) == NULL)
                {
                    printf("Error: incorrect [tclip(-eclip)]\n");
                    exit(-1);
                }
                if((tok  = strtok(NULL, "-")) != NULL) EClip = atof(tok);
            }
        }

    }

    ItrNum = 1;
    fix_perm = 0;
    if (argc >= 7)
    {
        tok  = strtok(argv[6], "-");
        ItrNum = atoi(tok);
        if((tok  = strtok(NULL, "-")) != NULL)
        {
            if (strcmp(tok, "1") != 0)
            {
                printf("Error: incorrect [#itr(-1)]\n");
                exit(-1);
            }
            else fix_perm = 1;
        }
    }

    Alg = 0;
    if (argc >= 8) Alg = atoi(argv[7]);
    if (Alg <= 0 || Alg > 10)
    {
        printf("Error: incorrect [Alg]\n");
        exit(-1);
    }

    for(i=0; i<2; i++)
    {
        Balloc[i] = 1.0;
        Balloc_s[i] = str_1;
    }
    if (argc >= 9)
    {
        if((Balloc_s[0] = strtok(argv[8], "-")) == NULL)
        {
            printf("Error: incorrect [Balloc]\n");
            exit(-1);
        }
        Balloc[0] = atof(Balloc_s[0]);
        if((Balloc_s[1] = strtok(NULL, "-")) == NULL)
        {
            printf("Error: incorrect [Balloc]\n");
            exit(-1);
        }
        Balloc[1] = atof(Balloc_s[1]);
    }

    // Privacy budget allocation
    for(i=0; i<2; i++) Balloc_sum += Balloc[i];
    if(NSType == -1 || NSType == 0 || NSType == 1)
    {
        EpsNsDeg = 0.0;				// Epsilon for calculating the noisy degree
        EpsAllbutNsDeg = Eps;
        Eps1st = Eps * Balloc[0] / Balloc_sum;
        Eps2ndTrSt = Eps * Balloc[1] / Balloc_sum;
    }
    else if(NSType == 2 || NSType == 3)
    {
        EpsNsDeg = Eps / 10;		// Epsilon for calculating the noisy degree
        EpsAllbutNsDeg = Eps - EpsNsDeg;
        Eps1st = EpsAllbutNsDeg * Balloc[0] / Balloc_sum;
        Eps2ndTrSt = EpsAllbutNsDeg * Balloc[1]  / Balloc_sum;
    }

    // mu/mu^2/mu^3 = exp(Eps1st)/(exp(Eps1st)+1) when it is set to 1
    if(Mu == 1.0)
    {
        Mu = exp(Eps1st) / (exp(Eps1st) + 1);
    }
    // When Alg == 3 (efficient interactive local II), calculate mu from mu^2
    if (Alg == 3)
    {
        Mu = pow(Mu, 1.0/2.0);
    }
    // When Alg == 4 (efficient interactive local III), calculate mu from mu^3
    if (Alg == 4)
    {
        Mu = pow(Mu, 1.0/3.0);
    }

    // Total number of nodes --> all_node_num
    fp = FileOpen(EdgeFile, "r");
    for(i=0; i<2; i++) fgets(s, 1024, fp);
    all_node_num = atoi(s);
    fclose(fp);

    // malloc
    malloc2D(&node_order, ItrNum, all_node_num);

    // Use all nodes
    if (NodeNum == -1)
    {
        NodeNum = all_node_num;
        for(j=0; j<NodeNum; j++) node_order[0][j] = j;
    }
    // Randomly generate the order of nodes --> node_order
    else
    {
        i = EdgeFile.find_last_of("/");
        outdir = EdgeFile.substr(0, i+1);
        outfile = outdir + "node-order_itr" + to_string(ItrNum) + ".csv";
        if(checkFileExistence(outfile))
        {
            fp = FileOpen(outfile, "r");
            for(j=0; j<all_node_num; j++)
            {
                fgets(s, 1024, fp);
                strtok(s, ",");
                for(i=0; i<ItrNum; i++)
                {
                    node_order[i][j] = atoi(strtok(NULL, ","));
                }
            }
            fclose(fp);
        }
        else
        {
            for(i=0; i<ItrNum; i++)
            {
                MakeRndPerm(node_order[i], all_node_num, all_node_num);
            }
            fp = FileOpen(outfile, "w");
            for(j=0; j<all_node_num; j++)
            {
                fprintf(fp, "%d,", j);
                for(i=0; i<ItrNum; i++) fprintf(fp, "%d,", node_order[i][j]);
                fprintf(fp, "\n");
            }
            fclose(fp);
        }

        // Use only the first permutation
        if (fix_perm)
        {
            for(j=0; j<all_node_num; j++)
            {
                for(i=1; i<ItrNum; i++) node_order[i][j] = node_order[0][j];
            }
        }
    }

    // #triplet --> triplet_num
    triplet_num = NodeNum*(NodeNum-1)*(NodeNum-2)/6;

    // Initialization
    malloc1D(&deg, NodeNum);
    malloc1D(&deg_lower, NodeNum);
    butt_l2_ns_avg=butt_re_ns_avg=0.0;
    butt_l2_ns_avg1=butt_re_ns_avg1=0.0;
    butt_l2_ns_avg2=butt_re_ns_avg2=0.0;
    butt_l2_ns_avg3=butt_re_ns_avg3=0.0;
    butt_l2_ns_avg4=butt_re_ns_avg4=0.0;
    tri_re_ns_avg = tri_l2_ns_avg = 0.0;
    st2_re_ns_avg = st2_l2_ns_avg = 0.0;
    st3_re_ns_avg = st3_l2_ns_avg = 0.0;
    clst_re_ns_avg = clst_l2_ns_avg = 0.0;

    // Output the header
    i = EdgeFile.find_last_of("/");
    outdir = EdgeFile.substr(0, i+1);
    for(i=0; i<3; i++)
    {
        if(fix_perm) outfile = outdir + "res_n" + to_string(NodeNum) + "_alg" + to_string(Alg) + "_eps" + Eps_s + "-" + Mu_s + "_ns" + to_string(NSType) + "_cl" + Clip_s + "_ba" + Balloc_s[0] + "-" + Balloc_s[1] + "_itr" + to_string(ItrNum) + "-1.csv";
        else outfile = outdir + "res_n" + to_string(NodeNum) + "_alg" + to_string(Alg) + "_eps" + Eps_s + "-" + Mu_s + "_ns" + to_string(NSType) + "_cl" + Clip_s + "_ba" + Balloc_s[0] + "-" + Balloc_s[1] + "_itr" + to_string(ItrNum) + ".csv";
        fp = FileOpen(outfile, "w");
        fprintf(fp, "#butt0(true),#butt0(est),#butt0(rel-err),#butt0(l2-loss),#butt1(true),#butt1(est),#butt1(rel-err),#butt1(l2-loss),#butt2(true),#butt2(est),#butt2(rel-err),#butt2(l2-loss),#butt3(true),#butt3(est),#butt3(rel-err),#butt3(l2-loss),#butt4(true),#butt4(est),#butt4(rel-err),#butt4(l2-loss)\n");
        fclose(fp);
    }

    // For each iteration
    for(itr=0; itr<ItrNum; itr++)
    {
        // Read edges for each iteration when NodeNum < all_node_num
        if(NodeNum < all_node_num || itr == 0)
        {
            // Initialization
            a_mat = new map<int, int>[NodeNum];

            // Read edges from the edge file --> a_mat
            ReadEdges(a_mat, node_order[itr]);
            // Degree --> deg
            for(i=0; i<NodeNum; i++) deg[i] = 0;
            for(i=0; i<NodeNum; i++)
            {
                for (aitr = a_mat[i].begin(); aitr != a_mat[i].end(); aitr++) deg[i] += 1;
            }

            // Degree --> deg_lower
            for(i=0; i<NodeNum; i++) deg_lower[i] = 0;
            for(i=0; i<NodeNum; i++)
            {
                for (aitr = a_mat[i].begin(); aitr != a_mat[i].end(); aitr++)
                {
                    if(aitr->first < i) deg_lower[i] += 1;
                }
            }

            // max(deg) --> max_deg
            max_deg = 0;
            for(i=0; i<NodeNum; i++)
            {
                if(max_deg < deg[i]) max_deg = deg[i];
            }

            // Total number of edges --> tot_edge_num
            tot_edge_num = 0;
            for(i=0; i<NodeNum; i++) tot_edge_num += (long long)deg[i];
            tot_edge_num /= 2;

            // (p,q)-clique -> butt_num
            butt_num=0;
            for(i=0; i<NodeNum; i++)
            {
                for (aitr = a_mat[i].begin(); aitr != a_mat[i].end(); aitr++)
                {
                    j = aitr->first;
                    for (aitr2 = a_mat[i].begin(); aitr2 != a_mat[i].end(); aitr2++)
                    {
                        k = aitr2->first;
                        if (j >= k || i==j || i==k) continue;
                        for (aitr3 = a_mat[k].begin(); aitr3 != a_mat[k].end(); aitr3++)
                        {
                            int kk = aitr3->first;
                            if (k==kk || kk>=i || j==kk) continue;
                            if (a_mat[j].count(kk) >0) butt_num++;
                        }
                    }
                    double p_sp=0.7;
                }
            }

            tri_num = 0;

            // #2-stars, #3-stars --> st2_num, st3_num
            st2_num = st3_num = 0;
            for(i=0; i<NodeNum; i++)
            {
                st2_num += ((long long)deg[i] * ((long long)deg[i]-1)) / 2;
                st3_num += ((long long)deg[i] * ((long long)deg[i]-1) * ((long long)deg[i]-2)) / 6;
            }
        }

        /************************ Calculate sub-graph counts ************************/
        eclip_sum = tclip_sum = 0.0;
        eclip_num = tclip_num = 0;
        // Interactive (2-rounds) local
        
        butt_num_ns = butt_num_ns4=0.0;
        CalcBaseline(a_mat, deg_lower, outfile, butt_num_ns, sen_tri);
        
        CalcKstars(a_mat,deg_lower, outfile, butt_num_ns4, sen_tri);


        /**************************** Evaluate the loss *****************************/
        // relative error -> butt_re_ns
        butt_re_ns = fabs(butt_num_ns-(double)butt_num)/max((double)butt_num, 0.001*NodeNum);
        butt_re_ns_avg+=butt_re_ns;

        butt_l2_ns = (butt_num_ns - (double)butt_num)*(butt_num_ns - (double)butt_num);
        butt_l2_ns_avg+=butt_l2_ns;

        butt_re_ns4 = fabs(butt_num_ns4-(double)butt_num)/max((double)butt_num, 0.001*NodeNum);
        butt_re_ns_avg4+=butt_re_ns4;
        butt_l2_ns4 = (butt_num_ns4 - (double)butt_num)*(butt_num_ns4 - (double)butt_num);
        butt_l2_ns_avg4+=butt_l2_ns4;

        /**************************** Output the results ****************************/
        fp = FileOpen(outfile, "a");
        fprintf(fp, "%e,%e,%e,%e,%e,%e,%e,%e\n",
                (double)butt_num, butt_num_ns, butt_re_ns, butt_l2_ns,
                (double)butt_num, butt_num_ns4, butt_re_ns4, butt_l2_ns4);
        fclose(fp);

        if(NodeNum < all_node_num || itr == ItrNum - 1)
        {
            delete[] a_mat;
        }
    }

    /************************* Output the results (AVG) *************************/
    butt_l2_ns_avg /= (double)ItrNum;
    butt_re_ns_avg /= (double)ItrNum;
    
    butt_l2_ns_avg4 /= (double)ItrNum;
    butt_re_ns_avg4 /= (double)ItrNum;

    fp = FileOpen(outfile, "a");
    fprintf(fp, "function,AVG(rel-err),AVG(l2-loss)\n");
    fprintf(fp, "Baseline,%e,%e\n", butt_re_ns_avg, butt_l2_ns_avg);
    fprintf(fp, "k-stars LDP,%e,%e\n", butt_re_ns_avg4, butt_l2_ns_avg4);
    fclose(fp);

    // free
    free2D(node_order, ItrNum);
    free1D(deg);
    free1D(deg_lower);

    return 0;
}
