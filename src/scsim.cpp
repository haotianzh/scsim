// Not exactly a script, but I have to do it this way since AWK does not work
#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <string>
#include <sys/types.h>
#include <time.h>
// #include <unistd.h>
#include <cmath>
#include <map>
#include <set>
#include <queue>
#include <random>
#ifdef _WIN32
    #include <io.h>
    #include <process.h>
    #include <windows.h> // Often needed for Sleep() if you used sleep()
    
    // Windows usually maps 'unistd' functions to underscore versions
    // If you get errors about 'access' or 'F_OK', uncomment these:
    // #define access _access
    // #define F_OK 0
#else
    #include <unistd.h>
#endif


using namespace std;


// parameters
const int MIN_NUM_PARAMS = 10;
const int MAX_NUM_PARAMS = 17;
static double rateCNVInc = 0.00;
static double rateCNVDec = 0.00;
static double rateDropout = 0.1;
static double rateDropoutstd = 0.0;
static double rateErr = 0.0;
static int numDoublet = 0;
static double fracRecurrentMut = 0.0;
static int numAllelesPerCopy = 5;
static int minCNV = 1;
static int maxCNV = 1;
static int defCNV = 1;
//static double prob0Const = 0.8;
//static double prob1Const = 1.0;
static double fracMissing = 0.0;
static bool seedRndUserDef = true;
static int seedRndUser = 0xF1CD1321;    // dec: 4056748833
//static double boostFac = 1.0;       // for dropout (only) one, tweak the prob to hint on the change; for unchanged one, do the opposite. That is, if a 0 is a real 0, inc the prob of being zero; if 0 is a dropout, dec the prob of being zero
static double minAlleFreq = 0.0;
const double MIN_BR_LEN = 0.0001;
const double SMALL_NUM = 0.00001;
const double SMALL_PROB = 0.000000001;
const double LARGE_PROB = 0.999999999;
static bool fBinary = true;
// the following refers to the single-cell sequencing realted aspects
static double aveReadDepth = 4.0;
static double stdReadDepth = 2.0;
static double aveStrandBias = 0.0;
static double stdStrandBias = 0.0;
static vector<double> listCellDropoutRates;
static bool fBetaBinomial=false;
static int betaparama=2;
static int betaparamb=2;

//*******************************************************************************************************
// Utilies

static int QSortCompareInt( const void *arg1, const void *arg2 )
{
    /* Compare all of both strings: */
    // assume sorting in accending order
    int n1 = *((int *) arg1);
    int n2 = *((int *) arg2);
    //cout <<"arg1 = " << n1 << ", arg2 = " << n2 << endl;
    if( n1 > n2)
    {
        return 1;
    }
    else if( n1 < n2)
    {
        return -1;
    }
    else
    {
        return 0;
    }
}
void SortIntVec( vector<int> &vecVals, int start, int end )
{
    //#if 0
    if( vecVals.size() <= 1)
    {
        // do nothing
        return;
    }
    if (end < 0 )
    {
        end = vecVals.size() - 1;
    }
    int sortLen = end - start +1;
    int *array = new int[sortLen];
    for(int i=start; i<= end; ++i)
    {
        array[i-start] = vecVals[i];
    }
    qsort( (void *)array, sortLen, sizeof( int ), QSortCompareInt );
    // Now write back
    for(int i=start; i<=end; ++i)
    {
        vecVals[i] = array[i-start];
    }
    delete [] array;
}
static int FindTwoClusterBoundary( const vector<int> &listVals )
{
    // goal: find the cutoff point s.t. the two parts are as similar as possible
    // return the largest value in part1 (i.e. anything <= this value is the second part)
    vector<int> listValsUse = listVals;
    SortIntVec(listValsUse, 0, listValsUse.size()-1);
    //
    int res = -1;
    double minCenterSum = 1000000000000000000.0;
    for(int i=1; i<(int)listValsUse.size()-1; ++i)
    {
        // compute the partition
        int szPart1 = i, szPart2 = (int)listValsUse.size()-szPart1;
        if(szPart2 == 0)
        {
            szPart2 = 1;
        }
        int sumPart1 = 0, sumPart2 = 0;
        for(int j=0; j<i; ++j)
        {
            sumPart1 += listValsUse[j];
        }
        for(int j=i; j<(int)listValsUse.size(); ++j)
        {
            sumPart2 += listValsUse[j];
        }
        double partCenter1 = ((double)sumPart1/szPart1);
        double partCenter2 = ((double)sumPart2/szPart2);
        double sumSqDiff1 = 0.0, sumSqDiff2 = 0.0;
        for(int j=0; j<i; ++j)
        {
            double diff = listValsUse[j] - partCenter1;
            sumSqDiff1 += diff*diff;
        }
        for(int j=i; j<(int)listValsUse.size(); ++j)
        {
            double diff = listValsUse[j] - partCenter2;
            sumSqDiff2 += diff*diff;
        }
        double meanSqDiff1 = sumSqDiff1/szPart1;
        double meanSqDiff2 = sumSqDiff2/szPart2;
        double meanSqDiffSum = meanSqDiff1 + meanSqDiff2;
        if( meanSqDiffSum < minCenterSum )
        {
            minCenterSum = meanSqDiffSum;
            res = listValsUse[i-1];
        }
    }
    return res;
}


//******************************************************************************************************
// define cell tree

class CellTreeNode
{
public:
    CellTreeNode() : cellId(-1), pParent(NULL), brLen(-1.0) {}
    ~CellTreeNode() {
        for(int i=0; i<(int)listChildren.size(); ++i) {
            delete listChildren[i];
        }
    }
    
    int GetCellId() const { return cellId; }
    void SetCellId(int cid) { cellId = cid; }
    CellTreeNode *GetParent() const { return pParent; }
    void SetParent(CellTreeNode *pParentIn) { pParent = pParentIn; }
    int GetNumChildren() const { return listChildren.size(); }
    CellTreeNode *GetChild(int i) const { return listChildren[i]; }
    void AddChild(CellTreeNode *pChild) { listChildren.push_back(pChild); pChild->SetParent( const_cast<CellTreeNode *>( this )); }
    bool IsLeaf() const { return listChildren.size() == 0; }
    string GetNewick() const {
        string sres;
        if( IsLeaf() ) {
            sres = std::to_string( GetCellId() );
        }
        else{
            string res = "(";
            for(int i=0; i<(int)listChildren.size(); ++i) {
                string strStep = listChildren[i]->GetNewick( );
                res += strStep;
                if( i <(int)listChildren.size()-1) {
                    res += ",";
                }
            }
            res += ")";
            sres = res;
        }
        if( GetLength() >= 0.0)
        {
            sres += ":" + std::to_string(GetLength());
        }
        return sres;
    }
    double GetLength() const { return brLen; }
    void Setlength(double len) { brLen = len; }
    void GetAllNodesBelow( vector<CellTreeNode *> &listNodes ) const {
        listNodes.push_back( const_cast<CellTreeNode *>(this) );
        for(int i=0; i<GetNumChildren(); ++i) {
            GetChild(i)->GetAllNodesBelow( listNodes );
        }
    }
    int GetNumLeaves() const {
        vector<CellTreeNode *> listNodes;
        GetAllNodesBelow(listNodes);
        int res = 0;
        for(int i=0; i<(int)listNodes.size(); ++i) {
            if( listNodes[i]->IsLeaf() ) {
                ++res;
            }
        }
        return res;
    }
    
private:
    int cellId;
    CellTreeNode *pParent;
    vector<CellTreeNode *> listChildren;
    double brLen;
};


// construct cell tree recurisvely
CellTreeNode * ProcTreeStrPos(string &strTree, int pos, int &posEnd)
{
//cout << "Construct subtree at: " << pos << endl;
    
    CellTreeNode *pres = NULL;
    // if a subtree, do it recurisvely
    if( strTree[pos] != '(')
    {
        // find the : or )
        std::size_t found1 = strTree.find_first_of(":", pos);
        if( found1 == std::string::npos)
        {
            std::size_t found2 = strTree.find_first_of(",", pos);
            std::size_t found3 = strTree.find_first_of(")", pos);
            if( found2 == std::string::npos)
            {
                found1 = found3;
            }
            else if( found3 == std::string::npos )
            {
                found1 = found2;
            }
            else
            {
                found1 = std::min(found2, found3);
            }
        }

        if( found1 == std::string::npos)
        {
            cout << "Wrong format. Stop\n";
            exit(1);
        }
        //
        std::size_t fposTaxon = found1;
        std::size_t found4 = strTree.find_first_of( ":", pos );
        if( found4 != std::string::npos && found4 < fposTaxon)
        {
            fposTaxon = found4;
        }
        string sid = strTree.substr( pos, fposTaxon-pos+1 );
        posEnd = found1-1;
        int cid = atoi( sid.c_str() );
        CellTreeNode *pleaf = new CellTreeNode;
        pleaf->SetCellId(cid);
//cout << "Leaf: posend is " << posEnd << endl;
        pres = pleaf;
    }
    else
    {
        // construt a subtree, first find the matching other ")"
        int netParenth = 1;
        posEnd = -1;
        for( int p=pos+1; p<strTree.length(); ++p )
        {
            if( strTree[p] == '(')
            {
                ++netParenth;
            }
            else if( strTree[p] == ')' )
            {
                --netParenth;
            }
            
            if( netParenth == 0 )
            {
                posEnd = p;
                break;
            }
        }
        
        if(posEnd < 0 )
        {
            cout << "Wrong format2. Stop\n";
            exit(1);
        }
//cout << "Internal: posEnd: " << posEnd << endl;
        //
        std::size_t posCur = pos;
        if( strTree[posCur] == '(')
        {
            ++posCur;
        }
        CellTreeNode *pNode = new CellTreeNode;
        netParenth = 1;
        while( posCur < posEnd  )
        {
            if( strTree[posCur] == ',')
            {
                ++posCur;
                continue;
            }
            if( strTree[posCur] == ':' )
            {
                break;
            }
            
//cout << "Reconstruct subtree at: " << posCur << endl;
            // read one underlying subtree
            int posEndChild= -1;
            CellTreeNode *pChild = ProcTreeStrPos(strTree, posCur, posEndChild);
            if( posEndChild == std::string::npos)
            {
                break;
            }
            pNode->AddChild(pChild);
            posCur = posEndChild+1;
        }
        
        pres = pNode;
    }
    
    // set up br length if needed
    if( posEnd < strTree.length()-1 )
    {
        std::size_t posCur = posEnd+1;
        // read in the br len if any
        if( strTree[posCur] == ':' || strTree[posCur] == ')' )
        {
            std::size_t found2 = strTree.find_first_of(",", posCur);
            std::size_t found3 = strTree.find_first_of(")", posCur);
            std::size_t pe = std::min(found2, found3);
            string slen = strTree.substr( posCur+1, pe-posCur );
            double brLen = atof( slen.c_str() );
            pres->Setlength(brLen);
    
            posEnd = pe;
        }
    }
    
    
    return pres;
}


// process a string from a particular position
CellTreeNode *ProcTreeStr(string &strTree)
{
    // find the first (
    int pos = 0;
    while( strTree[pos] != '(')
    {
        ++pos;
    }
    
    int posEnd = -1;
    CellTreeNode *pTreeRoot = ProcTreeStrPos(strTree, pos, posEnd);
    return pTreeRoot;
}

//******************************************************************************************************
// CNV

class CellSeqCopy
{
public:
    CellSeqCopy() {}
    CellSeqCopy *GetParent() const { return pParent; }
    void SetParent(CellSeqCopy *pp) { pParent = pp; }
    void AddChild(CellSeqCopy *pc) { listChildren.push_back(pc); }
    int GetNumChildren() const { return listChildren.size(); }
    CellSeqCopy *GetChild(int i) const { return listChildren[i]; }
    CellSeqCopy *MakeCopy() {
        //
        CellSeqCopy *pNewCopy = new CellSeqCopy;
        pNewCopy->SetParent( this );
        ;
        return pNewCopy;
    }
    
private:
    CellSeqCopy *pParent;
    vector<CellSeqCopy *> listChildren;
};

//******************************************************************************************************
// Simulation

static void InitRand(long seedRnd)
{
    if( seedRndUserDef == true || seedRnd >= 0 )
    {
        srand(seedRnd);
    }
    else
    {
        srand (time(NULL));
    }
}
static double RndFrac()
{
    return ((double) rand() / (RAND_MAX));
}
static bool RndEvt(double rate)
{
    double rnd = RndFrac();
    if( rnd<rate )
    {
        return true;
    }
    else
    {
        return false;
    }
}
static int RndInt(int rangeMin, int rangeMax)
{
//cout << "RndInt: min: " << rangeMin << ", rangeMax: " << rangeMax << endl;
    return rangeMin + (int)( (rangeMax-rangeMin+1)*RndFrac() );
}
// This functionreturn a weighted uniformly item index from the list
static int GetWeightedRandItemIndex( const vector<double> &itemWeights )
{
//cout << "GetWeightedRandItemIndex: ";
//for(int i=0; i<(int)itemWeights.size(); ++i)
//{
//cout << itemWeights[i] << " ";
//}
//cout << endl;
    double accum = 0.0;
    for( unsigned int i=0; i<itemWeights.size(); ++i )
    {
        //cout << "one weight = " << itemWeights[i] << endl;
        accum += itemWeights[i];
    }
    //YW_ASSERT_INFO( accum > 0.0000001, "3. Can not be too small" );
    double frac = RndFrac();
    double curFract = 0.0;
    for( unsigned int i=0; i<itemWeights.size(); ++i )
    {
        curFract += itemWeights[i]/accum;
        if( curFract >= frac )
        {
            return i;
        }
    }
    // Can not come here
    cout << "Something wrong here\n";
    return -1;      // should nothappen
}

static double SampleDropoutRateForCell()
{
// HACK: Oct 10, 2025, YW
// for now, use the same droput rate
return rateDropout;
    
    static bool fInit = false;
    //std::random_device rd{};
    static std::mt19937 gen(time(0));
    if( fInit == false)
    {
        //
        const int RAND_SEED1 = 18973;
        gen.seed( RAND_SEED1 );
        fInit = true;
    }
    
    // values near the mean are the most likely
    // standard deviation affects the dispersion of generated values from the mean
    std::normal_distribution<double> d(rateDropout, rateDropoutstd);
    
    double res = d(gen);
    if( res < SMALL_PROB)
    {
        res = SMALL_PROB;
    }
    else if( res > LARGE_PROB)
    {
        res = LARGE_PROB;
    }
    return res;
}


void SimulateCNVOnTree( CellTreeNode *pTreeRoot, map<CellTreeNode *, vector<CellSeqCopy *> > &listCopies )
{
    // root node has 2 copies by default
    listCopies.clear();
    CellSeqCopy *pr1 = new CellSeqCopy;
    CellSeqCopy *pr2 = new CellSeqCopy;
    listCopies[pTreeRoot].push_back(pr1);
    listCopies[pTreeRoot].push_back(pr2);
    
    //
    queue<CellTreeNode *> queueNodesProc;
    for(int i=0; i<pTreeRoot->GetNumChildren(); ++i)
    {
        queueNodesProc.push(pTreeRoot->GetChild(i));
    }
    while( queueNodesProc.empty() == false )
    {
        CellTreeNode *pNode = queueNodesProc.front();
        queueNodesProc.pop();
        for(int i=0; i<pNode->GetNumChildren(); ++i)
        {
            queueNodesProc.push(pNode->GetChild(i));
        }
        
//cout << "Processing node in CNV simulaton: " << pNode->GetNewick() << endl;
        CellTreeNode *pPar = pNode->GetParent();
        if( pPar == NULL )
        {
            cout << "WRONG.\n";
            exit(1);
        }
        
        // first copy all
        for(int i=0; i<(int)listCopies[pPar].size(); ++i)
        {
            listCopies[pNode].push_back(listCopies[pPar][i]->MakeCopy());
        }
        
        //
        double len = pNode->GetLength();
//cout << "Len: " << len << ", number of copies at this node: " << listCopies[pNode].size() << endl;
        bool fChange = RndEvt( len*(rateCNVInc+rateCNVDec) );
        if( fChange)
        {
//cout << "Change: \n";
            bool fInc = RndEvt( rateCNVInc/(rateCNVInc+rateCNVDec) );
            if( fInc )
            {
//cout << "Inc\n";
                // random inc based a copy of parent
                if( listCopies[pPar].size() > 0 )
                {
                    int indexCpy = RndInt(0, listCopies[pPar].size()-1);
//cout << "indexCpy: " << indexCpy << endl;
                    CellSeqCopy *pc = listCopies[pPar][indexCpy]->MakeCopy();
                    listCopies[ pNode ].push_back( pc );
                }
            }
            else
            {
//cout << "Dec\n";
                if( listCopies[pPar].size() > 0 )
                {
                    int indexCpy = RndInt(0, listCopies[pPar].size()-1);
//cout << "indexCpy: " << indexCpy << endl;
                    listCopies[ pNode ].erase( listCopies[pNode].begin()+indexCpy );
                }
            }
        }
    }
    
    // print out current copies
    //cout << "Copy numbers of extant cells\n";
    map<int,int> mapCellCopyNums;
    for( map<CellTreeNode *, vector<CellSeqCopy *> > :: iterator it = listCopies.begin(); it != listCopies.end(); ++it )
    {
        if( it->first->IsLeaf() )
        {
            mapCellCopyNums[ it->first->GetCellId() ] = it->second.size();
            //cout << "Cell " << it->first->GetCellId() << ": number of copies is " << it->second.size() << endl;
        }
    }
    for(map<int,int> :: iterator it = mapCellCopyNums.begin(); it != mapCellCopyNums.end(); ++it)
    {
        //cout << "Site: " << site << "Cell " << it->first << ": number of copies is " << it->second << endl;
    }
    
    //for( map<CellTreeNode *, vector<CellSeqCopy *> > :: iterator it = listCopies.begin(); it != listCopies.end(); ++it )
    //{
    //    cout << "Node " << it->first->GetNewick() << ": number of copies is " << it->second.size() << endl;
    //}
}

static bool IsRecurrMut()
{
    double rf = RndFrac();
    if( rf < fracRecurrentMut )
    {
        return true;
    }
    else
    {
        return false;
    }
}

bool SimulatePointMutsOnRecur( CellTreeNode *pTreeRoot, const map<CellTreeNode *, vector<CellSeqCopy *> > &listCopies, map<int, pair<int,int> > &mapCellGenos );

static void SimulatePointMutsOnSingleNoRecur( CellTreeNode *pTreeRoot, const map<CellTreeNode *, vector<CellSeqCopy *> > &listCopies, map<int, pair<int,int> > &mapCellGenos )
{
    map<CellTreeNode *, vector<CellSeqCopy *> > &listCopiesUse = const_cast<map<CellTreeNode *, vector<CellSeqCopy *> > &>(listCopies);
    
    mapCellGenos.clear();
    map<CellSeqCopy *, int> mapMuts;
    
    vector<CellTreeNode *> listNodes;
    pTreeRoot->GetAllNodesBelow( listNodes );
    vector<double> listLens;
    bool fZeroLen = true;
    for(int i=0; i<(int)listNodes.size(); ++i)
    {
        double len = listNodes[i]->GetLength();
        if(len <= 0.0 )
        {
            len = 0.0;
        }
        else
        {
            fZeroLen = false;
        }
        
        // skip all leaf
        if( listNodes[i]->IsLeaf() )
        {
            //    len = 0.0;
        }
        
        listLens.push_back(len);
    }
    
    // setup uniform length if length is not specified
    if( fZeroLen )
    {
        for(int i=0; i<(int)listLens.size(); ++i)
        {
            listLens[i] = 1.0;
        }
        fZeroLen = false;
    }
    
    int indexMut = GetWeightedRandItemIndex( listLens );
    //cout << "indexMut: " << indexMut << endl;
    CellTreeNode *pNodeMut = listNodes[indexMut];
    //cout << "Mutation at: " << pNodeMut->GetNewick() << endl;
    
    // now generate muts; root always has zero
    CellTreeNode *pNodeCur = pTreeRoot;
    for(int i=0; i<(int)listCopiesUse[pTreeRoot].size(); ++i)
    {
        mapMuts[ listCopiesUse[pTreeRoot][i] ]  = 0;
    }
    queue<CellTreeNode *> queueNodesProc;
    for(int i=0; i<pTreeRoot->GetNumChildren(); ++i)
    {
        queueNodesProc.push(pTreeRoot->GetChild(i));
    }
    while( queueNodesProc.empty() == false )
    {
        CellTreeNode *pNode = queueNodesProc.front();
        queueNodesProc.pop();
        for(int i=0; i<pNode->GetNumChildren(); ++i)
        {
            queueNodesProc.push(pNode->GetChild(i));
        }
        
        //cout << "Processing node in point mutation simulaton: " << pNode->GetNewick() << endl;
        CellTreeNode *pPar = pNode->GetParent();
        
        // copy parent if no mutations occur
        if( pNode != pNodeMut )
        {
            for(int i=0; i<(int)listCopiesUse[pNode].size(); ++i )
            {
                CellSeqCopy *pCpyPar = listCopiesUse[pNode][i]->GetParent();
                int mut = mapMuts[ pCpyPar ];
                mapMuts[ listCopiesUse[pNode][i] ] = mut;
            }
        }
        else
        {
            // pick a random copy to mut;
            int mutCopyIndex = RndInt( 0, listCopiesUse[pNode].size()-1 );
            for(int i=0; i<(int)listCopiesUse[pNode].size(); ++i )
            {
                CellSeqCopy *pCpyPar = listCopiesUse[pNode][i]->GetParent();
                int mutPar = mapMuts[ pCpyPar ];
                if( i != mutCopyIndex )
                {
                    mapMuts[ listCopiesUse[pNode][i] ] = mutPar;
                }
                else
                {
                    int mut = 1;
                    if(mutPar == 1)
                    {
                        mut = 0;
                    }
                    mapMuts[ listCopiesUse[pNode][i] ] = mut;
                }
            }
        }
    }
    //cout << "Now collect genotypes...\n";
    //
    for(int i=0; i<(int)listNodes.size(); ++i )
    {
        if( listNodes[i]->IsLeaf() == false )
        {
            continue;
        }
        int n0 = 0, n1 = 0;
        for( int j=0; j<(int)listCopiesUse[ listNodes[i] ].size(); ++j )
        {
            if(mapMuts[ listCopiesUse[  listNodes[i] ] [j] ] == 0   )
            {
                ++n0;
            }
            else
            {
                ++n1;
            }
        }
        int nid = listNodes[i]->GetCellId();
        mapCellGenos[nid].first = n0;
        mapCellGenos[nid].second = n1;
    }
    
    // dump genotypes
    //for( map<int, pair<int,int> > :: iterator it = mapCellGenos.begin(); it != mapCellGenos.end(); ++it )
    //{
    //    cout << "Cell " << it->first << ": " << it->second.first << "  " << it->second.second << endl;
    //}
}

static void SimulatePointMutsOnNoRecur( CellTreeNode *pTreeRoot, const map<CellTreeNode *, vector<CellSeqCopy *> > &listCopies, map<int, pair<int,int> > &mapCellGenos )
{
    // enforce minor allele frequency
    while(true)
    {
        map<int, pair<int,int> > mapCellGenosStep;
        SimulatePointMutsOnSingleNoRecur( pTreeRoot, listCopies, mapCellGenosStep );
        
        // is this passing the allele frequency threshold test?
        int numMuts = 0;
        for(map<int,pair<int,int> > :: iterator it = mapCellGenosStep.begin(); it != mapCellGenosStep.end(); ++it)
        {
            if( it->second.second > 0 )
            {
                ++numMuts;
            }
        }
        if( numMuts >= (int)( mapCellGenosStep.size() * minAlleFreq ) )
        {
            mapCellGenos = mapCellGenosStep;
            break;
        }
    }
}


bool SimulatePointMutsOn( CellTreeNode *pTreeRoot, const map<CellTreeNode *, vector<CellSeqCopy *> > &listCopies, map<int, pair<int,int> > &mapCellGenos )
{
    // return FALSE if recurrent mut
    bool fRec = IsRecurrMut();
    if( fRec )
    {
//cout << "RECURRENT MUTATION\n";
        return SimulatePointMutsOnRecur( pTreeRoot, listCopies, mapCellGenos );
    }
    
    SimulatePointMutsOnNoRecur(pTreeRoot, listCopies, mapCellGenos);
    return true;
}

bool SimulatePointMutsOnRecur( CellTreeNode *pTreeRoot, const map<CellTreeNode *, vector<CellSeqCopy *> > &listCopies, map<int, pair<int,int> > &mapCellGenos )
{
    // recurrent mutation; for now support only binary genotypes
    map<int, pair<int,int> > mapCellGenos1, mapCellGenos2;
    SimulatePointMutsOnNoRecur( pTreeRoot, listCopies, mapCellGenos1 );
    SimulatePointMutsOnNoRecur( pTreeRoot, listCopies, mapCellGenos2 );
    
    // combine to form a recurrent site
    for( map<int,pair<int,int> > :: iterator it1 = mapCellGenos1.begin(); it1 != mapCellGenos1.end(); ++it1)
    {
        map<int,pair<int,int> > :: iterator it2 = mapCellGenos2.find(it1->first);
        if( it2 == mapCellGenos2.end() )
        {
            cout << "WRONG\n";
            exit(1);
        }
        pair<int,int> pp;
        pp.first = 2;
        pp.second = 0;
        if( (it1->second.second > 0 && it2->second.second == 0) ||  (it1->second.second == 0 && it2->second.second > 0)  )
        {
            pp.first = 1;
            pp.second = 1;
        }
        
        mapCellGenos[it1->first] = pp;
    }
    
    return false;
}

static double SimulateNoiseForGenoSingleAllele( int cell, int site, bool fAllele1, pair<int,int> &pairAlleles, bool &fDropIn )
{
    // one allele: update map of genotype prob
    int allele = 0;
    if(fAllele1)
    {
        allele = 1;
    }
    
    if( cell >=(int)listCellDropoutRates.size() )
    {
        cout << "Cell dropout rates are not initialized\n";
        exit(1);
    }
    double rateDropoutCell = listCellDropoutRates[cell];
    bool fdrop = RndEvt(rateDropoutCell);
    if( fdrop )
    {
        fDropIn = true;
        return rateDropoutCell;
    }
    bool ferr = RndEvt(rateErr);
    if( ferr )
    {
        if( allele == 0 )
        {
            ++pairAlleles.second;
        }
        else
        {
            ++pairAlleles.first;
        }
        return (1.0-rateDropout)*rateErr;
    }
    else
    {
        if( allele == 0 )
        {
            ++pairAlleles.first;
        }
        else
        {
            ++pairAlleles.second;
        }
        return (1.0-rateDropout) * (1.0-rateErr);
    }
}
static double SimulateNoiseForGenoSingleGeno( int cell, int site, const pair<int,int> &geno, pair<int,int> &genoSim, map<pair<int,int>, pair<int,int> > &listDropPositions)
{
    // assume two hapotypes for now
    int numHaps = 0;
    genoSim.first = 0;
    genoSim.second = 0;
    double prob = 1.0;
    pair<int,int> pp(cell, site);
    for(int i=0; i<geno.first; ++i)
    {
        bool fDrop0 = false;
        double p = SimulateNoiseForGenoSingleAllele( cell, site, false, genoSim, fDrop0 );
        prob *= p;
        
        if(fDrop0)
        {
            if( listDropPositions.find(pp) == listDropPositions.end() )
            {
                listDropPositions[pp].first = 0;
                listDropPositions[pp].second = 0;
            }
            if( numHaps == 0)
            {
                listDropPositions[pp].first = 1;
            }
            else
            {
                listDropPositions[pp].second = 1;
            }
cout << "Cell: " << cell << ", site: " << site << " has 0-drop\n";
        }
        ++numHaps;
    }
    for(int i=0; i<geno.second; ++i)
    {
        bool fDrop1 = false;
        double p = SimulateNoiseForGenoSingleAllele( cell, site, true, genoSim, fDrop1 );
        prob *= p;
        
        if( fDrop1)
        {
            if( listDropPositions.find(pp) == listDropPositions.end() )
            {
                listDropPositions[pp].first = 0;
                listDropPositions[pp].second = 0;
            }
            if( numHaps == 0)
            {
                listDropPositions[pp].first = 1;
            }
            else
            {
                listDropPositions[pp].second = 1;
            }
cout << "Cell: " << cell << ", site: " << site << " has 1-drop\n";
        }
        ++numHaps;
    }
    return prob;
}

static void MakeProbInRange(double &prob)
{
    if(prob < SMALL_PROB)
    {
        prob = SMALL_PROB;
    }
    if(prob> LARGE_PROB)
    {
        prob = LARGE_PROB;
    }
}

static double CalcProb0(int cell)
{
    // allele 0: may be a result of dropout or genotype error
    double rateDropoutCell = listCellDropoutRates[cell];
    double res = 1.0-rateDropoutCell-rateErr;
    MakeProbInRange(res);
    return res;
}
static double CalcProb1()
{
    // allele 1: error is only cause
    return 1.0-rateErr;
}
#if 0
static double CalcProb0Boost(double prob)
{
    // 1 droput to 0; so dec prob of being zero
    return prob/boostFac;
}
static double CalcProb1Boost(double prob)
{
    // allele 1 remain 1 (not droput)
    double probUse = prob*boostFac;
    if( probUse >= 1.0)
    {
        probUse = 1.0;
    }
    return probUse;
}
#endif

//******************************************************************************************************
// beta-binomial related

static int SampleBetaBinomial(int numTrials, int a, int b)
{
    // beta-binomail: among n trials (with initial parameters a and b)
    // approach: use URN model; assume a and b are positive
    // assume random is initialized
    int numReds=a, numBlacks=b;
    int res = 0;
    for(int i=1; i<=numTrials; ++i)
    {
        double probReds = ((double)numReds)/(numReds+numBlacks);
        bool fRed = RndEvt(probReds);
        if(fRed)
        {
            ++numReds;
            ++res;
        }
        else
        {
            ++numBlacks;
        }
    }
    return res;
}

static void TestBetaBinomial()
{
    int a=2, b=2;
    int n=10;
    
    // repeats certain number of times
    int numTests=100;
    map<int,int> listResults;
    for(int i=0; i<numTests; ++i)
    {
        int res = SampleBetaBinomial(n, a, b);
        if(listResults.find(res)==listResults.end())
        {
            listResults[res]=0;
        }
        ++listResults[res];
    }
    cout << "TestBetaBinomial: list of results: ";
    for( map<int,int> :: iterator it = listResults.begin(); it != listResults.end(); ++it )
    {
        cout << it->first << " " << it->second << endl;
    }
    cout << endl;
}


//******************************************************************************************************
// sequence reads sampling from genotypes

static int SampleReadDepth()
{
    static bool fInit = false;
    //std::random_device rd{};
    static std::mt19937 gen(time(0));
    if( fInit == false)
    {
        //
        const int RAND_SEED1 = 18973;
        gen.seed( RAND_SEED1 );
        fInit = true;
    }
    
    int depth=1;
    
    if( fBetaBinomial == false )
    {
        // use default normal distribution
        // values near the mean are the most likely
        // standard deviation affects the dispersion of generated values from the mean
        std::normal_distribution<double> d(aveReadDepth, stdReadDepth);
        
        depth = std::round(d(gen));
        if( depth < 0 )
        {
            depth = 0;
        }
//cout << "default normal distribution: ave read depth: " << aveReadDepth << ", std depth " << stdReadDepth << ", depth = "  << depth << endl;
    }
    else
    {
        int numBetaTests = 2*((int)aveReadDepth);
        depth = SampleBetaBinomial( numBetaTests, betaparama, betaparamb );
//cout << "read depth (by beta-binomial): " << depth << endl;
    }

    return depth;
}

static void SampleStrandProb(int numCopies, vector<double> &listStrandProbs)
{
    // suppose there are different copies at a site. A copy may not have all the same (say 50/50 chance)
    // this adds noise for things like dropout
    listStrandProbs.clear();
    if( numCopies == 0 )
    {
        return;
    }
    //init to all uniform
    for(int i=0; i<numCopies; ++i)
    {
        listStrandProbs.push_back(1.0/numCopies);
    }
    // add random noise
    //std::random_device rd{};
    static bool fInit = false;
    static std::mt19937 gen(time(0));
    if( fInit == false )
    {
        const int RAND_SEED2 = 258703;
        gen.seed(RAND_SEED2);
        fInit = true;
    }
    std::normal_distribution<double> d(aveStrandBias, stdStrandBias);
    double sumVal = 0.0;
    for(int i=0; i<numCopies; ++i)
    {
        double bias = d(gen);
        listStrandProbs[i] += bias;
        if( listStrandProbs[i] <= SMALL_NUM )
        {
            listStrandProbs[i] = SMALL_NUM;
        }
        sumVal += listStrandProbs[i];
    }
    // normalize to be 1.0
    for(int i=0; i<numCopies; ++i)
    {
        listStrandProbs[i] = listStrandProbs[i]/numCopies;
    }
}

static int SampleStrand(const vector<double> &listStrandProbs)
{
    // decide which copy/strand to sample
    double rndFrac = RndFrac();
    int res = (int)listStrandProbs.size()-1;
    double sum = 0.0;
    for(int i=0; i<(int)listStrandProbs.size(); ++i)
    {
        sum += listStrandProbs[i];
        if( rndFrac < sum)
        {
            res = i;
            break;
        }
    }
    return res;
}

static void SampleReadsForGenotype(int cell, int site, int geno, const map<pair<int,int>, pair<int,int> > &listDropPositions, pair<int,int> &readsCount01 )
{
    // for a genotype
    int numCopies = 2;
    //vector<double> listStrandProbs;
    //SampleStrandProb(numCopies, listStrandProbs);
    //int numReads = SampleReadDepth();
    
    // test if dropout occurs; if so, only go with allele0
    //double rndDrop = RndFrac();
    //if( rndDrop < rateDropout )
    //{
//cout << "Dropout occurs for genotype: ";
//cout << "genotype: [" << geno.first << "," << geno.second << "] \n";
        // droput
    //    readsCount01.first = numReads;
    //    readsCount01.second = 0;
    //    return;
    //}
    
    pair<int,int> pp(cell, site);
    map<pair<int,int>,pair<int,int> > :: const_iterator it = listDropPositions.find(pp);
    
    
    readsCount01.first = 0;
    readsCount01.second = 0;
    for(int strand=0; strand<numCopies; ++strand)
    {
        int numReads = SampleReadDepth();
        if( strand == 0)
        {
            if( it != listDropPositions.end() && it->second.first ==1)
            {
                // dropout
//cout << "cell: " << cell << ", site: " << site << " has a 0-droput.\n";
                continue;
            }
            
#if 0
            if( geno != 2 )
            {
                readsCount01.first += numReads;
            }
            else
            {
                readsCount01.second += numReads;
            }
#endif
//#if 0
            for( int kk = 0; kk<numReads; ++kk)
            {
                // simulate reads error
                bool fErr = RndEvt(rateErr);
                if( (fErr == false && geno != 2) || (fErr == true && geno == 2) )
                {
                    ++readsCount01.first;
                }
                else
                {
                    ++readsCount01.second;
                }
            }
//#endif
        }
        else
        {
            if( it != listDropPositions.end() && it->second.second ==1)
            {
//cout << "cell: " << cell << ", site: " << site << " has a 1-droput.\n";
                // dropout
                continue;
            }
            
//#if 0
            for( int kk = 0; kk<numReads; ++kk)
            {
                // simulate reads error
                bool fErr = RndEvt(rateErr);
                if( (fErr == false && geno == 0) || (fErr == true && geno != 0) )
                {
                    ++readsCount01.first;
                }
                else
                {
                    ++readsCount01.second;
                }
                
            }
//#endif
#if 0
            if( geno == 0)
            {
                readsCount01.first += numReads;
            }
            else
            {
                readsCount01.second += numReads;
            }
#endif
        }
        
#if 0
        int strand = SampleStrand(listStrandProbs);
        if( strand == 0)
        {
            if( geno != 2 )
            {
                ++readsCount01.first;
            }
            else
            {
                ++readsCount01.second;
            }
        }
        else
        {
            if( geno == 0)
            {
                ++readsCount01.first;
            }
            else
            {
                ++readsCount01.second;
            }
        }
#endif
    }
    
//cout << "genotype: " << geno << ": ";
//cout << "Read count: " << readsCount01.first << "," << readsCount01.second << endl;
}

static void SampleReadsForGenotype2(int cell, int site, int nAllele0, int nAllele1, pair<int,int> &readsCount01 )
{
//cout << "SampleReadsForGenotype2: cell " << cell << ", site << " << site << " (" << nAllele0 << "," << nAllele1 << ")" << "Cell dropout rate: " << listCellDropoutRates[cell] << endl;
    // for a genotype
    //vector<double> listStrandProbs;
    //SampleStrandProb(numCopies, listStrandProbs);
    
    // create a generic list of alleles
    vector<int> listAlleles;
    for(int i=0; i<nAllele0; ++i)
    {
        listAlleles.push_back(0);
    }
    for(int i=0; i<nAllele1; ++i)
    {
        listAlleles.push_back(1);
    }
    //cout << "List of alleles: ";
    //for(unsigned int i=0; i<listAlleles.size(); ++i)
    //{
    //    cout << listAlleles[i];
    //}
    //cout << endl;
    // test if dropout occurs; if so, only go with allele0
    //double rndDrop = RndFrac();
    //if( rndDrop < rateDropout )
    //{
//cout << "Dropout occurs for genotype: ";
//cout << "genotype: [" << geno.first << "," << geno.second << "] \n";
        // droput
    //    readsCount01.first = numReads;
    //    readsCount01.second = 0;
    //    return;
    //}
    
    readsCount01.first = 0;
    readsCount01.second = 0;
    for(unsigned int i=0; i<listAlleles.size(); ++i)
    {
        // first decide if we drop out this allele altogether
        bool fDrop = RndEvt(  listCellDropoutRates[cell] );
        if( fDrop )
        {
//cout << "One allele dropout..\n";
            continue;
        }
        
        int allele = listAlleles[i];
        int numReads = SampleReadDepth();
//cout << "Num reads: " << numReads << endl;
        
        // now consider each read separately to accomondate reads error
        for(int kk=0; kk<numReads; ++kk)
        {
            // simulate reads error
            bool fErr = RndEvt(rateErr);
            int rd = allele;
            if( fErr )
            {
                rd = 1 - rd;
            }
//cout << "read: " << rd << endl;
            if( rd == 0 )
            {
                ++readsCount01.first;
            }
            else
            {
                ++readsCount01.second;
            }
        }
        
    }
//cout << "Read count: " << readsCount01.first << "," << readsCount01.second << endl;
}



static double CalcReadAlleleForGeno(int allele, int g)
{
    // allele = 0 or 1, g=0,1,2
    // use errRate as the reads error rate, also drop rate matters too
    if( g == 0 )
    {
        if( allele == 0 )
        {
            return 1.0 - rateErr;
        }
        else
        {
            return rateErr;
        }
    }
    else if(g == 2)
    {
        if( allele == 1 )
        {
            return 1.0 - rateErr;
        }
        else
        {
            return rateErr;
        }
    }
    else
    {
        return 0.5;
    }
}

static double CalcReadAlleleForGenoWithOneDrop0(int n0, int n1)
{
    double p0 = 1.0 - rateErr;
    double p1 = rateErr;
    return pow(p0, n0) * pow(p1, n1);
}
static double CalcReadAlleleForGenoWithOneDrop2(int n0, int n1)
{
    double p1 = 1.0 - rateErr;
    double p0 = rateErr;
    return pow(p0, n0) * pow(p1, n1);
}

static double CalcReadAlleleForGenoWithOneDrop(int g, int n0, int n1)
{
    // droput can occur at either allele
    if( g== 0 )
    {
        return CalcReadAlleleForGenoWithOneDrop0(n0, n1);
    }
    else if(g == 2)
    {
        return CalcReadAlleleForGenoWithOneDrop2(n0, n1);
    }
    else
    {
        // assume dropout occurs w/ 50/50 chance at one allele
        double p00 = CalcReadAlleleForGenoWithOneDrop0(n0, n1);
        double p11 = CalcReadAlleleForGenoWithOneDrop2(n0, n1);
        return (p00+p11)*0.5;
    }
}

static void CalcGenotypeProbOfReads(const pair<int,int> &readsCount01, double allele0Freq, vector<double> &listGenoProbs)
{
#if 0
//cout << "rateDrop: " << rateDropout << ", Allele0 freq: " << allele0Freq << ", read count: [" << readsCount01.first << "," << readsCount01.second << "]" << endl;
    // allele0Freq: freq of allele 0 at this site
    // YW: assume 0/1 genotypes
    double p00 = allele0Freq, p10 = 1.0-allele0Freq;
    double pa0g0 = CalcReadAlleleForGeno(0, 0);
    double pa0g1 = CalcReadAlleleForGeno(0, 1);
    double pa0g2 = CalcReadAlleleForGeno(0, 2);
    double pa1g0 = CalcReadAlleleForGeno(1, 0);
    double pa1g1 = CalcReadAlleleForGeno(1, 1);
    double pa1g2 = CalcReadAlleleForGeno(1, 2);
    double pa0g12 = pa0g1+pa0g2, pa1g12 = pa1g1+pa1g2;
    double p0 = pow(pa0g0, readsCount01.first) * pow(pa1g0, readsCount01.second);
    double p1 = pow(pa0g12, readsCount01.first) * pow(pa1g12, readsCount01.second);
    // now considering the droput
    double pdrop = pow(1.0-rateErr,readsCount01.first)*pow(rateErr,readsCount01.second);
    //cout << "pdrop: " << pdrop << endl;
    p0 = (1.0-rateDropout)*p0 + rateDropout*pdrop;
    p1 = (1.0-rateDropout)*p1 + rateDropout*pdrop;
    p0 *= p00;
    p1 *= p10;
//cout << "With droput: p0: " << p0 << ", p1: " << p1 << endl;
    
    // normalize
    double sum = p0+p1;
    p0 = p0/sum;
    p1 = p1/sum;
    
    listGenoProbs.clear();
    listGenoProbs.push_back(p0);
    listGenoProbs.push_back(p1);
    listGenoProbs.push_back(0.0);
    
    // make sure the prob is not too small
    MakeProbInRange( listGenoProbs[0] );
    MakeProbInRange( listGenoProbs[1] );
    MakeProbInRange( listGenoProbs[2] );
#endif
    
//#if 0
    //double p0 = rateDropout + (1.0-rateDropout) * allele0Freq*allele0Freq, p1 = (1.0-rateDropout) * 2*allele0Freq*(1.0-allele0Freq), p2 = (1.0-rateDropout) * (1.0-allele0Freq)*(1.0-allele0Freq);
    double p00 = allele0Freq*allele0Freq, p10 = 2*allele0Freq*(1.0-allele0Freq), p20 = (1.0-allele0Freq)*(1.0-allele0Freq);
    double p0 = 1.0, p1 = 1.0, p2 = 1.0;
    double pa0g0 = CalcReadAlleleForGeno(0, 0);
    double pa0g1 = CalcReadAlleleForGeno(0, 1);
    double pa0g2 = CalcReadAlleleForGeno(0, 2);
    double pa1g0 = CalcReadAlleleForGeno(1, 0);
    double pa1g1 = CalcReadAlleleForGeno(1, 1);
    double pa1g2 = CalcReadAlleleForGeno(1, 2);
//cout << "p00:" << p00 << ", p10:" << p10 << ",p20:" << p20 << ",pa0g0: " << pa0g0 << ", pa0g1: " << pa0g1 << ", pa0g2: " << pa0g2 << ", pa1g0:" << pa1g0 << ", pa1g1:" << pa1g1 << ", pa1g2:" << pa1g2 << endl;
    p0 = pow(pa0g0, readsCount01.first) * pow(pa1g0, readsCount01.second) * p00;
    p1 = pow(pa0g1, readsCount01.first) * pow(pa1g1, readsCount01.second) * p10;
    p2 = pow(pa0g2, readsCount01.first) * pow(pa1g2, readsCount01.second) * p20;
//cout << "Before normalize: p0: " << p0 << ", p1: " << p1 << ", p2: " << p2 << endl;
    double sum = p0+p1+p2;
    p0 = p0/sum;
    p1 = p1/sum;
    p2 = p2/sum;
    
#if 0
    // add droput
    double p0all = rateDropout + (1.0-rateDropout)*p0;
    double p1all = (1.0-rateDropout)*p1;
    double p2all = (1.0-rateDropout)*p2;
#endif
    
    double p0d = CalcReadAlleleForGenoWithOneDrop(0, readsCount01.first, readsCount01.second);
    double p1d = CalcReadAlleleForGenoWithOneDrop(1, readsCount01.first, readsCount01.second);
    double p2d = CalcReadAlleleForGenoWithOneDrop(2, readsCount01.first, readsCount01.second);
    double p0all = rateDropout*p0d + (1.0-rateDropout)*p0;
    double p1all = rateDropout*p1d + (1.0-rateDropout)*p1;
    double p2all = rateDropout*p2d + (1.0-rateDropout)*p2;
    
    listGenoProbs.clear();
    listGenoProbs.push_back(p0all);
    listGenoProbs.push_back(p1all);
    listGenoProbs.push_back(p2all);
    
    // make sure the prob is not too small
    MakeProbInRange( listGenoProbs[0] );
    MakeProbInRange( listGenoProbs[1] );
    MakeProbInRange( listGenoProbs[2] );
//#endif
    
//cout <<"List of simulated genotype probability (after normalize): " << listGenoProbs[0] << "," << listGenoProbs[1] << "," << listGenoProbs[2] << endl;
}


static double CalcProbReadDepthNormDist(double numReads, double aveReadDepthSite, double stdReadDepthSite )
{
    // assume normal distribution
    double diff = numReads - aveReadDepthSite;
    return exp( -1.0*diff*diff/(2.0*stdReadDepthSite*stdReadDepthSite) );
}
static double GammaFunc(int n)
{
    double res=1.0;
    for(int i=1; i<=n-1; ++i)
    {
        res *= i;
    }
    return res;
}
static double CalcProbReadDepthBetaBinom(int numReads, int aveReadDepthSite)
{
    int n=aveReadDepthSite*(betaparama+betaparamb)/betaparama;
    double res = (GammaFunc(n+1)/( GammaFunc(numReads+1)*GammaFunc(n+1-numReads) ) ) * (GammaFunc(numReads+betaparama) * GammaFunc(n-numReads+betaparamb)/GammaFunc(n+betaparama+betaparamb) ) * GammaFunc(betaparama+betaparamb) / ( GammaFunc(betaparama)*GammaFunc(betaparamb) );
    return res;
}
static double CalcGenotypeProbForGenoSingleStrand( const pair<int,int> &readsCount01,  double aveReadDepthSite, double stdReadDepthSite, int geno )
{
    
    // single strand: no droput,
    // if the geno is different from the read, must be errors
    double probErr = 1.0;
    if( geno == 0 )
    {
        probErr = pow(rateErr, readsCount01.second);
    }
    else
    {
        probErr = pow(rateErr, readsCount01.first);
    }
    int readNum = readsCount01.first+readsCount01.second;
    double probDepth = 1.0;
    if( fBetaBinomial == false )
    {
        probDepth = CalcProbReadDepthNormDist(readNum, aveReadDepthSite, stdReadDepthSite);
    }
    else
    {
        int aveReadDepthSiteUse = (double)aveReadDepthSite;
        if(aveReadDepthSiteUse <= 0 )
        {
            aveReadDepthSiteUse = 1;
        }
        probDepth = CalcProbReadDepthBetaBinom(readNum, aveReadDepthSiteUse);
    }
    double res = probDepth * probErr;
//cout << "readNum: " << readNum << ", probDepth: " << probDepth << ", probErr: " << probErr << endl;
//cout << "CalcGenotypeProbForGenoSingleStrand: read: " << readsCount01.first << "," << readsCount01.second << ", geno: " << geno << ", aveRead: " << aveReadDepthSite << ", stdRead: " << stdReadDepthSite << ", prob=" << res << endl;
    return res;
}
static double CalcGenotypeProbForGenoDrop( const pair<int,int> &readsCount01,  double aveReadDepthSite, double stdReadDepthSite, int geno0, int geno1, bool fdrop0, bool fdrop1 )
{
    if( fdrop0 == true && fdrop1 == true )
    {
        // all errors
        double res = pow(rateErr, readsCount01.first+readsCount01.second);
        return res;
    }
    else if( fdrop0 == true && fdrop1 == false )
    {
        return CalcGenotypeProbForGenoSingleStrand( readsCount01, aveReadDepthSite, stdReadDepthSite, geno1 );
    }
    else if( fdrop0 == false && fdrop1 == true )
    {
        return CalcGenotypeProbForGenoSingleStrand( readsCount01, aveReadDepthSite, stdReadDepthSite, geno0 );;
    }
    else
    {
        // consider cases: if homozygote, then take half
        if( geno0 == geno1)
        {
            pair<int,int> pp1(readsCount01.first/2, readsCount01.second/2);
            pair<int,int> pp2(readsCount01.first-pp1.first, readsCount01.second -pp1.second);
            double p1=CalcGenotypeProbForGenoSingleStrand( pp1, aveReadDepthSite, stdReadDepthSite, geno0 );
            double p2=CalcGenotypeProbForGenoSingleStrand( pp2, aveReadDepthSite, stdReadDepthSite, geno1 );
            return p1*p2;
        }
        else
        {
            // heterozygote: just use the read count
            pair<int,int> pp1(readsCount01.first, 0), pp2(0, readsCount01.second);
            double p1=CalcGenotypeProbForGenoSingleStrand( pp1, aveReadDepthSite, stdReadDepthSite, 0 );
            double p2=CalcGenotypeProbForGenoSingleStrand( pp2, aveReadDepthSite, stdReadDepthSite, 1 );
            return p1*p2;
        }
    }
}

static void CalcGenotypeProbOfReadsNew( int cell, double aveReadDepthSite, double stdReadDepthSite, const pair<int,int> &readsCount01, double allele0Freq, int numInfDrops, vector<double> &listGenoProbs)
{
//cout << "CalcGenotypeProbOfReadsNew: cell: " << cell << ", aveReadDepthSite: " << aveReadDepthSite << ", stdReadDepthSite: " << stdReadDepthSite << ", read count: " << readsCount01.first << ", " << readsCount01.second  << ", allele0Freq: " << allele0Freq << ", numInfDrops: " << numInfDrops << endl;
    
    
#if 0  // add droput imputation
    // assume normal distribution of read counts; calculate prob of 0
    double p00d00 = CalcGenotypeProbForGenoDrop( readsCount01, aveReadDepthSite, stdReadDepthSite, 0, 0, false, false );
    double p00d01 = CalcGenotypeProbForGenoDrop( readsCount01, aveReadDepthSite, stdReadDepthSite, 0, 0, false, true );
    double p00d10 = CalcGenotypeProbForGenoDrop( readsCount01, aveReadDepthSite, stdReadDepthSite, 0, 0, true, false );
    double p00d11 = CalcGenotypeProbForGenoDrop( readsCount01, aveReadDepthSite, stdReadDepthSite, 0, 0, true, true );
    double p01d00 = CalcGenotypeProbForGenoDrop( readsCount01, aveReadDepthSite, stdReadDepthSite, 0, 1, false, false );
    double p01d01 = CalcGenotypeProbForGenoDrop( readsCount01, aveReadDepthSite, stdReadDepthSite, 0, 1, false, true );
    double p01d10 = CalcGenotypeProbForGenoDrop( readsCount01, aveReadDepthSite, stdReadDepthSite, 0, 1, true, false );
    double p01d11 = CalcGenotypeProbForGenoDrop( readsCount01, aveReadDepthSite, stdReadDepthSite, 0, 1, true, true );
    double p11d00 = CalcGenotypeProbForGenoDrop( readsCount01, aveReadDepthSite, stdReadDepthSite, 1, 1, false, false );
    double p11d01 = CalcGenotypeProbForGenoDrop( readsCount01, aveReadDepthSite, stdReadDepthSite, 1, 1, false, true );
    double p11d10 = CalcGenotypeProbForGenoDrop( readsCount01, aveReadDepthSite, stdReadDepthSite, 1, 1, true, false );
    double p11d11 = CalcGenotypeProbForGenoDrop( readsCount01, aveReadDepthSite, stdReadDepthSite, 1, 1, true, true );
    //cout << "p00d00:" << p00d00 << ", p00d01:" << p00d01 << ", p00d10:" << p00d10 << ", p00d11: " << p00d11 << endl;
    //cout << "p01d00:" << p01d00 << ", p01d01:" << p01d01 << ", p01d10:" << p01d10 << ", p01d11: " << p01d11 << endl;
    //cout << "p11d00:" << p11d00 << ", p11d01:" << p11d01 << ", p11d10:" << p11d10 << ", p11d11: " << p11d11 << endl;
    double pprior00 = allele0Freq*allele0Freq;
    double pprior01 = 2.0*(1.0-allele0Freq)*allele0Freq;
    double pprior11 = (1.0-allele0Freq)*(1.0-allele0Freq);
    //cout << "facNorm: " << facNorm << endl;
    listGenoProbs.clear();
    double p00 = p00d00*pprior00;
    if( numInfDrops == 1 )
    {
        p00 = (p00d01 + p00d10)*pprior00;
    }
    else if( numInfDrops == 2 )
    {
        p00 = p00d11*pprior00;
    }
    double p01 = p01d00*pprior01;
    if( numInfDrops == 1 )
    {
        p01 = (p01d01 + p01d10)*pprior01;
    }
    else if( numInfDrops == 2 )
    {
        p01 = p01d11*pprior01;
    }
    double p11 = p11d00*pprior11;
    if( numInfDrops == 1 )
    {
        p11 = (p11d01 + p11d10)*pprior11;
    }
    else if( numInfDrops == 2 )
    {
        p11 = p11d11*pprior11;
    }
    double facNorm = p00+p01+p11;
    listGenoProbs.push_back(p00/facNorm);
    listGenoProbs.push_back(p01/facNorm);
    listGenoProbs.push_back(p11/facNorm);
    MakeProbInRange( listGenoProbs[0] );
    MakeProbInRange( listGenoProbs[1] );
    MakeProbInRange( listGenoProbs[2] );
//cout << "probability: " << listGenoProbs[0] << ", " << listGenoProbs[1] << ", " << listGenoProbs[2] << endl;
    
#endif
    
    
//#if 0
    // assume normal distribution of read counts; calculate prob of 0
    double p00d00 = CalcGenotypeProbForGenoDrop( readsCount01, aveReadDepthSite, stdReadDepthSite, 0, 0, false, false );
    double p00d01 = CalcGenotypeProbForGenoDrop( readsCount01, aveReadDepthSite, stdReadDepthSite, 0, 0, false, true );
    double p00d10 = CalcGenotypeProbForGenoDrop( readsCount01, aveReadDepthSite, stdReadDepthSite, 0, 0, true, false );
    double p00d11 = CalcGenotypeProbForGenoDrop( readsCount01, aveReadDepthSite, stdReadDepthSite, 0, 0, true, true );
    double p01d00 = CalcGenotypeProbForGenoDrop( readsCount01, aveReadDepthSite, stdReadDepthSite, 0, 1, false, false );
    double p01d01 = CalcGenotypeProbForGenoDrop( readsCount01, aveReadDepthSite, stdReadDepthSite, 0, 1, false, true );
    double p01d10 = CalcGenotypeProbForGenoDrop( readsCount01, aveReadDepthSite, stdReadDepthSite, 0, 1, true, false );
    double p01d11 = CalcGenotypeProbForGenoDrop( readsCount01, aveReadDepthSite, stdReadDepthSite, 0, 1, true, true );
    double p11d00 = CalcGenotypeProbForGenoDrop( readsCount01, aveReadDepthSite, stdReadDepthSite, 1, 1, false, false );
    double p11d01 = CalcGenotypeProbForGenoDrop( readsCount01, aveReadDepthSite, stdReadDepthSite, 1, 1, false, true );
    double p11d10 = CalcGenotypeProbForGenoDrop( readsCount01, aveReadDepthSite, stdReadDepthSite, 1, 1, true, false );
    double p11d11 = CalcGenotypeProbForGenoDrop( readsCount01, aveReadDepthSite, stdReadDepthSite, 1, 1, true, true );
//cout << "p00d00:" << p00d00 << ", p00d01:" << p00d01 << ", p00d10:" << p00d10 << ", p00d11: " << p00d11 << endl;
//cout << "p01d00:" << p01d00 << ", p01d01:" << p01d01 << ", p01d10:" << p01d10 << ", p01d11: " << p01d11 << endl;
//cout << "p11d00:" << p11d00 << ", p11d01:" << p11d01 << ", p11d10:" << p11d10 << ", p11d11: " << p11d11 << endl;
    double rateDropoutCell = listCellDropoutRates[cell];
//cout << "rateDropoutCell: " << rateDropoutCell << endl;
    double p00d00prior = allele0Freq*allele0Freq*(1.0-rateDropoutCell)*(1.0-rateDropoutCell);
    double p00d01prior = allele0Freq*allele0Freq*rateDropoutCell*(1.0-rateDropoutCell);
    double p00d11prior = allele0Freq*allele0Freq*rateDropoutCell*rateDropoutCell;
    double p01d00prior = 2.0*(1.0-allele0Freq)*allele0Freq*(1.0-rateDropoutCell)*(1.0-rateDropoutCell);
    double p01d01prior = 2.0*(1.0-allele0Freq)*allele0Freq*rateDropoutCell*(1.0-rateDropoutCell);
    double p01d11prior = 2.0*(1.0-allele0Freq)*allele0Freq*rateDropoutCell*rateDropoutCell;
    double p11d00prior = (1.0-allele0Freq)*(1.0-allele0Freq)*(1.0-rateDropoutCell)*(1.0-rateDropoutCell);
    double p11d01prior = (1.0-allele0Freq)*(1.0-allele0Freq)*rateDropoutCell*(1.0-rateDropoutCell);
    double p11d11prior = (1.0-allele0Freq)*(1.0-allele0Freq)*rateDropoutCell*rateDropoutCell;
//cout << "p00d00prior:" << p00d00prior << ", p00d01prior:" << p00d01prior << ", p00d11prior:" << p00d11prior << endl;
//cout << "p01d00prior:" << p01d00prior << ", p01d01prior:" << p00d01prior << ", p01d11prior:" << p01d11prior << endl;
//cout << "p11d00prior:" << p11d00prior << ", p11d01prior:" << p11d01prior << ", p11d11prior:" << p11d11prior << endl;
    double facNorm = p00d00*p00d00prior + p00d01*p00d01prior + p00d10*p00d01prior + p00d11*p00d11prior + p01d00*p01d00prior + p01d01*p01d01prior + p01d10*p01d01prior + p01d11*p01d11prior + p11d00*p11d00prior + p11d01*p11d01prior + p11d10*p11d01prior + p11d11*p11d11prior;
//cout << "facNorm: " << facNorm << endl;
    listGenoProbs.clear();
    double p00 = (p00d00*p00d00prior + p00d01*p00d01prior + p00d10*p00d01prior + p00d11*p00d11prior )/ facNorm;
    double p01 = ( p01d00*p01d00prior + p01d01*p01d01prior + p01d10*p01d01prior + p01d11*p01d11prior )/facNorm;
    double p11 = ( p11d00*p11d00prior + p11d01*p11d01prior + p11d10*p11d01prior + p11d11*p11d11prior )/facNorm;
    listGenoProbs.push_back(p00);
    listGenoProbs.push_back(p01);
    listGenoProbs.push_back(p11);
    MakeProbInRange( listGenoProbs[0] );
    MakeProbInRange( listGenoProbs[1] );
    MakeProbInRange( listGenoProbs[2] );
//cout << "probability: " << listGenoProbs[0] << ", " << listGenoProbs[1] << ", " << listGenoProbs[2] << endl;
//#endif
}

static void ConvGenoProbToOututProb(const vector<double> &listGenoProbs, bool fBinary, vector<double> &listOutProbs)
{
    listOutProbs.clear();
    if( fBinary == false )
    {
        for(int i=0; i<(int)listGenoProbs.size()-1; ++i)
        {
            listOutProbs.push_back(listGenoProbs[i]);
        }
    }
    else
    {
        // only output the zero prob
        listOutProbs.push_back(listGenoProbs[0]);
    }
}

static int CalcAllele0FreqForGenos(const vector<int> &listGenos)
{
    //
    int res = (int)(listGenos.size())*2;
    //if( fBinary == false )
    //{
    //    res *= 2;
    //}
    for( int i=0; i<(int)listGenos.size(); ++i )
    {
        res -= listGenos[i];
    }
    return res;
}
static int CountTotAlleleForGenos(const map<int,pair<int,int> > &mapGenosOfCells)
{
    //
    int res = 0;
    
    //if( fBinary )
    //{
    //    res = mapGenosOfCells.size();
    //}
    //else
    //{
    for( map<int,pair<int,int> > :: const_iterator it = mapGenosOfCells.begin(); it != mapGenosOfCells.end(); ++it )
    {
        res += it->second.first;
        res += it->second.second;
    }
    // }
    return res;
}
static int CallGenoWithProb(const vector<double> &listProbs)
{
    int res = -1;
    double probMax = -1.0;
    double probSum = 0.0;
    for(int i=0; i<(int)listProbs.size(); ++i)
    {
        if( probMax < listProbs[i])
        {
            res = i;
            probMax = listProbs[i];
        }
        probSum += listProbs[i];
    }
    double probLast = 1.0-probSum;
    if( probLast < 0.0 )
    {
        probLast = 0.0;
//cout << "probSum: " << probSum << endl;
//        cout << "FATAL ERROR: genotype probability should be between 0 and 1, and all probability should sum to 1.0 for each genotype.\n";
//        exit(1);
    }
    if( probMax < probLast )
    {
        res = listProbs.size();
    }
    return res;
}
static void ConvReadBasedProbToGeno( const vector<vector<vector<double> > > &listReadBasedProb, vector<vector<int> > &listGenos)
{
    listGenos.resize( listReadBasedProb.size() );
    for(int i=0; i<(int)listReadBasedProb.size(); ++i)
    {
        listGenos[i].resize( listReadBasedProb[i].size() );
        for(int j=0; j<(int)listReadBasedProb[i].size(); ++j)
        {
            int geno = CallGenoWithProb( listReadBasedProb[i][j] );
            listGenos[i][j] = geno;
        }
    }
}

static void AdjustReadCountForDrop(int cell, int site, const map<pair<int,int>, pair<int,int> > &setDropPos, pair<int,int> &readCount)
{
cout << "At cell " << cell << ", site: " << site << ", read count: " << readCount.first << " " << readCount.second << endl;
    pair<int,int> pp(cell, site);
    if( setDropPos.find(pp) != setDropPos.end() )
    {
        map<pair<int,int>,pair<int,int> > :: const_iterator it = setDropPos.find(pp);
        // make read count of 1 to be 0 if dropout
//cout << "Setting read count to be 0 at site " << site << ", cell " << cell << endl;
        if( it->second.first == 1 )
        {
cout << "Setting read count of 0-allele to be 0 at site " << site << ", cell " << cell << endl;
            readCount.first = 0;
        }
        if( it->second.second == 1 )
        {
cout << "Setting read count of 1-allele to be 0 at site " << site << ", cell " << cell << endl;
            readCount.second = 0;
        }
    }
}

// analyze reads counts
static void AnalyzeReadsAtSites( const vector<pair<int,int> > &listReadCnts, double &fracAllele0, double &aveReadDepthSite, double &stdReadDepthSite)
{
    //
    int numNonDoubleDrops = 0;
    int numZeros = 0, numAlleleTot = 0;
    double sumReadDepthSingleStrand = 0.0;
    vector<double> listReadSingle;
    for(int i=0; i<(int)listReadCnts.size(); ++i)
    {
        if( listReadCnts[i].first+listReadCnts[i].second>0)
        {
            int numReads = listReadCnts[i].first+listReadCnts[i].second;
            double rateDropoutCell = listCellDropoutRates[i];
            double readSingleExp = (1.0-rateDropoutCell)*numReads/2 + rateDropoutCell*numReads;
            sumReadDepthSingleStrand += readSingleExp;
            listReadSingle.push_back(readSingleExp);
            ++numNonDoubleDrops;
        }
        numZeros += listReadCnts[i].first;
        numAlleleTot += listReadCnts[i].first + listReadCnts[i].second;
    }
    fracAllele0 = ((double)numZeros+1)/((double)numAlleleTot+1);
    aveReadDepthSite = (sumReadDepthSingleStrand+1.0)/(numNonDoubleDrops+1.0);
    // std
    double sum = 0.0;
    for(int i=0; i<(int)listReadSingle.size(); ++i)
    {
        double diff = listReadSingle[i] - aveReadDepthSite;
        sum += diff*diff;
    }
    stdReadDepthSite = sqrt((sum+1.0)/(listReadSingle.size()+1.0));
}

static void InfDroputsByReadsAtSite( const vector<pair<int,int> > &listReadCnts, vector<int> &listNumDropouts)
{
    // based on the read counts on the genotypes of this site, infer the number of dropouts
    // approach: if read counts are zero, then double droput; otherwise, find the bipartiton
    vector<int> listTotNonZeroReadCounts;
    for(int i=0; i<(int)listReadCnts.size(); ++i)
    {
        int rc = listReadCnts[i].first + listReadCnts[i].second;
        if( rc > 0)
        {
            listTotNonZeroReadCounts.push_back(rc);
        }
    }
    int bound0and1 = FindTwoClusterBoundary(listTotNonZeroReadCounts);
    listNumDropouts.clear();
    for(int i=0; i<(int)listReadCnts.size(); ++i)
    {
        int rc = listReadCnts[i].first + listReadCnts[i].second;
        int ndrops = 0;
        if( rc == 0)
        {
            ndrops = 2;
        }
        else if( rc <= bound0and1)
        {
            ndrops = 1;
        }
        listNumDropouts.push_back(ndrops);
    }
//cout << "bound0and1: " << bound0and1 << ", and number of inferred dropouts: " << endl;
//cout << "InfDroputsByReadsAtSite: listReadCounts: ";
//for(int i=0; i<(int)listReadCnts.size(); ++i)
//{
//cout << "(" << listReadCnts[i].first << "," << listReadCnts[i].second << ") d: " << listNumDropouts[i] << " ";
//}
//cout << endl;
}

static void ChooseMissingPos(int numCells, int numSites, set<pair<int,int> > &setMissingPos)
{
    //
    int numMissings = numCells*numSites*fracMissing;
    for(int i=0; i<numMissings; ++i)
    {
        int rndCell = RndInt(0, numCells-1);
        int rndSite = RndInt(0, numSites-1);
        pair<int,int> pp(rndCell, rndSite);
        setMissingPos.insert(pp);
    }
}
static bool IsMissingAt(int cell, int site, const set<pair<int,int> > &setMissingPos)
{
    pair<int,int> pp(cell,site);
    return setMissingPos.find(pp) != setMissingPos.end();
}
static double GetMissingProb()
{
    if( fBinary)
    {
        return 0.5;
    }
    else
    {
        return 0.33;
    }
}

//******************************************************************************************************

void SimulateIndependentSitesNoisyGeno(CellTreeNode *pRoot, int numSites, int numVariantsPerSite)
{
    //const int SEED_RND = 0xF1CD1321;
    InitRand(seedRndUser);

// test
//TestBetaBinomial();
//exit(1);
    
    map<int, vector< pair<int,int> > > mapListReadCnts;
    vector< vector<double> > listGenoProbs;
    vector<vector<int> > listGenoSimPreDbl;
    vector<vector<int> > listGenoSim;
    vector<vector<pair<int,int> > > listGenoGenSim;  // [site][cell]: pair = (# of 0 alleles, # of 1 alleles)
    vector<vector<int> > listGenoNoise;
    vector<vector<bool> > listGenoDrop;
    vector<vector<pair<int,int> > > listReadCountsAtSites;      // [site, cell] = <read count 0, read count 1>
    vector<vector<pair<int,int> > > listReadCountsAtSites2;      // [site, cell] = <read count 0, read count 1>;  corrected @Oct 10, 2025
    vector<vector< vector<double> > > listReadBasedProb;        // [site,cell] = vec of genotypes of the reads <prob of g0, prob of g1..>
    vector< map<int,pair<int,int> > > listGenoSimCopiesPreDoubles;
    set<pair<int,int> > setMissingPos;
    
    // first simulate genotypes
    int indexSNV = 1;
    for(int s=0; s<numSites; ++s)
    {
//cout << "Site: " << s+1 << endl;
        map<CellTreeNode *, vector<CellSeqCopy *> > listCopies;
        SimulateCNVOnTree( pRoot, listCopies );
        
        for(int i=0; i<numVariantsPerSite; ++i)
        {
            map<int, pair<int,int> > mapCellGenos;
            bool fNoRec = SimulatePointMutsOn( pRoot, listCopies, mapCellGenos );
            
            if( fNoRec == false )
            {
                cout << "SNV " << indexSNV << ": recurrent mutation\n";
            }
            
            listGenoSimCopiesPreDoubles.push_back(mapCellGenos);
            
            ++indexSNV;
        }
        
        map<int,int> mapCellCopyNums;
        for( map<CellTreeNode *, vector<CellSeqCopy *> > :: iterator it = listCopies.begin(); it != listCopies.end(); ++it )
        {
            if( it->first->IsLeaf() )
            {
                mapCellCopyNums[ it->first->GetCellId() ] = it->second.size();
                //cout << "Cell " << it->first->GetCellId() << ": number of copies is " << it->second.size() << endl;
            }
        }
        for(map<int,int> :: iterator it = mapCellCopyNums.begin(); it != mapCellCopyNums.end(); ++it)
        {
            cout << "Site " << s << " cell: " << it->first << " number of copies = " << it->second << endl;
        }
        
        // cleanup
        for( map<CellTreeNode *, vector<CellSeqCopy *> > :: iterator it = listCopies.begin(); it != listCopies.end(); ++it )
        {
            for(int i=0; i<(int)it->second.size(); ++i)
            {
                delete it->second[i];
            }
        }
    }
    
    // print out true simulated genotypes before converting...
    cout << "LIST OF TRUE SIMULATED GENOTYPES: \n";
    for(int i=0; i<(int)listGenoSimCopiesPreDoubles.size(); ++i)
    {
        cout << "***SITE " << i << endl;
        for(map<int,pair<int,int> > :: iterator it = listGenoSimCopiesPreDoubles[i].begin(); it != listGenoSimCopiesPreDoubles[i].end(); ++it)
        {
            cout << "Cell " << it->first << ": Genotype = (" << it->second.first << "," << it->second.second << ") \n";
        }
    }
    
    
    // find out the list of true simulated genotypes prior to doublet merging
    for( int i=0; i<(int)listGenoSimCopiesPreDoubles.size(); ++i )
    {
        vector<int> vecGenoSim;
        
        // simulate reads
        map<int,pair<int,int> > &mapCellGenosUse = listGenoSimCopiesPreDoubles[i];
        for(map<int,pair<int,int> > :: iterator it = mapCellGenosUse.begin(); it != mapCellGenosUse.end(); ++it)
        {
            int hapSim = 0;
            if( it->second.second >= 1 )
            {
                hapSim = 1;
            }
            vecGenoSim.push_back(hapSim);
        }
        listGenoSimPreDbl.push_back( vecGenoSim );
    }
    
    
    vector< map<int,pair<int,int> > > listGenoSimCopiesPostDoubles = listGenoSimCopiesPreDoubles;
    int numTotCells = pRoot->GetNumLeaves();
    int numCellsOut = numTotCells - numDoublet;
    cout << "Number of cells to output: " << numCellsOut << endl;
    if( numCellsOut <= 0 )
    {
        cout << "Number of doublets is too large\n";
        exit(1);
    }
    
    if( numDoublet > 0)
    {
        // determine which one is merged to where
        map<int,int> mapCellDoublets;
        for(int h=numCellsOut+1; h<=numTotCells; ++h)
        {
            int rndCell = RndInt( 1, numCellsOut );
            mapCellDoublets[h] = rndCell;
            cout << "**Doublet cell " << h << " is merged with cell " << rndCell << endl;
        }
        
        // now deal with doublets
        // first merge the remaining ones
        for(int i=0; i<(int)listGenoSimCopiesPreDoubles.size(); ++i)
        {
            for( int h=numCellsOut+1; h<=numTotCells; ++h )
            {
                int cellDest = mapCellDoublets[h];
                listGenoSimCopiesPostDoubles[i][cellDest].first += listGenoSimCopiesPostDoubles[i][h].first;
                listGenoSimCopiesPostDoubles[i][cellDest].second += listGenoSimCopiesPostDoubles[i][h].second;
            }
        }
        
        // erase all the doublet cells
        for(int i=0; i<(int)listGenoSimCopiesPreDoubles.size(); ++i)
        {
            for(int h=numCellsOut+1; h<=numTotCells; ++h)
            {
                listGenoSimCopiesPostDoubles[i].erase( h );
            }
        }
    }
    
    // initialize droput vector
    if( listGenoSimCopiesPostDoubles.size() == 0 )
    {
        cout << "WRONG: nothing to simulate from\n";
        exit(1);
    }
    listCellDropoutRates.clear();
    for(map<int,pair<int,int> > :: iterator it = listGenoSimCopiesPostDoubles[0].begin(); it != listGenoSimCopiesPostDoubles[0].end(); ++it)
    {
        //int cell = it->first;
        double rateDropoutCell = SampleDropoutRateForCell();
cout << "Cell " << it->first-1 << ": droput rate: " << rateDropoutCell << endl;
        listCellDropoutRates.push_back(rateDropoutCell);
    }
    
    // HACK:
    // use identical droput
    listCellDropoutRates.clear();
    listCellDropoutRates.resize( numTotCells );
    for(int i=0; i<numTotCells; ++i)
    {
        listCellDropoutRates[i] = rateDropout;
    }
    // add one more
    listCellDropoutRates.push_back(rateDropout);
    
    // now simulate genotype probability directly from genotypes
    // remember where dropout occurs: <site, cell>
    map<pair<int,int>, pair<int,int> > listDropPositions;
    
    for( int i=0; i<(int)listGenoSimCopiesPostDoubles.size(); ++i )
    {
        vector<double> vecGenoProb;
        vector<int> vecGenoSim, vecGenoNoise;
        
        // simulate genotypes
        map<int,pair<int,int> > &mapCellGenosUse = listGenoSimCopiesPostDoubles[i];
        for(map<int,pair<int,int> > :: iterator it = mapCellGenosUse.begin(); it != mapCellGenosUse.end(); ++it)
        {
//cout << "Site: " << i << ", cell : " << it->first << " has genotype: " << it->second.first << "," << it->second.second << endl;
            pair<int,int> genoNoise;
            //bool fDrop = false;
            double prob =  SimulateNoiseForGenoSingleGeno(it->first-1, i, it->second, genoNoise, listDropPositions);
            //if(fDrop)
            //{
//cout << "Site: " << i << ", cell: " << it->first-1 << " has a droput\n";
            //    pair<int,int> pp(i, it->first-1);
            //    listDropPositions.insert(pp);
            //}
            //cout << "Cell: " << it->first << ": " << cntReads.first << ", " << cntReads.second << endl;
            int alleleSum =  genoNoise.second;
            if( alleleSum >= 2 )
            {
                //cout << "Warning: Only binary genotype is supported right now\n";
                alleleSum = 1;
            }
            double pp = 0.0;
            int hap = 0;
            if( alleleSum == 0 )
            {
                //cout << prob << " ";
                //pp = prob;
                
                // YW: 0 is assumed to be good
                pp = CalcProb0(it->first-1);
                
                hap = 0;
            }
            else
            {
                //cout << 1.0 - prob << " ";
                //pp = 1.0 - prob;
                //
                pp = 1.0 - CalcProb1();
                
                hap = 1;
            }
            
#if 0
            // factor into boosting effect
            bool fAllele1orig = it->second.second >= 1;
            if( fAllele1orig )
            {
                // for now only deal with allele 1 (and droput)
                if( alleleSum == 0 )
                {
                    pp = CalcProb0Boost(pp);
                }
                else
                {
                    pp = CalcProb1Boost(pp);
                }
            }
#endif
            
            vecGenoProb.push_back(pp);
            vecGenoNoise.push_back(hap);
            
            int hapSim = 0;
            if( it->second.second >= 1 )
            {
                hapSim = 1;
            }
            vecGenoSim.push_back(hapSim);
        }
        //cout << endl;
        listGenoProbs.push_back(vecGenoProb);
        listGenoNoise.push_back(vecGenoNoise);
        listGenoSim.push_back( vecGenoSim );
    }
    
    // simulate reads from noiay genotypes
    for( int i=0; i<(int)listGenoNoise.size(); ++i )
    {
//cout << "^^^Site " << i << endl;
        vector<pair<int,int> > vecReadCntsSite;
        vector<vector<double> > vecReadBasedProb;
        
        // simulate reads from noisy genotypes
        //int cntAllele0 = CalcAllele0FreqForGenos(listGenoNoise[i]);
        //int cntTot = 2*listGenoNoise[i].size();
        //if( fBinary == false)
        //{
        //    cntTot *= 2;
        //}
        // use pseudo count=1
        //double allele0Freq = ((double)cntAllele0)/(cntTot+1);
//cout << "cntAllele0: " << cntAllele0 << ", cntTot: " << cntTot << ", freq: " << allele0Freq << endl;
        for(int j=0; j<(int)listGenoNoise[i].size(); ++j)
        {
            pair<int,int> readCnt;
            //SampleReadsForGenotype(j, i, listGenoNoise[i][j], listDropPositions, readCnt);
            SampleReadsForGenotype(j, i, listGenoSim[i][j], listDropPositions, readCnt);
            vecReadCntsSite.push_back(readCnt);
//cout << "[" << i << "," << j << "]: prob=" << probGenoReadsOut[0] << ", readCnt: " << readCnt.first << "," << readCnt.second << endl;
        }
        
        // now calculate the probability
        double allele0Freq = 0.0, aveReadDepthSite = 0.0, stdReadDepthSite = 0.0;
        AnalyzeReadsAtSites( vecReadCntsSite, allele0Freq, aveReadDepthSite, stdReadDepthSite  );
        
        // also impute dropouts
        vector<int> listInfDropouts;
        InfDroputsByReadsAtSite(vecReadCntsSite, listInfDropouts);
        // TBDDDDDDDDDDD
        
        for(int j=0; j<(int)vecReadCntsSite.size(); ++j)
        {
            // simulate prob
            //cout << "Prob at site: " << i << " cell " << j << endl;
            pair<int,int> readCnt = vecReadCntsSite[j];
            vector<double> probGenoReads;
//cout << "Site: " << i << ", cell: " << j << endl;
            CalcGenotypeProbOfReadsNew( j, aveReadDepthSite, stdReadDepthSite, readCnt, allele0Freq, listInfDropouts[j], probGenoReads );
            
            //CalcGenotypeProbOfReads( readCnt, allele0Freq, probGenoReads );
            vector<double> probGenoReadsOut;
            ConvGenoProbToOututProb(probGenoReads, fBinary, probGenoReadsOut);
            
            vecReadBasedProb.push_back(probGenoReadsOut);
            //cout << "[" << i << "," << j << "]: prob=" << probGenoReadsOut[0] << ", readCnt: " << readCnt.first << "," << readCnt.second << endl;
        }
        
        
        listReadCountsAtSites.push_back(vecReadCntsSite);
        listReadBasedProb.push_back(vecReadBasedProb);
    }
    
    // YW: Oct 10, 2025. A new way to simulate reads with dropouts and copy number changes
    // simulate reads from noiay genotypes
    for( int i=0; i<(int)listGenoSimCopiesPreDoubles.size(); ++i )
    {
//cout << "^^^Site " << i << endl;
        vector<pair<int,int> > vecReadCntsSite;
        
        for(map<int,pair<int,int> > :: iterator it = listGenoSimCopiesPreDoubles[i].begin(); it != listGenoSimCopiesPreDoubles[i].end(); ++it)
        {
//cout << "Cell " << it->first << " (" << it->second.first << "," << it->second.second << ") \n";
            pair<int,int> readCnt;
            //SampleReadsForGenotype(j, i, listGenoNoise[i][j], listDropPositions, readCnt);
            SampleReadsForGenotype2(it->first, i, it->second.first, it->second.second, readCnt);
            vecReadCntsSite.push_back(readCnt);
//cout << "[" << i << "," << j << "]: prob=" << probGenoReadsOut[0] << ", readCnt: " << readCnt.first << "," << readCnt.second << endl;
        }
        listReadCountsAtSites2.push_back(vecReadCntsSite);
    }
    
    
    
    if( fracMissing > 0.0)
    {
        if( numDoublet > 0 )
        {
            cout << "ERROR: for now missing value is not supported with non-zero doublets\n";
            exit(1);
        }
        ChooseMissingPos( listGenoProbs[0].size(), listGenoProbs.size(), setMissingPos );
    }

#if 0
    // dump out TRUE genotypes
    if( numDoublet > 0 )
    {
        //
        cout << "True genotype (BEFORE doublet merging, if any)\n";
        for(int i=0; i<(int)listGenoSimPreDbl.size(); ++i)
        {
            for(int j=0; j<(int)listGenoSimPreDbl[i].size(); ++j)
            {
                cout << listGenoSimPreDbl[i][j] << " ";
            }
            cout << endl;
        }
    }
    cout << "True genotype (with doublets merged, if any)\n";
    for(int i=0; i<(int)listGenoSim.size(); ++i)
    {
        for(int j=0; j<(int)listGenoSim[i].size(); ++j)
        {
            cout << listGenoSim[i][j] << " ";
        }
        cout << endl;
    }
    cout << "Noisy genotypes\n";
    for(int i=0; i<(int)listGenoNoise.size(); ++i)
    {
        for(int j=0; j<(int)listGenoNoise[i].size(); ++j)
        {
            if( IsMissingAt( j, i, setMissingPos ) == false )
            {
                cout << listGenoNoise[i][j] ;
            }
            else
            {
                cout << "*";
            }
            cout << " ";
        }
        cout << endl;
    }
    // compare
    int numDiff = 0;
    for(int i=0; i<(int)listGenoSim.size(); ++i)
    {
        for(int j=0; j<(int)listGenoSim[i].size(); ++j)
        {
            if( listGenoSim[i][j] != listGenoNoise[i][j]  )
            {
                cout << "Difference at: [ " << i+1 << " " << j+1 << " ]: changed to " << listGenoNoise[i][j] << endl;
                ++numDiff;
            }
            else if( IsMissingAt(j, i, setMissingPos) == true )
            {
                cout << "Gentoype missing at: [ " << i+1 << " " << j+1 << " ]\n";
            }
        }
    }
    cout << "Percentage of difference is: " << (double)numDiff/( listGenoNoise.size() * listGenoNoise[0].size() ) << endl;
    
    // YW: assume precision is 5 for probability (to avoid 0.99999999)
    cout << "----------------------------------------------------------------------\n";
    cout << "Probablistic genotypes\n\n";
    cout << "HAPLOID " << listGenoProbs.size() << " " << listGenoProbs[0].size() << endl;
    cout.precision(10);
    for(int i=0; i<(int)listGenoProbs.size(); ++i)
    {
        for( int j=0; j<(int)listGenoProbs[i].size(); ++j)
        {
            if( IsMissingAt(j, i, setMissingPos) == false )
            {
                cout << listGenoProbs[i][j];
            }
            else
            {
                cout << GetMissingProb();
            }
            cout << " ";
        }
        cout << endl;
    }
    
    // YW: assume precision is 5 for probability (to avoid 0.99999999)
    cout << "----------------------------------------------------------------------\n";
    cout << "Read-based genotypes\n\n";
    if( fBinary )
    {
        cout << "HAPLOID ";
    }
    else
    {
        cout << "TERNARY ";
    }
    cout << listReadBasedProb.size() << " " << listGenoProbs[0].size() << endl;
    cout.precision(10);
    for(int i=0; i<(int)listReadBasedProb.size(); ++i)
    {
        for( int j=0; j<(int)listReadBasedProb[i].size(); ++j)
        {
            for(int k=0; k<(int)listReadBasedProb[i][j].size(); ++k)
            {
                if( IsMissingAt(j,i, setMissingPos) == false )
                {
                    cout << listReadBasedProb[i][j][k] ;
                }
                else
                {
                    cout << GetMissingProb();
                }
                cout << " ";
            }
        }
        cout << endl;
    }
    vector<vector<int> > listCalledGEnoFromReadBasedProb;
    ConvReadBasedProbToGeno( listReadBasedProb, listCalledGEnoFromReadBasedProb );
    cout << "\nCalled genotypes from read-based genotypes\n";
    for(int i=0; i<(int)listCalledGEnoFromReadBasedProb.size(); ++i)
    {
        for( int j=0; j<(int)listCalledGEnoFromReadBasedProb[i].size(); ++j)
        {
            if( IsMissingAt(j, i, setMissingPos) == false )
            {
                cout << listCalledGEnoFromReadBasedProb[i][j];
            }
            else
            {
                // by default output zero
                cout << "*";
            }
            cout << " ";
        }
        cout << endl;
    }

    cout << "----------------------------------------------------------------------\n";
    cout << "Read Counts for Each Site and Cell\n";
    for (int i = 0; i < (int)listReadCountsAtSites.size(); ++i) {
        cout << "Site " << i + 1 << ":\n";
        for (int j = 0; j < (int)listReadCountsAtSites[i].size(); ++j) {
            cout << "  Cell " << j + 1 << ": Read Count (0,1) = (" 
                 << listReadCountsAtSites[i][j].first << ", " 
                 << listReadCountsAtSites[i][j].second << ")\n";
        }
    }
    
#endif
    
    cout << "----------------------------------------------------------------------\n";
    cout << "Read Counts for Each Site and Cell\n";
    for (int i = 0; i < (int)listReadCountsAtSites2.size(); ++i) {
        cout << "Site " << i + 1 << ":\n";
        for (int j = 0; j < (int)listReadCountsAtSites2[i].size(); ++j) {
            cout << "  Cell " << j + 1 << ": Read Count (0,1) = ("
                 << listReadCountsAtSites2[i][j].first << ", "
                 << listReadCountsAtSites2[i][j].second << ")\n";
        }
    }

    cout << "----------------------------------------------------------------------\n";
    cout << "Copy Numbers for Each Site and Cell\n";
    for (int i = 0; i < (int)listGenoSimCopiesPostDoubles.size(); ++i) {
        cout << "Site " << i + 1 << ":\n";
        for (const auto &entry : listGenoSimCopiesPostDoubles[i]) {
            int cell = entry.first;
            int copyNumber = entry.second.first + entry.second.second;
            cout << "  Cell " << cell << ": Copy Number = " << copyNumber << "\n";
        }
    }



}


//******************************************************************************************************
// Usage: ./xxxx tree-file numsites error-rate dropout-rate doublet-rate rate-cnv prob0 prob1
// tree-file: contain a single newick format

const char *VERSION = "scsim ver 1.0, August 27, 2018";

//
int main(int argc, char **argv)
{
    const int DEF_NUM_SITES = 11;
    int numSites = DEF_NUM_SITES;
    int numVariantsPerSite = 1;
    
	if(  argc < MIN_NUM_PARAMS || argc > MAX_NUM_PARAMS )
	{
		cout << "Arguments: wrong number. Usage: ./scsim tree-file numsites numVariantsPerSite error-rate dropout-rate doublet cnv-rate-inc cnv-rate-dec rec-rate binary droput-cell-var ave-read-depth std-read-depth rnd-seed fbetabinomial\n";
        cout << "\t tree-file: in Newick format w/ branch length\n";
        cout << "\t numsites: number of somatic single nucleiotide variants\n";
        cout << "\t numVariantsPerSite: number of variants per site";
        cout << "\t error-rate: error rate of genotypes\n";
        cout << "\t dropout-rate: rate of genotype dropout\n";
        cout << "\t doublet: number of doublet error (the number of output genotypes will subtract this number)\n";
        cout << "\t cnv-rate-inc: rate of copy number increase\n";
        cout << "\t cnv-rate-dec: rate of copy number decrease\n";
        cout << "\t rec-rate: rate of recurrent mutations\n";
        //cout << "\t prob1: confidence of genotype 1\n";
        //cout << "\t min-allele-freq: minimum mutant allele frequency (between 0 to 1, optional)\n";
        cout << "\t binary: 0,  ternary: 1 (optional)\n";
        cout << "\t droput cell variation (optional)\n";
        cout << "\t average read depth (optional)\n";
        cout << "\t read depth standard deviation (optional)\n";
        cout << "\t missing value (optional)\n";
        //cout << "\t cell copy bias standard deviation (optional)\n";
        cout << "\t rnd-seed: random seed (must be integer, optional)\n";
        cout << "\t fbetabinomial: turn on beta-binomial for read counts if set to 1 (0: use default normal distribution), optional\n";
        //cout << "\t boost-fac: boosting to make dropout more visible (at least 1.0, optional)\n";
		return -1;
	}
    
    sscanf( argv[2], "%d", &numSites);
    cout << "Number of sites: " << numSites << endl;
    
    sscanf( argv[3], "%d", &numVariantsPerSite);
    cout << "Number of SNVs per site: " << numVariantsPerSite << endl;
    
    float rateErrf;
    sscanf( argv[4], "%f", &rateErrf );
    rateErr = rateErrf;
    cout << "Rate of genotype error: " << rateErr << endl;
    
    float rateDropf;
    sscanf( argv[5], "%f", &rateDropf );
    rateDropout = rateDropf;
    cout << "Rate of dropout: " << rateDropout << endl;
    
    // assume the last several cells in each data are doublets since they are purely random anyway
    sscanf( argv[6], "%d", &numDoublet );
    cout << "Number of doublets: " << numDoublet << endl;
    
    float rateCNVIncf;
    sscanf( argv[7], "%f", &rateCNVIncf );
    rateCNVInc = rateCNVIncf;
    cout << "Rate of CNV increase: " << rateCNVInc << endl;
    
    float rateCNVDecf;
    sscanf( argv[8], "%f", &rateCNVDecf );
    rateCNVDec = rateCNVDecf;
    cout << "Rate of CNV decrease: " << rateCNVDec << endl;

    //
    float fracRecMf;
    sscanf( argv[9], "%f", &fracRecMf );
    fracRecurrentMut = fracRecMf;
    cout << "Fraction of recurrent mutations: " << fracRecMf << endl;
    
    //float prob1f;
    //sscanf( argv[8], "%f", &prob1f );
    //prob1Const = prob1f;
    //cout << "Confidence of genotype 1: " << prob1Const << endl;
    
    if( argc >= 11 )
    {
        // flag = 0: binary, 1: ternary
        int flagBin;
        sscanf(argv[10], "%d", &flagBin);
        if( flagBin == 0 )
        {
            fBinary = true;
        }
        else
        {
            cout << "Outputing ternary genotypes. Note: use read-based genotypes in this cases\n";
            fBinary = false;
        }
    }
    
    if( argc >= 12 )
    {
        float stdCellDropf = 0.0;
        sscanf(argv[11], "%f", &stdCellDropf);
        rateDropoutstd = stdCellDropf;
        cout << "Cell droput rate variation: " << rateDropoutstd << endl;
    }
    
    if( argc >= 13 )
    {
        float aveReadDepthf;
        sscanf(argv[12], "%f", &aveReadDepthf);
        aveReadDepth = aveReadDepthf;
        cout << "Average read depth: " << aveReadDepth << endl;
#if 0
        float minAlleFreqUse;
        sscanf(argv[10], "%f", &minAlleFreqUse);
        if( minAlleFreqUse < 0.0 || minAlleFreqUse >= 1.0 )
        {
            cout << "Minor allele frequency must be between 0.0 and 1.0\n";
            exit(1);
        }
        minAlleFreq = minAlleFreqUse;
        cout << "Minor mutant allele frequency: " << minAlleFreq << endl;

        // the previous must also be present
        float facBoostUse;
        sscanf(argv[10], "%f", &facBoostUse);
        if( facBoostUse < 1.0 )
        {
            cout << "Boosting factor must be at least 1.0\n";
            exit(1);
        }
        boostFac = facBoostUse;
        cout << "Boosting factor: " << boostFac << endl;
#endif
    }
    if( argc >= 14 )
    {
        float stdReadDepthf;
        sscanf(argv[13], "%f", &stdReadDepthf);
        stdReadDepth = stdReadDepthf;
        cout << "STD of read depth: " << stdReadDepth << endl;
    }

    if( argc >= 15 )
    {
        float fracMissingf;
        sscanf(argv[14], "%f", &fracMissingf);
        fracMissing = fracMissingf;
        cout << "Missing value rate: " << fracMissing << endl;
    }
    
    if( argc >= 16 )
    {
        seedRndUserDef = false;
        sscanf( argv[15], "%d", &seedRndUser);
        cout << "Random seed: " << seedRndUser << endl;
    }
    
    if( argc >= 17 )
    {
        int fBeta = 0;
        sscanf( argv[16], "%d", &fBeta);
        if(fBeta == 1 )
        {
            fBetaBinomial = true;
            cout << "Turn on beta-binomial mode for allele counts distribution.\n";
        }
    }
    
    // read the file
    
    CellTreeNode *pRoot;
    ifstream inFile(argv[1]);
    while( !inFile.eof() )
    {
        char buf[102400];
        inFile.getline(buf, 102400);
        //inFile.close();
        string strTreeFile = buf;

        if( strTreeFile.length() == 0 )
        {
            break;
        }

        // now process the tree and output the
        pRoot = ProcTreeStr( strTreeFile );
cout << "Read tree: " << pRoot->GetNewick() << endl;
        //delete pRoot;

        // right now only process one tree
        break;
    }

    inFile.close();
    

    // start simulation
    SimulateIndependentSitesNoisyGeno(pRoot, numSites, numVariantsPerSite);
    

    delete pRoot;

    return 0;
}

